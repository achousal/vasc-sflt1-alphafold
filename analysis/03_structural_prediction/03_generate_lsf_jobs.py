"""03_generate_lsf_jobs.py -- Generate orchestrated LSF submission for AlphaFold Multimer.

Produces per-candidate LSF job scripts with sentinel-based completion tracking,
a shared wrapper script, a job manifest, and a barrier-orchestrator that monitors
all jobs in parallel with fail-fast error detection.

Target cluster: Minerva (Mount Sinai), LSF scheduler.

Architecture adapted from CeliacRisks HPC orchestration pattern:
  - Job manifest (JSON) with metadata for all jobs
  - Wrapper script with EXIT trap sentinel writes
  - Orchestrator LSF job: submits all AF jobs, polls sentinels, detects failures
  - Chunked batch submission to prevent scheduler flooding
  - Structured logging (orchestrator_state.jsonl)
"""

import base64
import json
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Orchestrator tuning
BATCH_CHUNK_SIZE = 10
POLL_INTERVAL_SECONDS = 60
AF_JOB_TIMEOUT_SECONDS = 36 * 3600  # 36h for all AF jobs to finish
ORCHESTRATOR_WALLTIME = "48:00"
ORCHESTRATOR_MEM = 2000  # MB

# AF job resources (per job)
AF_QUEUE = "gpu"
AF_CORES = 4
AF_MEM = 16000  # MB per core
AF_WALLTIME = "24:00"
AF_GPU = 1

# HPC defaults
DEFAULT_PROJECT_ACCOUNT = "acc_vascbrain"
DEFAULT_GPU_TYPE = ""


def _safe_name(target: str) -> str:
    """Sanitize a target name for filesystem and shell use."""
    return target.replace("/", "-").replace(" ", "_").replace(":", "_")


def _encode_b64(command: str) -> str:
    """Base64-encode a command for safe shell transport."""
    return base64.b64encode(command.encode("utf-8")).decode("ascii")


# ---------------------------------------------------------------------------
# Component generators
# ---------------------------------------------------------------------------


def _build_af_command(
    fasta_path: str,
    output_dir: str,
    hpc_root: str = "",
) -> str:
    """Build the AlphaFold multimer CLI command (singularity container).

    Parameters
    ----------
    fasta_path : str
        Path to the two-chain FASTA file (relative to project root).
    output_dir : str
        Path for AF output (relative to project root).
    hpc_root : str
        If set, cd to this directory before running. Makes relative paths
        resolve correctly on HPC.
    """
    lines = []
    if hpc_root:
        lines.append(f'cd "{hpc_root}"')
    lines.extend([
        "module purge",
        "module load alphafold/2.3.2-singularity",
        "",
        f'FASTA_INPUT="{fasta_path}"',
        f'OUTPUT_DIR="{output_dir}"',
        "",
        'mkdir -p "$OUTPUT_DIR"',
        "",
        "# Resolve to absolute paths for singularity bind mounts",
        'WORK_DIR="$(pwd)"',
        'if [[ "$FASTA_INPUT" != /* ]]; then',
        '    FASTA_ABS="${WORK_DIR}/${FASTA_INPUT}"',
        "else",
        '    FASTA_ABS="$FASTA_INPUT"',
        "fi",
        'if [[ "$OUTPUT_DIR" != /* ]]; then',
        '    OUTPUT_ABS="${WORK_DIR}/${OUTPUT_DIR}"',
        "else",
        '    OUTPUT_ABS="$OUTPUT_DIR"',
        "fi",
        "",
        'echo "[$(date)] Starting AlphaFold Multimer (singularity)"',
        'echo "FASTA: $FASTA_ABS"',
        'echo "Output: $OUTPUT_ABS"',
        'echo "Container: $AF2IMAGE"',
        'echo "Data: $AF2DATA"',
        "",
        "singularity run --nv \\",
        '    --bind "${AF2DATA}":/data \\',
        "    --bind /sc/arion:/sc/arion \\",
        '    "$AF2IMAGE" \\',
        '    --fasta_paths="$FASTA_ABS" \\',
        '    --output_dir="$OUTPUT_ABS" \\',
        "    --data_dir=/data \\",
        "    --uniref90_database_path=/data/uniref90/uniref90.fasta \\",
        "    --mgnify_database_path=/data/mgnify/mgy_clusters_2022_05.fa \\",
        "    --template_mmcif_dir=/data/pdb_mmcif/mmcif_files \\",
        "    --obsolete_pdbs_path=/data/pdb_mmcif/obsolete.dat \\",
        "    --pdb_seqres_database_path=/data/pdb_seqres/pdb_seqres.txt \\",
        "    --uniprot_database_path=/data/uniprot/uniprot.fasta \\",
        "    --uniref30_database_path=/data/uniref30/UniRef30_2021_03 \\",
        "    --bfd_database_path=/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\",
        "    --model_preset=multimer \\",
        "    --db_preset=full_dbs \\",
        "    --max_template_date=2024-01-01 \\",
        "    --num_multimer_predictions_per_model=5 \\",
        "    --use_gpu_relax",
        "",
        'echo "[$(date)] AlphaFold Multimer complete"',
    ])
    return "\n".join(lines)


def _build_wrapper_script() -> str:
    """Build the shared wrapper script.

    Every AF job runs through this wrapper, which:
    1. Validates required environment variables
    2. Decodes the base64-encoded command
    3. Executes it with error trapping
    4. Writes job name to sentinel completed.log on EXIT (success or failure)

    The EXIT trap ensures the sentinel is always written, so the orchestrator
    can detect both successes and failures. check_upstream_failures() then
    distinguishes success from failure via bjobs/bhist.
    """
    return """#!/bin/bash
set -euo pipefail

if [ -z "${AF_JOB_COMMAND_B64:-}" ]; then
    echo "[$(date '+%F %T')] FATAL: AF_JOB_COMMAND_B64 not set"
    exit 1
fi
if [ -z "${AF_JOB_NAME:-}" ]; then
    echo "[$(date '+%F %T')] FATAL: AF_JOB_NAME not set"
    exit 1
fi
if [ -z "${AF_SENTINEL_DIR:-}" ]; then
    echo "[$(date '+%F %T')] FATAL: AF_SENTINEL_DIR not set"
    exit 1
fi

# Write sentinel on exit (success or failure) so orchestrator always sees completion.
# check_upstream_failures() distinguishes EXIT/TERM from success via bjobs/bhist.
_af_rc=0
trap 'echo "${AF_JOB_NAME}" >> "$AF_SENTINEL_DIR/completed.log"' EXIT

COMMAND=$(python3 - "$AF_JOB_COMMAND_B64" <<'PY'
import base64
import sys
print(base64.b64decode(sys.argv[1]).decode("utf-8"), end="")
PY
)

eval "$COMMAND" || _af_rc=$?
exit "$_af_rc"
"""


def _build_job_manifest(
    candidates: pd.DataFrame,
    fasta_dir: Path,
    results_base_dir: Path,
    hpc_root: str = "",
) -> dict:
    """Build JSON manifest mapping job_key -> job metadata.

    Parameters
    ----------
    candidates : pd.DataFrame
        Candidates with target, uniprot columns.
    fasta_dir : Path
        Directory containing two-chain FASTA files.
    results_base_dir : Path
        Base directory where AF outputs go (e.g. results/).
    hpc_root : str
        If set, injected into AF commands as working directory.

    Returns
    -------
    dict
        Mapping of job_key -> {job_name, target, uniprot, queue, cores,
        mem_per_core, walltime, fasta_path, output_dir, command_b64}.
    """
    manifest = {}
    for _, row in candidates.iterrows():
        target = row["target"]
        uniprot = row["uniprot"]
        sname = _safe_name(target)

        fasta_path = fasta_dir / f"sflt1_vs_{sname}.fasta"
        if not fasta_path.exists():
            logger.warning("No FASTA for %s, skipping manifest entry", target)
            continue

        output_dir = results_base_dir / sname
        command = _build_af_command(
            str(fasta_path), str(output_dir), hpc_root=hpc_root
        )

        job_key = f"af2_{sname}"
        manifest[job_key] = {
            "job_name": f"af2_{sname}",
            "target": target,
            "uniprot": uniprot,
            "queue": AF_QUEUE,
            "cores": AF_CORES,
            "mem_per_core": AF_MEM,
            "gpu": AF_GPU,
            "walltime": AF_WALLTIME,
            "fasta_path": str(fasta_path),
            "output_dir": str(output_dir),
            "command_b64": _encode_b64(command),
        }

    return manifest


def _build_orchestrator_functions() -> str:
    """Build shared bash functions for the orchestrator.

    Contains:
    - manifest_job_tsv(): read manifest and return TSV for a job key
    - check_upstream_failures(): fail-fast on EXIT/TERM via bjobs/bhist
    - barrier_wait(): poll sentinels with timeout
    - submit_and_track(): submit one job via bsub, capture ID
    - submit_batch(): submit jobs in chunks with pauses
    """
    return r"""
# ---------------------------------------------------------------------------
# Orchestrator bash functions
# ---------------------------------------------------------------------------

manifest_job_tsv() {
    local job_key="$1"
    python3 - "$MANIFEST_PATH" "$job_key" <<'PY'
import json
import sys

manifest_path, job_key = sys.argv[1], sys.argv[2]
with open(manifest_path, encoding="utf-8") as f:
    jobs = json.load(f)

job = jobs.get(job_key)
if job is None:
    raise SystemExit(1)

fields = [
    job["job_name"],
    str(job["queue"]),
    str(job["cores"]),
    str(job["mem_per_core"]),
    str(job["gpu"]),
    str(job["walltime"]),
    job["command_b64"],
]
print("\t".join(fields), end="")
PY
}

check_upstream_failures() {
    local -a job_ids=("$@")
    local jid
    for jid in "${job_ids[@]}"; do
        [ -z "$jid" ] && continue

        local raw
        raw=$(bjobs -noheader -o "stat" "$jid" 2>&1 || true)
        local stat
        stat=$(echo "$raw" | awk 'NF {print $1; exit}')

        if [ "$stat" = "EXIT" ] || [ "$stat" = "TERM" ]; then
            local jname
            jname=$(bjobs -noheader -o "job_name" "$jid" 2>/dev/null | awk 'NF {print $1; exit}')
            echo "[$(date '+%F %T')] FATAL: upstream job $jid ($jname) $stat (bjobs)"
            printf '{"event":"job_failed","job_id":"%s","job_name":"%s","status":"%s","ts":"%s"}\n' \
                "$jid" "$jname" "$stat" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
            exit 1
        fi

        if [ -z "$stat" ] || echo "$raw" | grep -qi "not found"; then
            local hist_exit
            hist_exit=$(bhist -l "$jid" 2>/dev/null | awk '/Completed <exit>|Exited with exit code/ {print; exit}' || true)
            if [ -n "$hist_exit" ]; then
                echo "[$(date '+%F %T')] FATAL: upstream job $jid EXIT (bhist): $hist_exit"
                printf '{"event":"job_failed","job_id":"%s","status":"EXIT","ts":"%s"}\n' \
                    "$jid" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
                exit 1
            fi
            local hist_term
            hist_term=$(bhist -l "$jid" 2>/dev/null | awk '/TERM/ {print; exit}' || true)
            if [ -n "$hist_term" ]; then
                echo "[$(date '+%F %T')] FATAL: upstream job $jid TERM (bhist): $hist_term"
                printf '{"event":"job_failed","job_id":"%s","status":"TERM","ts":"%s"}\n' \
                    "$jid" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
                exit 1
            fi
        fi
    done
}

barrier_wait() {
    local label="$1"; shift
    local timeout="$1"; shift
    local poll="$1"; shift
    local -a job_names=("$@")
    local total=${#job_names[@]}
    local elapsed=0

    printf '{"event":"barrier_start","stage":"%s","n_jobs":%d,"ts":"%s"}\n' \
        "$label" "$total" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
    echo "[$(date '+%F %T')] Waiting for $label ($total jobs, timeout=${timeout}s)..."

    if [ "$total" -eq 0 ]; then
        printf '{"event":"barrier_done","stage":"%s","ts":"%s"}\n' \
            "$label" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
        return 0
    fi

    while true; do
        if [ ${#UPSTREAM_IDS[@]} -gt 0 ]; then
            check_upstream_failures "${UPSTREAM_IDS[@]}"
        fi

        local missing=0
        local done_count=0
        local name
        for name in "${job_names[@]}"; do
            if grep -qx "${name}" "$SENTINEL_DIR/completed.log" 2>/dev/null; then
                done_count=$((done_count + 1))
            else
                missing=$((missing + 1))
            fi
        done

        if [ "$missing" -eq 0 ]; then
            echo "[$(date '+%F %T')] $label complete ($total/$total)."
            printf '{"event":"barrier_done","stage":"%s","ts":"%s"}\n' \
                "$label" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
            return 0
        fi

        if [ "$elapsed" -ge "$timeout" ]; then
            echo "[$(date '+%F %T')] TIMEOUT: $label after ${timeout}s ($done_count/$total done)"
            echo "[$(date '+%F %T')] Missing jobs:"
            for name in "${job_names[@]}"; do
                grep -qx "${name}" "$SENTINEL_DIR/completed.log" 2>/dev/null || echo "  $name"
            done
            printf '{"event":"barrier_timeout","stage":"%s","done":%d,"total":%d,"ts":"%s"}\n' \
                "$label" "$done_count" "$total" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
            exit 1
        fi

        echo "[$(date '+%F %T')] Progress: $done_count/$total done, polling in ${poll}s..."
        sleep "$poll"
        elapsed=$((elapsed + poll))
    done
}

submit_and_track() {
    local job_key="$1"
    local label="$2"
    local id_file="$3"

    local job_tsv
    if ! job_tsv=$(manifest_job_tsv "$job_key"); then
        echo "[$(date '+%F %T')] FATAL: no manifest entry for $job_key"
        exit 1
    fi

    local job_name queue cores mem_per_core gpu walltime command_b64
    IFS=$'\t' read -r job_name queue cores mem_per_core gpu walltime command_b64 <<< "$job_tsv"

    local job_script
    local bsub_directive="#BSUB"
    local account_directive=""
    if [ -n "$PROJECT_ACCOUNT" ]; then
        account_directive="${bsub_directive} -P $PROJECT_ACCOUNT"
    fi
    local gpu_type_directive=""
    if [ -n "$GPU_TYPE" ]; then
        gpu_type_directive="${bsub_directive} -R $GPU_TYPE"
    fi
    job_script=$(cat <<EOF
#!/bin/bash
${bsub_directive} -L /bin/bash
${bsub_directive} -J $job_name
${account_directive:+$account_directive}
${bsub_directive} -q $queue
${bsub_directive} -n $cores
${bsub_directive} -R "rusage[mem=$mem_per_core] span[hosts=1]"
${bsub_directive} -gpu "num=$gpu"
${gpu_type_directive:+$gpu_type_directive}
${bsub_directive} -W $walltime
${bsub_directive} -o $LOG_DIR/${job_name}_%J.out
${bsub_directive} -e $LOG_DIR/${job_name}_%J.err

set -euo pipefail
export AF_JOB_COMMAND_B64="$command_b64"
export AF_JOB_NAME="$job_name"
export AF_SENTINEL_DIR="$SENTINEL_DIR"
"$WRAPPER_SCRIPT"
EOF
)

    local output
    output=$(echo "$job_script" | bsub 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "[$(date '+%F %T')] FATAL: bsub failed for $label (rc=$rc): $output"
        printf '{"event":"submit_failed","job":"%s","rc":%d,"ts":"%s"}\n' \
            "$label" "$rc" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
        exit 1
    fi

    local job_id
    job_id=$(echo "$output" | sed -n 's/.*Job <\([0-9]*\)>.*/\1/p')
    if [ -z "$job_id" ]; then
        echo "[$(date '+%F %T')] FATAL: cannot parse job ID for $label: $output"
        exit 1
    fi

    echo "[$(date '+%F %T')] Submitted $label: Job $job_id"
    printf '{"event":"submitted","job":"%s","job_id":"%s","ts":"%s"}\n' \
        "$label" "$job_id" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
    echo "$job_id" >> "$id_file"
}

submit_batch() {
    local chunk_size="$1"; shift
    local id_file="$1"; shift
    local -a job_keys=("$@")
    local i
    for ((i=0; i<${#job_keys[@]}; i++)); do
        submit_and_track "${job_keys[$i]}" "${job_keys[$i]}" "$id_file"
        if (( (i + 1) % chunk_size == 0 && i + 1 < ${#job_keys[@]} )); then
            echo "[$(date '+%F %T')] Submitted $((i+1))/${#job_keys[@]}, pausing..."
            sleep 2
        fi
    done
}
"""


def _build_orchestrator_script(
    manifest: dict,
    manifest_path: Path,
    wrapper_path: Path,
    sentinel_dir: Path,
    log_dir: Path,
    state_file: Path,
    project_account: str = "",
    gpu_type: str = "",
) -> str:
    """Build the complete orchestrator LSF script.

    The orchestrator runs as a lightweight LSF job that:
    1. Submits all AF jobs in batches (prevent scheduler flooding)
    2. Polls sentinel files for completion every POLL_INTERVAL_SECONDS
    3. Detects job failures via bjobs/bhist (fail-fast)
    4. Logs state transitions to orchestrator_state.jsonl
    5. Exits successfully when all jobs complete
    """
    job_keys = sorted(manifest.keys())
    job_names = [manifest[k]["job_name"] for k in job_keys]

    # Build bash arrays
    keys_array = "JOB_KEYS=(\n" + "\n".join(f'  "{k}"' for k in job_keys) + "\n)"
    names_array = "JOB_NAMES=(\n" + "\n".join(f'  "{n}"' for n in job_names) + "\n)"

    account_line = f"\n#BSUB -P {project_account}" if project_account else ""
    gpu_type_export = f'\nGPU_TYPE="{gpu_type}"' if gpu_type else '\nGPU_TYPE=""'

    return f"""#!/bin/bash
#BSUB -L /bin/bash
#BSUB -J af2_orchestrator{account_line}
#BSUB -q premium
#BSUB -n 1
#BSUB -R "rusage[mem={ORCHESTRATOR_MEM}] span[hosts=1]"
#BSUB -W {ORCHESTRATOR_WALLTIME}
#BSUB -oo {log_dir}/orchestrator_%J.log
#BSUB -eo {log_dir}/orchestrator_%J.log

# AlphaFold Multimer Orchestrator
# Submits {len(manifest)} jobs, monitors via sentinels, fail-fast on errors.
# Generated: {datetime.now().isoformat(timespec="seconds")}

set -euo pipefail

MANIFEST_PATH="{manifest_path.resolve()}"
WRAPPER_SCRIPT="{wrapper_path.resolve()}"
SENTINEL_DIR="{sentinel_dir.resolve()}"
STATE_FILE="{state_file.resolve()}"
LOG_DIR="{log_dir.resolve()}"
POLL_INTERVAL={POLL_INTERVAL_SECONDS}
BATCH_CHUNK={BATCH_CHUNK_SIZE}
AF_TIMEOUT={AF_JOB_TIMEOUT_SECONDS}
EXPECTED_JOBS={len(manifest)}
PROJECT_ACCOUNT="{project_account}"{gpu_type_export}

{_build_orchestrator_functions()}

# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

mkdir -p "$SENTINEL_DIR" "$LOG_DIR"
touch "$STATE_FILE"
touch "$SENTINEL_DIR/completed.log"

echo "[$(date '+%F %T')] Orchestrator started: $EXPECTED_JOBS AlphaFold jobs"
printf '{{"event":"orchestrator_start","n_jobs":%d,"ts":"%s"}}\\n' \
    "$EXPECTED_JOBS" "$(date -u '+%FT%TZ')" >> "$STATE_FILE"

# Job arrays
{keys_array}
{names_array}

# Submit all AF jobs in batches
IDS_FILE=$(mktemp "$SENTINEL_DIR/job_ids.XXXXXX")
submit_batch "$BATCH_CHUNK" "$IDS_FILE" "${{JOB_KEYS[@]}}"

mapfile -t UPSTREAM_IDS < "$IDS_FILE"
echo "[$(date '+%F %T')] All ${{#UPSTREAM_IDS[@]}} jobs submitted."

# Wait for all to complete
barrier_wait "alphafold" "$AF_TIMEOUT" "$POLL_INTERVAL" "${{JOB_NAMES[@]}}"

echo "[$(date '+%F %T')] All AlphaFold jobs complete."
printf '{{"event":"orchestrator_done","ts":"%s"}}\\n' "$(date -u '+%FT%TZ')" >> "$STATE_FILE"
echo ""
echo "Next steps:"
echo "  python analysis/03_structural_prediction/run_structural.py --step parse --af-output-dir results/"
echo "  python analysis/03_structural_prediction/run_structural.py --step plot"
"""


def _build_submit_orchestrator_script(
    orchestrator_path: Path,
) -> str:
    """Build the convenience script that submits the orchestrator job."""
    return f"""#!/bin/bash
# submit_orchestrator.sh -- Submit the AlphaFold orchestrator job
# Generated: {datetime.now().isoformat(timespec="seconds")}
#
# Usage: bash submit_orchestrator.sh
#   or:  bash submit_orchestrator.sh --dry-run
#
# The orchestrator job:
#   1. Submits all AlphaFold Multimer jobs in parallel (batches of {BATCH_CHUNK_SIZE})
#   2. Monitors completion via sentinel files
#   3. Detects failures via bjobs/bhist (fail-fast)
#   4. Logs state to orchestrator_state.jsonl
#
# Monitor: bjobs -w | grep af2_
# Logs:    tail -f logs/orchestrator_*.log
# State:   cat logs/sentinels/orchestrator_state.jsonl | python3 -m json.tool --no-ensure-ascii

set -euo pipefail

DRY_RUN=false
if [[ "${{1:-}}" == "--dry-run" ]]; then
    DRY_RUN=true
fi

ORCH_SCRIPT="{orchestrator_path}"

if [[ "$DRY_RUN" == "true" ]]; then
    echo "DRY RUN: would submit orchestrator:"
    echo "  bsub < $ORCH_SCRIPT"
    echo ""
    echo "Orchestrator will submit these AlphaFold jobs:"
    python3 - "{orchestrator_path.parent / 'manifest.json'}" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    m = json.load(f)
for k in sorted(m):
    print(f'  {{k}}: {{m[k]["target"]}} ({{m[k]["uniprot"]}})')
print(f'\\nTotal: {{len(m)}} jobs')
PY
    exit 0
fi

echo "Submitting AlphaFold orchestrator..."
bsub < "$ORCH_SCRIPT"
echo ""
echo "Monitor with:"
echo "  bjobs -w | grep af2_"
echo "  tail -f logs/orchestrator_*.log"
echo "  cat logs/sentinels/orchestrator_state.jsonl"
"""


# ---------------------------------------------------------------------------
# Per-job LSF scripts (standalone, for manual resubmission)
# ---------------------------------------------------------------------------


def _build_standalone_lsf_script(
    job_entry: dict,
    wrapper_path: Path,
    sentinel_dir: Path,
    log_dir: Path,
    timestamp: str,
    project_account: str = "",
    gpu_type: str = "",
) -> str:
    """Build a standalone LSF script for one AF job.

    These can be submitted individually for reruns:
        bsub < jobs/af2_VEGFA.lsf
    """
    job_name = job_entry["job_name"]
    account_line = f"\n#BSUB -P {project_account}" if project_account else ""
    gpu_type_line = f"\n#BSUB -R {gpu_type}" if gpu_type else ""
    return f"""#!/bin/bash
#BSUB -L /bin/bash
#BSUB -J {job_name}{account_line}
#BSUB -q {job_entry["queue"]}
#BSUB -n {job_entry["cores"]}
#BSUB -R "rusage[mem={job_entry["mem_per_core"]}] span[hosts=1]"
#BSUB -gpu "num={job_entry["gpu"]}"{gpu_type_line}
#BSUB -W {job_entry["walltime"]}
#BSUB -o {log_dir}/{job_name}_%J.out
#BSUB -e {log_dir}/{job_name}_%J.err

# AlphaFold Multimer: sFLT1 vs {job_entry["target"]}
# UniProt: {job_entry["uniprot"]}
# Generated: {timestamp}
#
# Standalone script for manual resubmission.
# For parallel submission of all jobs, use submit_all.sh instead.

set -euo pipefail
export AF_JOB_COMMAND_B64="{job_entry["command_b64"]}"
export AF_JOB_NAME="{job_name}"
export AF_SENTINEL_DIR="{sentinel_dir.resolve()}"
"{wrapper_path.resolve()}"
"""


def _build_submit_all_script(jobs_dir: Path) -> str:
    """Build a simple direct-submission script for all AF jobs.

    Loops over all af2_*.lsf files and submits them directly via bsub,
    with a short delay between submissions to avoid scheduler flooding.
    """
    return f"""#!/bin/bash
# submit_all.sh -- Submit all AlphaFold jobs directly (parallel)
# Generated: {datetime.now().isoformat(timespec="seconds")}
#
# Usage: bash submit_all.sh
#   or:  bash submit_all.sh --dry-run

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
DRY_RUN=false
if [[ "${{1:-}}" == "--dry-run" ]]; then
    DRY_RUN=true
fi

count=0
for lsf in "$SCRIPT_DIR"/af2_*.lsf; do
    [ -f "$lsf" ] || continue
    name=$(basename "$lsf" .lsf)
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] bsub < $lsf"
    else
        echo "Submitting $name..."
        bsub < "$lsf"
        sleep 2
    fi
    count=$((count + 1))
done

if [[ "$DRY_RUN" == "true" ]]; then
    echo ""
    echo "Would submit $count jobs. Remove --dry-run to submit."
else
    echo ""
    echo "Submitted $count jobs."
    echo "Monitor: bjobs -w | grep af2_"
fi
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_lsf_scripts(
    candidates_path: Path,
    fasta_dir: Path,
    jobs_dir: Path,
    project_account: str = "",
    gpu_type: str = "",
    hpc_root: str = "",
) -> list[Path]:
    """Generate orchestrated LSF infrastructure for all candidates.

    Produces:
      - manifest.json: job metadata for all candidates
      - wrapper.sh: shared command wrapper with sentinel writes
      - orchestrator.lsf: barrier-orchestrator LSF script
      - submit_orchestrator.sh: convenience submission script
      - af2_*.lsf: per-job standalone scripts (for manual resubmission)

    Parameters
    ----------
    candidates_path : Path
        Path to step03_candidates.csv.
    fasta_dir : Path
        Directory containing two-chain FASTA files.
    jobs_dir : Path
        Output directory for all generated scripts.
    project_account : str
        LSF project account (e.g., 'acc_vascbrain'). Added as #BSUB -P.
    gpu_type : str
        GPU resource constraint (e.g., 'a100'). Added as #BSUB -R.
    hpc_root : str
        Absolute path to project root on HPC. If set, AF commands cd here
        first and all paths resolve relative to it.

    Returns
    -------
    list[Path]
        Paths to generated per-job LSF scripts.
    """
    candidates = pd.read_csv(candidates_path)
    jobs_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().isoformat(timespec="seconds")

    # Use hpc_root for path resolution if provided
    if hpc_root:
        hpc_base = Path(hpc_root)
        log_dir = hpc_base / "logs"
        sentinel_dir = log_dir / "sentinels"
        state_file = sentinel_dir / "orchestrator_state.jsonl"
        results_base = Path("results")  # relative, resolved via cd in command
    else:
        log_dir = Path("logs")
        sentinel_dir = log_dir / "sentinels"
        state_file = sentinel_dir / "orchestrator_state.jsonl"
        results_base = Path("results")

    # 1. Build manifest
    manifest = _build_job_manifest(
        candidates, fasta_dir, results_base, hpc_root=hpc_root
    )
    manifest_path = jobs_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    logger.info("Generated: manifest.json (%d jobs)", len(manifest))

    # 2. Write wrapper script
    wrapper_path = jobs_dir / "wrapper.sh"
    wrapper_path.write_text(_build_wrapper_script())
    wrapper_path.chmod(0o750)
    logger.info("Generated: wrapper.sh")

    # 3. Build orchestrator script
    # For the orchestrator, use HPC paths for manifest/wrapper/sentinel
    if hpc_root:
        hpc_jobs_dir = hpc_base / jobs_dir.name
        orch_manifest_path = hpc_jobs_dir / "manifest.json"
        orch_wrapper_path = hpc_jobs_dir / "wrapper.sh"
    else:
        orch_manifest_path = manifest_path
        orch_wrapper_path = wrapper_path

    orchestrator_path = jobs_dir / "orchestrator.lsf"
    orchestrator_content = _build_orchestrator_script(
        manifest=manifest,
        manifest_path=orch_manifest_path,
        wrapper_path=orch_wrapper_path,
        sentinel_dir=sentinel_dir,
        log_dir=log_dir,
        state_file=state_file,
        project_account=project_account,
        gpu_type=gpu_type,
    )
    orchestrator_path.write_text(orchestrator_content)
    orchestrator_path.chmod(0o750)
    logger.info("Generated: orchestrator.lsf")

    # 4. Generate per-job standalone LSF scripts
    scripts = []
    for job_key, job_entry in manifest.items():
        script_content = _build_standalone_lsf_script(
            job_entry, orch_wrapper_path, sentinel_dir, log_dir, timestamp,
            project_account=project_account, gpu_type=gpu_type,
        )
        script_path = jobs_dir / f"{job_key}.lsf"
        script_path.write_text(script_content)
        scripts.append(script_path)
        logger.info("Generated: %s", script_path.name)

    # 5. Generate submission convenience script
    if hpc_root:
        orch_submit_path_ref = hpc_jobs_dir / "orchestrator.lsf"
    else:
        orch_submit_path_ref = orchestrator_path
    submit_path = jobs_dir / "submit_orchestrator.sh"
    submit_path.write_text(
        _build_submit_orchestrator_script(orch_submit_path_ref)
    )
    submit_path.chmod(0o755)
    logger.info("Generated: submit_orchestrator.sh")

    # 6. Generate submit_all.sh for direct parallel submission
    submit_all_path = jobs_dir / "submit_all.sh"
    submit_all_path.write_text(_build_submit_all_script(jobs_dir))
    submit_all_path.chmod(0o755)
    logger.info("Generated: submit_all.sh")

    logger.info(
        "Generated %d LSF scripts + orchestrator infrastructure", len(scripts)
    )
    return scripts
