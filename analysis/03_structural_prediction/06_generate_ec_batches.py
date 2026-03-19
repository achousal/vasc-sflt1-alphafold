#!/usr/bin/env python3
"""06_generate_ec_batches.py -- Generate corrected D1-D3 + D1-D6 + D1-D7 AF2 job batches.

Produces three independent job batches with construct-appropriate walltimes:
  Batch A: D1-D3 corrected (27-330, 304aa) -- all 24 original targets
  Batch B: D1-D6 (27-657, 631aa) -- top 5 hits + bottom 5 misses + NRP1
  Batch C: D1-D7 (27-747, 721aa) -- same 11 targets, deprioritized

Each batch gets its own fasta/, jobs/, and results/ subdirectories under
results/03_structural_prediction/{batch_name}/.

Walltime heuristics (A100, AF2 2.3.2, 25 models):
  - <800 total residues: 48h
  - 800-1200: 72h
  - 1200-1800: 96h
  - >1800: 144h (max practical on Minerva)

Usage:
    python 06_generate_ec_batches.py [--dry-run]
    python 06_generate_ec_batches.py --hpc-root /sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import importlib

_mod02 = importlib.import_module("02_fetch_sequences")
_mod03 = importlib.import_module("03_generate_lsf_jobs")

CONSTRUCTS = _mod02.CONSTRUCTS
CONSTRUCT_D1D3 = _mod02.CONSTRUCT_D1D3
CONSTRUCT_D1D6 = _mod02.CONSTRUCT_D1D6
CONSTRUCT_D1D7 = _mod02.CONSTRUCT_D1D7
fetch_all_sequences = _mod02.fetch_all_sequences

# ---------------------------------------------------------------------------
# Target lists
# ---------------------------------------------------------------------------

# All 24 original D1-D3 targets (from step03_candidates.csv)
ALL_24_TARGETS = [
    "VEGFA", "PCDH9", "NOE1", "SLIK4", "Amyloid-like protein 1",
    "MER", "APLP2.1", "STMN3", "Dtk", "PLXA4",
    "FSTL4", "DSCAM", "BAI1", "NRP1", "PLXA1",
    "NRX1B.1", "NRP2", "Calcineurin B a", "BASI", "Contactin-5",
    "SEMA3A", "SLIT2", "NGL1", "ROBO2",
]

# EC domain targets: top 5 + bottom 5 + NRP1
EC_TARGETS = [
    # Top 5 by D1-D3 ipTM
    "VEGFA",                    # 0.816 (positive control)
    "PCDH9",                    # 0.806
    "NOE1",                     # 0.708
    "SLIK4",                    # 0.602
    "Amyloid-like protein 1",   # 0.561
    # Bottom 5 by D1-D3 ipTM
    "NGL1",                     # 0.173
    "BASI",                     # 0.177
    "SLIT2",                    # 0.202
    "SEMA3A",                   # 0.205
    "Contactin-5",              # 0.206
    # Known false negative (VEGF-bridged ternary)
    "NRP1",                     # 0.222
]

# D1-D7: skip largest partners to avoid A100 OOM
# PCDH9 (~1300aa + 721 = 2021), SLIT2 (~1500aa + 721 = 2221) risky
D1D7_SKIP = {"PCDH9", "SLIT2"}


def _walltime_for_residues(total_residues: int) -> str:
    """Assign walltime based on total complex residues."""
    if total_residues < 800:
        return "48:00"
    if total_residues < 1200:
        return "72:00"
    if total_residues < 1800:
        return "96:00"
    return "144:00"


def _make_candidates_csv(targets: list[str], source_csv: Path, output_path: Path) -> Path:
    """Filter source candidates CSV to a target subset."""
    import pandas as pd

    all_candidates = pd.read_csv(source_csv)
    filtered = all_candidates[all_candidates["target"].isin(targets)].copy()

    missing = set(targets) - set(filtered["target"])
    if missing:
        logger.warning("Targets not found in source CSV: %s", missing)

    filtered.to_csv(output_path, index=False)
    logger.info("Wrote %d candidates to %s", len(filtered), output_path)
    return output_path


def generate_batch(
    batch_name: str,
    construct,
    targets: list[str],
    source_candidates_csv: Path,
    base_dir: Path,
    hpc_root: str = "",
    project_account: str = "acc_vascbrain",
    gpu_type: str = '"a100"',
    dry_run: bool = False,
) -> dict:
    """Generate FASTAs and LSF jobs for one batch.

    Returns summary dict with counts and paths.
    """
    batch_dir = base_dir / batch_name
    fasta_dir = batch_dir / "fasta"
    jobs_dir = batch_dir / "jobs"
    candidates_csv = batch_dir / "candidates.csv"

    batch_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Batch: %s (%s, %d targets)", batch_name, construct.name, len(targets))
    logger.info("  Construct: %s (%d aa, residues %d-%d)",
                construct.name, construct.length, construct.start, construct.end)
    logger.info("=" * 60)

    # 1. Filter candidates
    _make_candidates_csv(targets, source_candidates_csv, candidates_csv)

    # 2. Fetch sequences + write FASTAs
    results = fetch_all_sequences(candidates_csv, fasta_dir, construct=construct)

    ok_targets = {k: v for k, v in results.items() if v["status"] == "ok"}
    failed = {k: v for k, v in results.items() if v["status"] != "ok"}
    if failed:
        logger.warning("Failed fetches: %s", list(failed.keys()))

    # 3. Compute per-target walltimes
    walltime_map = {}
    for target, info in ok_targets.items():
        total = info["total_residues"]
        wt = _walltime_for_residues(total)
        walltime_map[target] = {"total_residues": total, "walltime": wt}
        logger.info("  %s: %d total residues -> %s walltime", target, total, wt)

    if dry_run:
        logger.info("[DRY RUN] Would generate %d LSF jobs for batch %s",
                     len(ok_targets), batch_name)
        return {"batch": batch_name, "construct": construct.name,
                "n_targets": len(ok_targets), "walltime_map": walltime_map}

    # 4. Generate LSF jobs with per-target walltimes
    #    We need to patch the walltime per job. The existing generator uses a
    #    single walltime. We'll generate with max walltime, then patch each
    #    .lsf file individually.
    hpc_batch_root = ""
    if hpc_root:
        hpc_batch_root = f"{hpc_root}/results/03_structural_prediction/{batch_name}"

    # Set module-level walltime to max for this batch
    max_wt = max(
        (v["walltime"] for v in walltime_map.values()),
        default="72:00",
    )
    original_walltime = _mod03.AF_WALLTIME
    _mod03.AF_WALLTIME = max_wt

    scripts = _mod03.generate_lsf_scripts(
        candidates_csv, fasta_dir, jobs_dir,
        project_account=project_account,
        gpu_type=gpu_type,
        hpc_root=hpc_batch_root if hpc_batch_root else "",
    )

    _mod03.AF_WALLTIME = original_walltime  # restore

    # 5. Patch per-job walltimes in .lsf files
    for script_path in scripts:
        target_name = script_path.stem.replace("af2_", "")
        # Match against walltime_map keys (need safe_name -> original name mapping)
        matched_wt = None
        for orig_target, wt_info in walltime_map.items():
            safe = orig_target.replace("/", "-").replace(" ", "_").replace(":", "_")
            if safe == target_name:
                matched_wt = wt_info["walltime"]
                break

        if matched_wt and matched_wt != max_wt:
            content = script_path.read_text()
            content = content.replace(f"-W {max_wt}", f"-W {matched_wt}")
            script_path.write_text(content)
            logger.info("  Patched %s walltime: %s -> %s",
                        script_path.name, max_wt, matched_wt)

    # Also patch walltime in manifest.json
    manifest_path = jobs_dir / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for job_key, job_entry in manifest.items():
            target_name = job_key.replace("af2_", "")
            for orig_target, wt_info in walltime_map.items():
                safe = orig_target.replace("/", "-").replace(" ", "_").replace(":", "_")
                if safe == target_name:
                    job_entry["walltime"] = wt_info["walltime"]
                    job_entry["total_residues"] = wt_info["total_residues"]
                    break
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))

    # 6. Write batch summary
    summary = {
        "batch": batch_name,
        "construct": construct.name,
        "construct_residues": f"{construct.start}-{construct.end} ({construct.length} aa)",
        "n_targets": len(ok_targets),
        "n_scripts": len(scripts),
        "walltime_map": walltime_map,
        "hpc_root": hpc_batch_root,
    }
    summary_path = batch_dir / "batch_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Batch summary: %s", summary_path)

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Generate corrected D1-D3 + D1-D6 + D1-D7 AF2 job batches"
    )
    parser.add_argument(
        "--hpc-root",
        type=str,
        default="/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold",
        help="Absolute project root on HPC",
    )
    parser.add_argument(
        "--project-account",
        type=str,
        default="acc_vascbrain",
    )
    parser.add_argument(
        "--gpu-type",
        type=str,
        default='"a100"',
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    project_root = SCRIPT_DIR.parent.parent  # vasc-sflt1-alphafold/
    results_dir = project_root / "results" / "03_structural_prediction"
    source_csv = results_dir / "step03_candidates.csv"

    if not source_csv.exists():
        logger.error("Source candidates not found: %s", source_csv)
        logger.info("Run step 3a first (--step select)")
        sys.exit(1)

    summaries = []

    # Batch A: D1-D3 corrected -- all 24 targets
    summaries.append(generate_batch(
        batch_name="d1d3_corrected",
        construct=CONSTRUCT_D1D3,
        targets=ALL_24_TARGETS,
        source_candidates_csv=source_csv,
        base_dir=results_dir,
        hpc_root=args.hpc_root,
        project_account=args.project_account,
        gpu_type=args.gpu_type,
        dry_run=args.dry_run,
    ))

    # Batch B: D1-D6 -- 11 EC targets
    summaries.append(generate_batch(
        batch_name="d1d6",
        construct=CONSTRUCT_D1D6,
        targets=EC_TARGETS,
        source_candidates_csv=source_csv,
        base_dir=results_dir,
        hpc_root=args.hpc_root,
        project_account=args.project_account,
        gpu_type=args.gpu_type,
        dry_run=args.dry_run,
    ))

    # Batch C: D1-D7 -- EC targets minus VRAM-risky large partners
    d1d7_targets = [t for t in EC_TARGETS if t not in D1D7_SKIP]
    summaries.append(generate_batch(
        batch_name="d1d7",
        construct=CONSTRUCT_D1D7,
        targets=d1d7_targets,
        source_candidates_csv=source_csv,
        base_dir=results_dir,
        hpc_root=args.hpc_root,
        project_account=args.project_account,
        gpu_type=args.gpu_type,
        dry_run=args.dry_run,
    ))

    # Print summary
    logger.info("")
    logger.info("=" * 60)
    logger.info("BATCH GENERATION COMPLETE")
    logger.info("=" * 60)
    total_jobs = 0
    for s in summaries:
        n = s["n_targets"]
        total_jobs += n
        logger.info("  %s: %d jobs (%s, %d aa construct)",
                     s["batch"], n, s["construct"],
                     CONSTRUCTS[s["construct"].lower().replace("-", "")].length
                     if s["construct"].lower().replace("-", "") in CONSTRUCTS else 0)
    logger.info("  TOTAL: %d jobs", total_jobs)
    logger.info("")
    logger.info("Next steps:")
    logger.info("  1. git add + commit + push")
    logger.info("  2. ssh minerva 'cd <project> && git pull'")
    logger.info("  3. Submit batch A first (priority): bash results/03_structural_prediction/d1d3_corrected/jobs/submit_all.sh")
    logger.info("  4. Submit batch B: bash results/03_structural_prediction/d1d6/jobs/submit_all.sh")
    logger.info("  5. Submit batch C (deprioritized): bash results/03_structural_prediction/d1d7/jobs/submit_all.sh")


if __name__ == "__main__":
    main()
