#!/usr/bin/env python3
"""check_af2_inputs.py -- Verify what AlphaFold will consume at runtime.

Decodes the base64 command from each LSF script and checks:
  1. AF2 flags: model_preset, db_preset, max_template_date, num_predictions
  2. Database paths referenced in the command
  3. FASTA path referenced in the command exists locally (and on HPC if accessible)
  4. Output directory path is inside the correct batch
  5. Consistent settings across all jobs in a batch

Run locally:  python analysis/checks/check_af2_inputs.py
Run on HPC:   python analysis/checks/check_af2_inputs.py --check-hpc-paths

The --check-hpc-paths flag verifies that FASTA files and database dirs exist
on the filesystem (only works when run on Minerva).
"""

import argparse
import base64
import json
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "03_structural_prediction"

BATCHES = ["d1d3", "d1d6", "d1d7"]

# Expected AF2 settings
EXPECTED_FLAGS = {
    "model_preset": "multimer",
    "db_preset": "reduced_dbs",
    "max_template_date": "2024-01-01",
    "num_multimer_predictions_per_model": "5",
}

# Expected databases (container paths)
EXPECTED_DBS = [
    "/data/uniref90/uniref90.fasta",
    "/data/mgnify/mgy_clusters_2022_05.fa",
    "/data/pdb_mmcif/mmcif_files",
    "/data/pdb_mmcif/obsolete.dat",
    "/data/pdb_seqres/pdb_seqres.txt",
    "/data/uniprot/uniprot.fasta",
    "/data/small_bfd/bfd-first_non_consensus_sequences.fasta",
]


def decode_lsf_command(lsf_path: Path) -> str | None:
    """Extract and decode the base64 command from an LSF script."""
    text = lsf_path.read_text()
    match = re.search(r'AF_JOB_COMMAND_B64="([A-Za-z0-9+/=\n]+)"', text)
    if not match:
        return None
    b64 = match.group(1).replace("\n", "")
    return base64.b64decode(b64).decode("utf-8")


def extract_flag(command: str, flag_name: str) -> str | None:
    """Extract --flag=value from a command string."""
    pattern = rf"--{flag_name}=(\S+)"
    match = re.search(pattern, command)
    return match.group(1).strip('"').strip("'") if match else None


def extract_fasta_path(command: str) -> str | None:
    """Extract the FASTA_INPUT assignment from the command."""
    match = re.search(r'FASTA_INPUT="([^"]+)"', command)
    return match.group(1) if match else None


def extract_output_dir(command: str) -> str | None:
    """Extract the OUTPUT_DIR assignment from the command."""
    match = re.search(r'OUTPUT_DIR="([^"]+)"', command)
    return match.group(1) if match else None


def extract_hpc_root(command: str) -> str | None:
    """Extract the cd target directory from the command."""
    match = re.search(r'^cd "([^"]+)"', command, re.MULTILINE)
    return match.group(1) if match else None


def check_batch(batch_name: str, check_hpc: bool) -> list[str]:
    """Check all LSF scripts in one batch. Returns failure messages."""
    failures = []
    jobs_dir = RESULTS_DIR / batch_name / "jobs"

    if not jobs_dir.exists():
        failures.append(f"{batch_name}: jobs/ directory missing")
        return failures

    lsf_files = sorted(jobs_dir.glob("af2_*.lsf"))
    if not lsf_files:
        failures.append(f"{batch_name}: no af2_*.lsf files found")
        return failures

    print(f"\n  {batch_name}: {len(lsf_files)} LSF scripts")

    # Collect settings across all jobs to check consistency
    all_settings = {}  # flag -> set of values seen
    all_dbs_seen = {}  # db_path -> count of jobs referencing it
    fasta_paths = {}   # target -> fasta relative path
    output_dirs = {}   # target -> output relative path
    hpc_roots = set()

    for lsf_path in lsf_files:
        target = lsf_path.stem.replace("af2_", "")
        command = decode_lsf_command(lsf_path)

        if command is None:
            failures.append(f"{batch_name}/{target}: cannot decode base64 command")
            continue

        # Extract and collect flags
        for flag, expected in EXPECTED_FLAGS.items():
            actual = extract_flag(command, flag)
            if actual is None:
                failures.append(f"{batch_name}/{target}: missing --{flag}")
            elif actual != expected:
                failures.append(
                    f"{batch_name}/{target}: --{flag}={actual}, expected {expected}"
                )
            all_settings.setdefault(flag, set()).add(actual)

        # Check database paths are referenced
        for db_path in EXPECTED_DBS:
            if db_path in command:
                all_dbs_seen[db_path] = all_dbs_seen.get(db_path, 0) + 1
            else:
                failures.append(f"{batch_name}/{target}: missing database ref {db_path}")

        # Extract paths
        fasta_rel = extract_fasta_path(command)
        output_rel = extract_output_dir(command)
        hpc_root = extract_hpc_root(command)

        if fasta_rel:
            fasta_paths[target] = fasta_rel
        if output_rel:
            output_dirs[target] = output_rel
        if hpc_root:
            hpc_roots.add(hpc_root)

        # Check FASTA path points into fasta/ subdir (not root fasta/)
        if fasta_rel and "fasta/sflt1_vs_" not in fasta_rel:
            failures.append(
                f"{batch_name}/{target}: FASTA path doesn't match expected pattern: {fasta_rel}"
            )

        # Check output dir points into results/ subdir
        if output_rel and not output_rel.startswith("results/"):
            failures.append(
                f"{batch_name}/{target}: output dir doesn't start with results/: {output_rel}"
            )

        # If checking HPC paths, verify files exist
        if check_hpc and hpc_root and fasta_rel:
            full_fasta = Path(hpc_root) / fasta_rel
            if not full_fasta.exists():
                failures.append(f"{batch_name}/{target}: FASTA not found on HPC: {full_fasta}")

    # Consistency checks
    for flag, values in all_settings.items():
        if len(values) > 1:
            failures.append(
                f"{batch_name}: inconsistent --{flag} across jobs: {values}"
            )

    if len(hpc_roots) > 1:
        failures.append(f"{batch_name}: inconsistent HPC roots: {hpc_roots}")

    # Print summary
    print(f"    HPC root: {next(iter(hpc_roots), 'not set')}")
    print(f"    AF2 flags: {', '.join(f'{k}={next(iter(v))}' for k, v in all_settings.items() if v)}")
    print(f"    Databases: {len(all_dbs_seen)}/{len(EXPECTED_DBS)} referenced")
    print(f"    FASTA paths: {len(fasta_paths)}")
    print(f"    Output dirs: {len(output_dirs)}")

    # Check database path consistency with db_preset
    db_preset = next(iter(all_settings.get("db_preset", set())), None)
    if db_preset == "reduced_dbs":
        # reduced_dbs should use small_bfd, NOT full bfd
        for lsf_path in lsf_files:
            command = decode_lsf_command(lsf_path)
            if command and "bfd_uniclust" in command:
                target = lsf_path.stem.replace("af2_", "")
                failures.append(
                    f"{batch_name}/{target}: reduced_dbs but references full BFD"
                )

    return failures


def main():
    parser = argparse.ArgumentParser(description="Verify AF2 runtime inputs")
    parser.add_argument(
        "--check-hpc-paths",
        action="store_true",
        help="Verify FASTA files exist on HPC filesystem (run on Minerva only)",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("AlphaFold Input Preflight Check")
    print("=" * 60)
    print(f"\n  Expected AF2 settings:")
    for k, v in EXPECTED_FLAGS.items():
        print(f"    --{k}={v}")
    print(f"\n  Expected databases ({len(EXPECTED_DBS)}):")
    for db in EXPECTED_DBS:
        print(f"    {db}")

    all_failures = []

    for batch_name in BATCHES:
        batch_dir = RESULTS_DIR / batch_name
        if not batch_dir.exists():
            print(f"\n  {batch_name}: SKIP (directory not found)")
            continue
        failures = check_batch(batch_name, args.check_hpc_paths)
        all_failures.extend(failures)

    # Summary
    print("\n" + "=" * 60)
    if all_failures:
        print(f"FAIL: {len(all_failures)} issues found\n")
        for f in all_failures:
            print(f"  FAIL: {f}")
        sys.exit(1)
    else:
        print("PASS: all AF2 inputs verified")
        if not args.check_hpc_paths:
            print("  (run with --check-hpc-paths on Minerva to verify filesystem paths)")
        sys.exit(0)


if __name__ == "__main__":
    main()
