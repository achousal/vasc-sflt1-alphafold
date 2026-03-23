#!/usr/bin/env python3
"""check_candidates_complete.py -- Verify all candidates have FASTA + LSF in every batch.

Checks:
  1. Every target in step03_candidates.csv has a FASTA file in each batch it belongs to
  2. Every FASTA file has a matching .lsf script
  3. Every .lsf script has a matching manifest.json entry
  4. No orphan files (FASTA without candidate, LSF without FASTA)
  5. Batch target lists match expected subsets

Run: python analysis/checks/check_candidates_complete.py
"""

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "03_structural_prediction"
CANDIDATES_PATH = RESULTS_DIR / "step03_candidates.csv"

# Expected batch membership (from 06_generate_ec_batches.py)
EC_TARGETS = {
    "VEGFA", "PCDH9", "NOE1", "SLIK4", "Amyloid-like protein 1",
    "NGL1", "BASI", "SLIT2", "SEMA3A", "Contactin-5", "NRP1",
}
D1D7_SKIP = {"PCDH9", "SLIT2"}

BATCH_EXPECTED_TARGETS = {
    "d1d3": None,  # all 24
    "d1d6": EC_TARGETS,
    "d1d7": EC_TARGETS - D1D7_SKIP,
}


def safe_name(target: str) -> str:
    return target.replace("/", "-").replace(" ", "_").replace(":", "_")


def check_batch(batch_name: str, candidates: pd.DataFrame, expected_targets: set | None) -> list[str]:
    """Check one batch for completeness. Returns failure messages."""
    failures = []
    batch_dir = RESULTS_DIR / batch_name
    fasta_dir = batch_dir / "fasta"
    jobs_dir = batch_dir / "jobs"
    manifest_path = jobs_dir / "manifest.json"

    if not batch_dir.exists():
        failures.append(f"{batch_name}: batch directory missing")
        return failures

    # Determine expected targets
    if expected_targets is None:
        target_set = set(candidates["target"])
    else:
        target_set = expected_targets

    print(f"\n  {batch_name}: expecting {len(target_set)} targets")

    # Check FASTA files
    fasta_targets = set()
    if fasta_dir.exists():
        for f in fasta_dir.glob("sflt1_vs_*.fasta"):
            name = f.stem.replace("sflt1_vs_", "")
            fasta_targets.add(name)

    # Check LSF files
    lsf_targets = set()
    if jobs_dir.exists():
        for f in jobs_dir.glob("af2_*.lsf"):
            name = f.stem.replace("af2_", "")
            lsf_targets.add(name)

    # Check manifest
    manifest_targets = set()
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for key in manifest:
            name = key.replace("af2_", "")
            manifest_targets.add(name)

    # Build safe_name -> original name mapping
    safe_to_original = {safe_name(t): t for t in target_set}
    expected_safe = set(safe_to_original.keys())

    # 1. Missing FASTA files
    missing_fasta = expected_safe - fasta_targets
    for m in sorted(missing_fasta):
        failures.append(f"{batch_name}: missing FASTA for {safe_to_original.get(m, m)}")

    # 2. Missing LSF scripts
    missing_lsf = expected_safe - lsf_targets
    for m in sorted(missing_lsf):
        failures.append(f"{batch_name}: missing LSF for {safe_to_original.get(m, m)}")

    # 3. Missing manifest entries
    missing_manifest = expected_safe - manifest_targets
    for m in sorted(missing_manifest):
        failures.append(f"{batch_name}: missing manifest entry for {safe_to_original.get(m, m)}")

    # 4. Orphan FASTA files (not in expected set)
    orphan_fasta = fasta_targets - expected_safe - {"sflt1_d1d3_reference", "sflt1_d1d6_reference", "sflt1_d1d7_reference"}
    for o in sorted(orphan_fasta):
        failures.append(f"{batch_name}: orphan FASTA file sflt1_vs_{o}.fasta (not in candidate list)")

    # 5. Orphan LSF files
    orphan_lsf = lsf_targets - expected_safe
    for o in sorted(orphan_lsf):
        failures.append(f"{batch_name}: orphan LSF file af2_{o}.lsf (not in candidate list)")

    # 6. FASTA without matching LSF
    fasta_no_lsf = fasta_targets - lsf_targets - {"sflt1_d1d3_reference", "sflt1_d1d6_reference", "sflt1_d1d7_reference"}
    for f in sorted(fasta_no_lsf):
        failures.append(f"{batch_name}: FASTA exists but no LSF for {f}")

    # 7. LSF without matching FASTA
    lsf_no_fasta = lsf_targets - fasta_targets
    for f in sorted(lsf_no_fasta):
        failures.append(f"{batch_name}: LSF exists but no FASTA for {f}")

    # Print counts
    print(f"    FASTA files: {len(fasta_targets)}")
    print(f"    LSF scripts: {len(lsf_targets)}")
    print(f"    Manifest entries: {len(manifest_targets)}")

    return failures


def main():
    print("=" * 60)
    print("Candidate Completeness Preflight Check")
    print("=" * 60)

    if not CANDIDATES_PATH.exists():
        print(f"\nFAIL: candidates file not found: {CANDIDATES_PATH}")
        sys.exit(1)

    candidates = pd.read_csv(CANDIDATES_PATH)
    print(f"\n  Candidates CSV: {len(candidates)} targets")

    # Verify candidates CSV has required columns
    required_cols = {"target", "uniprot", "rank"}
    missing_cols = required_cols - set(candidates.columns)
    if missing_cols:
        print(f"\nFAIL: candidates CSV missing columns: {missing_cols}")
        sys.exit(1)

    # Verify no duplicate targets
    dupes = candidates[candidates["target"].duplicated()]["target"].tolist()
    if dupes:
        print(f"\nFAIL: duplicate targets in candidates CSV: {dupes}")
        sys.exit(1)

    # Verify UniProt accessions are non-empty
    empty_uniprot = candidates[candidates["uniprot"].isna() | (candidates["uniprot"] == "")]["target"].tolist()
    if empty_uniprot:
        print(f"\nFAIL: targets with missing UniProt: {empty_uniprot}")
        sys.exit(1)

    all_failures = []

    for batch_name, expected in BATCH_EXPECTED_TARGETS.items():
        failures = check_batch(batch_name, candidates, expected)
        all_failures.extend(failures)

    # Summary
    print("\n" + "=" * 60)
    if all_failures:
        print(f"FAIL: {len(all_failures)} issues found\n")
        for f in all_failures:
            print(f"  FAIL: {f}")
        sys.exit(1)
    else:
        print("PASS: all candidates complete across all batches")
        sys.exit(0)


if __name__ == "__main__":
    main()
