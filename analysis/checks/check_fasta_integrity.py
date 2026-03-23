#!/usr/bin/env python3
"""check_fasta_integrity.py -- Verify FASTA files match expected sequences and lengths.

For each batch (d1d3_corrected, d1d6, d1d7):
  1. sFLT1 chain has correct length (304/631/721 aa) and no signal peptide
  2. Target chain lengths match batch_summary.json total_residues - construct_length
  3. No signal peptide residues in target chains (cross-ref UniProt Chain boundaries)
  4. Headers have correct residue range annotations

Run: python analysis/checks/check_fasta_integrity.py
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "03_structural_prediction"

CONSTRUCTS = {
    "d1d3_corrected": {"start": 27, "end": 330, "length": 304},
    "d1d6": {"start": 27, "end": 657, "length": 631},
    "d1d7": {"start": 27, "end": 747, "length": 721},
}

# Known signal peptide boundaries (from UniProt, verified 2026-03-22)
# Only for targets that are secreted/soluble (not TM-extracted)
SIGNAL_PEPTIDES = {
    "Q14563": {"name": "SEMA3A", "signal_end": 20, "mature_start": 21},
    "O94779": {"name": "Contactin-5", "signal_end": 18, "mature_start": 19, "chain_end": 1072},
    "Q6MZW2": {"name": "FSTL4", "signal_end": 22, "mature_start": 23},
    "Q99784": {"name": "NOE1", "signal_end": 16, "mature_start": 17},
    "O94813": {"name": "SLIT2", "signal_end": 30, "mature_start": 31},
}


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse a FASTA file into list of (header, sequence) tuples."""
    entries = []
    current_header = None
    current_seq = []

    for line in path.read_text().strip().split("\n"):
        if line.startswith(">"):
            if current_header is not None:
                entries.append((current_header, "".join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())

    if current_header is not None:
        entries.append((current_header, "".join(current_seq)))

    return entries


def check_batch(batch_name: str, construct: dict) -> list[str]:
    """Check all FASTA files in one batch. Returns list of failure messages."""
    failures = []
    batch_dir = RESULTS_DIR / batch_name
    fasta_dir = batch_dir / "fasta"
    summary_path = batch_dir / "batch_summary.json"

    if not fasta_dir.exists():
        failures.append(f"{batch_name}: fasta/ directory missing")
        return failures

    # Load batch summary for expected residue counts
    walltime_map = {}
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        walltime_map = summary.get("walltime_map", {})

    fasta_files = sorted(fasta_dir.glob("sflt1_vs_*.fasta"))
    if not fasta_files:
        failures.append(f"{batch_name}: no sflt1_vs_*.fasta files found")
        return failures

    print(f"\n  {batch_name}: {len(fasta_files)} FASTA files, construct {construct['length']} aa")

    for fasta_path in fasta_files:
        target_name = fasta_path.stem.replace("sflt1_vs_", "")
        entries = parse_fasta(fasta_path)

        if len(entries) != 2:
            failures.append(f"{batch_name}/{fasta_path.name}: expected 2 chains, got {len(entries)}")
            continue

        sflt1_header, sflt1_seq = entries[0]
        target_header, target_seq = entries[1]

        # Check 1: sFLT1 chain length
        expected_len = construct["length"]
        if len(sflt1_seq) != expected_len:
            failures.append(
                f"{batch_name}/{target_name}: sFLT1 length {len(sflt1_seq)} != expected {expected_len}"
            )

        # Check 2: sFLT1 header has correct residue range
        expected_range = f"residues_{construct['start']}-{construct['end']}"
        if expected_range not in sflt1_header:
            failures.append(
                f"{batch_name}/{target_name}: sFLT1 header missing '{expected_range}', "
                f"got '{sflt1_header}'"
            )

        # Check 3: sFLT1 does NOT start with signal peptide (Met-Val-Ser pattern)
        # Mature sFLT1 starts with Ser27 (S), not Met1 (M)
        if sflt1_seq and sflt1_seq[0] == "M" and len(sflt1_seq) > 300:
            failures.append(
                f"{batch_name}/{target_name}: sFLT1 starts with M -- possible signal peptide included"
            )

        # Check 4: target chain length matches batch_summary total_residues
        # Restore original target name (undo safe_name for lookup)
        original_name = target_name.replace("_", " ")
        if original_name not in walltime_map:
            # Try exact match
            original_name = target_name
        if original_name in walltime_map:
            expected_total = walltime_map[original_name]["total_residues"]
            actual_total = len(sflt1_seq) + len(target_seq)
            if actual_total != expected_total:
                failures.append(
                    f"{batch_name}/{target_name}: total residues {actual_total} != "
                    f"batch_summary {expected_total}"
                )

        # Check 5: target UniProt in header
        parts = target_header.split("|")
        if len(parts) < 2:
            failures.append(
                f"{batch_name}/{target_name}: target header missing UniProt accession"
            )
            continue

        uniprot = parts[1].strip()

        # Check 6: signal peptide not included for known secreted targets
        if uniprot in SIGNAL_PEPTIDES:
            sp = SIGNAL_PEPTIDES[uniprot]
            # If the sequence is full-length (no ectodomain/mature tag in header),
            # it might still have the signal peptide
            region_tag = parts[2].strip() if len(parts) > 2 else ""
            if not region_tag:
                # No region annotation -- suspicious for a secreted protein
                failures.append(
                    f"{batch_name}/{target_name}: secreted protein {sp['name']} ({uniprot}) "
                    f"has no region annotation in header -- may include signal peptide"
                )
            elif "mature" in region_tag or "ectodomain" in region_tag:
                # Has region tag -- verify start residue is past signal peptide
                # Parse "mature_21-771" or "ectodomain_22-856"
                try:
                    range_part = region_tag.split("_", 1)[1]
                    start_res = int(range_part.split("-")[0])
                    if start_res <= sp["signal_end"]:
                        failures.append(
                            f"{batch_name}/{target_name}: region starts at {start_res} "
                            f"but signal peptide ends at {sp['signal_end']}"
                        )
                except (IndexError, ValueError):
                    pass

    return failures


def main():
    print("=" * 60)
    print("FASTA Integrity Preflight Check")
    print("=" * 60)

    all_failures = []

    for batch_name, construct in CONSTRUCTS.items():
        batch_dir = RESULTS_DIR / batch_name
        if not batch_dir.exists():
            print(f"\n  {batch_name}: SKIP (directory not found)")
            continue
        failures = check_batch(batch_name, construct)
        all_failures.extend(failures)

    # Also check root fasta/ for stale artifacts
    print("\n  Stale artifact check (root fasta/):")
    root_fasta = RESULTS_DIR / "fasta"
    if root_fasta.exists():
        for f in sorted(root_fasta.glob("sflt1_vs_*.fasta")):
            entries = parse_fasta(f)
            if entries:
                sflt1_hdr = entries[0][0]
                if "residues_1-338" in sflt1_hdr:
                    print(f"    WARN: {f.name} has stale construct (residues_1-338)")

    # Summary
    print("\n" + "=" * 60)
    if all_failures:
        print(f"FAIL: {len(all_failures)} issues found\n")
        for f in all_failures:
            print(f"  FAIL: {f}")
        sys.exit(1)
    else:
        print("PASS: all FASTA files verified")
        sys.exit(0)


if __name__ == "__main__":
    main()
