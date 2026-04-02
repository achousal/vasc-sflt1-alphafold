#!/usr/bin/env python3
"""patch_missing_uniprot.py -- Fill missing UniProt accessions in consensus CSV.

Resolves 9 SomaScan target names to their canonical UniProt accessions.
Writes patched CSV alongside the original (does not overwrite).

Run: python analysis/checks/patch_missing_uniprot.py
     python analysis/checks/patch_missing_uniprot.py --dry-run
"""

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CONSENSUS_PATH = (
    PROJECT_ROOT / "results" / "01_cross_cohort_overlap" / "step01_consensus_proteins_pos.csv"
)

# Manually resolved UniProt accessions for 9 targets missing in the consensus CSV.
# Each mapping verified against UniProt Swiss-Prot (reviewed) entries for Homo sapiens.
#
# Disambiguation notes:
#   GFRP  -> GCHFR (P30047): GFRP is a protein alias for GTP cyclohydrolase 1
#            feedback regulatory protein; gene symbol is GCHFR.
#   GNA1  -> GNAI1 (P63096): SomaScan abbreviation for G(i) alpha-1 subunit.
#            Not to be confused with P63098 (PPP3R1, Calcineurin B a).
#   NRX3A -> NRXN3 (Q9Y4C0): SomaScan shorthand for neurexin-3-alpha isoform.
#            Q9HDB5 is the beta isoform; "A" suffix indicates alpha.
#   Lysosomal acid phosphatase -> ACP2 (P11117): descriptive SomaScan name.
UNIPROT_PATCHES = {
    "BID":                       {"uniprot": "P55957", "entrez_gene_symbol": "BID",    "entrez_gene_id": 637},
    "GFRP":                      {"uniprot": "P30047", "entrez_gene_symbol": "GCHFR",  "entrez_gene_id": 2644},
    "WBP2":                      {"uniprot": "Q969T9", "entrez_gene_symbol": "WBP2",   "entrez_gene_id": 23558},
    "GNA1":                      {"uniprot": "P63096", "entrez_gene_symbol": "GNAI1",  "entrez_gene_id": 2770},
    "Lysosomal acid phosphatase": {"uniprot": "P11117", "entrez_gene_symbol": "ACP2",  "entrez_gene_id": 53},
    "LSAMP":                     {"uniprot": "Q13449", "entrez_gene_symbol": "LSAMP",  "entrez_gene_id": 4045},
    "CBLN2":                     {"uniprot": "Q8IUK8", "entrez_gene_symbol": "CBLN2",  "entrez_gene_id": 147381},
    "NRX3A":                     {"uniprot": "Q9Y4C0", "entrez_gene_symbol": "NRXN3",  "entrez_gene_id": 9369},
    "ENPP6":                     {"uniprot": "Q6UWR7", "entrez_gene_symbol": "ENPP6",  "entrez_gene_id": 28978},
}


def main():
    parser = argparse.ArgumentParser(description="Patch missing UniProt accessions")
    parser.add_argument("--dry-run", action="store_true", help="Show changes without writing")
    args = parser.parse_args()

    if not CONSENSUS_PATH.exists():
        print(f"FAIL: consensus file not found: {CONSENSUS_PATH}")
        sys.exit(1)

    df = pd.read_csv(CONSENSUS_PATH)
    missing_before = df["uniprot"].isna().sum()
    print(f"Loaded {len(df)} rows, {missing_before} with missing UniProt")

    patched = 0
    for target, patch in UNIPROT_PATCHES.items():
        mask = df["target"] == target
        n_rows = mask.sum()
        if n_rows == 0:
            print(f"  WARNING: target '{target}' not found in CSV")
            continue

        for col, val in patch.items():
            df.loc[mask, col] = val
        patched += n_rows
        print(f"  PATCH: {target} ({n_rows} rows) -> {patch['uniprot']} ({patch['entrez_gene_symbol']})")

    missing_after = df["uniprot"].isna().sum()
    print(f"\nPatched {patched} rows. Missing UniProt: {missing_before} -> {missing_after}")

    if args.dry_run:
        print("\n(dry run -- no files written)")
        return

    # Backup original
    backup_path = CONSENSUS_PATH.with_suffix(".csv.bak")
    if not backup_path.exists():
        shutil.copy2(CONSENSUS_PATH, backup_path)
        print(f"Backup: {backup_path.name}")

    # Write patched CSV
    df.to_csv(CONSENSUS_PATH, index=False)
    print(f"Written: {CONSENSUS_PATH.name}")


if __name__ == "__main__":
    main()
