#!/usr/bin/env python3
"""10_merge_scores.py -- Merge per-target scores.json into interaction_scores.csv.

Walks the AF2 output directory, collects all scores.json sidecars written by
09_extract_and_cleanup.py, and produces the final step03_interaction_scores.csv.

Designed for HPC: stdlib only (no pandas/numpy dependency).

Usage:
    python3 10_merge_scores.py <af_output_dir> <candidates_csv> <output_csv>

    af_output_dir:   Root AF2 results directory (contains per-target subdirs)
    candidates_csv:  step03_candidates.csv (for target -> uniprot mapping)
    output_csv:      Output path for merged interaction scores
"""

import csv
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


def safe_name(target: str) -> str:
    return target.replace("/", "-").replace(" ", "_").replace(":", "_")


def merge_scores(
    af_output_dir: Path,
    candidates_csv: Path,
    output_csv: Path,
) -> None:
    """Merge all per-target scores.json into a single CSV.

    Parameters
    ----------
    af_output_dir : Path
        Root directory containing per-target AF2 output subdirs.
    candidates_csv : Path
        CSV with target, uniprot columns.
    output_csv : Path
        Output path for merged scores.
    """
    # Read candidates for target -> uniprot mapping
    with open(candidates_csv, newline="") as f:
        candidates = list(csv.DictReader(f))

    uniprot_map = {c["target"]: c.get("uniprot", "") for c in candidates}

    fieldnames = [
        "target", "uniprot", "iptm_best", "iptm_mean", "ptm_best",
        "mean_interchain_pae", "mean_interface_plddt", "interaction_call",
        "n_interface_residues", "n_models_parsed",
    ]

    rows = []
    n_found = 0
    n_missing = 0

    for cand in candidates:
        target = cand["target"]
        sname = safe_name(target)
        uniprot = uniprot_map.get(target, "")

        # Find scores.json in the AF2 output subdir
        target_dir = af_output_dir / sname
        nested = target_dir / f"sflt1_vs_{sname}"
        actual = nested if nested.exists() else target_dir

        scores_path = actual / "scores.json" if actual.exists() else None

        if scores_path and scores_path.exists():
            with open(scores_path) as f:
                scores = json.load(f)
            n_found += 1

            def fmt(v):
                return v if v is not None else ""

            rows.append({
                "target": target,
                "uniprot": uniprot,
                "iptm_best": fmt(scores.get("iptm_best")),
                "iptm_mean": fmt(scores.get("iptm_mean")),
                "ptm_best": fmt(scores.get("ptm_best")),
                "mean_interchain_pae": fmt(scores.get("mean_interchain_pae")),
                "mean_interface_plddt": fmt(scores.get("mean_interface_plddt")),
                "interaction_call": scores.get("interaction_call", "no_data"),
                "n_interface_residues": scores.get("n_interface_residues", 0),
                "n_models_parsed": scores.get("n_models_parsed", 0),
            })
        else:
            n_missing += 1
            rows.append({
                "target": target,
                "uniprot": uniprot,
                "iptm_best": "",
                "iptm_mean": "",
                "ptm_best": "",
                "mean_interchain_pae": "",
                "mean_interface_plddt": "",
                "interaction_call": "not_run",
                "n_interface_residues": 0,
                "n_models_parsed": 0,
            })
            logger.warning("No scores.json for %s", target)

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Merged %d targets (%d with scores, %d missing) -> %s",
                len(rows), n_found, n_missing, output_csv)


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <af_output_dir> <candidates_csv> <output_csv>",
              file=sys.stderr)
        sys.exit(1)

    merge_scores(
        af_output_dir=Path(sys.argv[1]),
        candidates_csv=Path(sys.argv[2]),
        output_csv=Path(sys.argv[3]),
    )


if __name__ == "__main__":
    main()
