#!/usr/bin/env python3
"""10_merge_scores.py -- Merge per-target scores.json into interaction_scores.csv.

Walks the AF2 output directory, collects all scores.json sidecars written by
09_extract_and_cleanup.py, and produces the final step03_interaction_scores.csv.

scores.json is the consolidated source depot for all AF2 metrics. This script
reads the "summary" block and flattens it into one CSV row per target.

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

# All summary fields to export, in column order
SUMMARY_FIELDS = [
    "iptm_best", "iptm_mean", "iptm_std", "ptm_best",
    "interchain_pae_best", "interchain_pae_mean", "interchain_pae_std",
    "interface_plddt_best", "interface_plddt_mean",
    "n_interface_residues_best",
    "ipsae_best", "ipsae_mean", "ipsae_std",
    "lis_best", "lis_mean", "lis_std",
    "interaction_call",
    "template_bias_flag", "template_coverage_a", "template_coverage_b", "n_templates",
    "n_predictions", "n_predictions_with_pae",
]


def safe_name(target: str) -> str:
    return target.replace("/", "-").replace(" ", "_").replace(":", "_")


def safe_name_alt(target: str) -> str:
    """Alternative safe_name convention (Minerva uses / -> _)."""
    return target.replace("/", "_").replace(" ", "_").replace(":", "_")


def _read_scores_json(scores_path: Path) -> dict:
    """Read scores.json and return the summary dict.

    Handles both old flat schema (pre-v2) and new nested schema (v2+).
    """
    with open(scores_path) as f:
        data = json.load(f)

    if "summary" in data:
        return data["summary"]

    # Legacy flat schema -- map old field names
    return {
        "iptm_best": data.get("iptm_best"),
        "iptm_mean": data.get("iptm_mean"),
        "iptm_std": None,
        "ptm_best": data.get("ptm_best"),
        "interchain_pae_best": data.get("mean_interchain_pae"),
        "interchain_pae_mean": data.get("mean_interchain_pae"),
        "interchain_pae_std": None,
        "interface_plddt_best": data.get("mean_interface_plddt"),
        "interface_plddt_mean": data.get("mean_interface_plddt"),
        "n_interface_residues_best": data.get("n_interface_residues"),
        "ipsae_best": None,
        "ipsae_mean": None,
        "ipsae_std": None,
        "lis_best": None,
        "lis_mean": None,
        "lis_std": None,
        "interaction_call": data.get("interaction_call", "no_data"),
        "template_bias_flag": None,
        "template_coverage_a": None,
        "template_coverage_b": None,
        "n_templates": None,
        "n_predictions": data.get("n_models_parsed", 0),
        "n_predictions_with_pae": data.get("n_models_parsed", 0),
    }


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
    with open(candidates_csv, newline="") as f:
        candidates = list(csv.DictReader(f))

    uniprot_map = {c["target"]: c.get("uniprot", "") for c in candidates}

    fieldnames = ["target", "uniprot"] + SUMMARY_FIELDS

    rows = []
    n_found = 0
    n_missing = 0

    for cand in candidates:
        target = cand["target"]
        sname = safe_name(target)
        uniprot = uniprot_map.get(target, "")

        # Find scores.json
        target_dir = af_output_dir / sname
        if not target_dir.exists():
            target_dir = af_output_dir / safe_name_alt(target)
        nested = target_dir / f"sflt1_vs_{sname}"
        if not nested.exists():
            nested = target_dir / f"sflt1_vs_{safe_name_alt(target)}"
        actual = nested if nested.exists() else target_dir

        scores_path = actual / "scores.json" if actual.exists() else None

        if scores_path and scores_path.exists():
            try:
                summary = _read_scores_json(scores_path)
            except Exception as e:
                logger.warning("Failed to parse scores.json for %s: %s", target, e)
                summary = {}

            n_found += 1

            row = {"target": target, "uniprot": uniprot}
            for field in SUMMARY_FIELDS:
                v = summary.get(field)
                row[field] = v if v is not None else ""
            rows.append(row)
        else:
            n_missing += 1
            row = {"target": target, "uniprot": uniprot}
            for field in SUMMARY_FIELDS:
                row[field] = ""
            row["interaction_call"] = "not_run"
            rows.append(row)
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
