#!/usr/bin/env python3
"""09_extract_and_cleanup.py -- Per-job post-AF2 score extraction and pkl cleanup.

Called automatically at the end of each AlphaFold job. Extracts scores from
result pkl files into a lightweight scores.json sidecar, then deletes the
pkl files to reclaim disk space.

Designed for HPC: stdlib + numpy only (no pandas dependency).

Usage (called by AF2 job wrapper, not manually):
    python3 09_extract_and_cleanup.py <af_subdir> <sflt1_length>

    af_subdir:     Path to the AF2 output subdirectory containing ranking_debug.json
    sflt1_length:  Chain A (sFLT1) residue count (304, 631, or 721)

Outputs:
    <af_subdir>/scores.json    -- extracted metrics
    Deletes: result_model_*.pkl (keeps features.pkl, PDBs, ranking_debug.json)
"""

import json
import logging
import pickle
import sys
from pathlib import Path

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# Interaction thresholds (match 04_parse_results.py)
IPTM_THRESHOLD = 0.6
PAE_THRESHOLD = 15.0
IPTM_HIGH_CONF = 0.8
PAE_HIGH_CONF = 10.0


def extract_scores(af_subdir: Path, sflt1_length: int) -> dict:
    """Extract all scores from AF2 result pkl files.

    Parameters
    ----------
    af_subdir : Path
        Directory containing ranking_debug.json and result_model_*.pkl.
    sflt1_length : int
        Chain A (sFLT1) residue count.

    Returns
    -------
    dict
        Extracted metrics ready for JSON serialization.
    """
    ranking_path = af_subdir / "ranking_debug.json"
    if not ranking_path.exists():
        return {"error": "no ranking_debug.json", "n_models_parsed": 0}

    with open(ranking_path) as f:
        ranking = json.load(f)

    order = ranking.get("order", [])
    iptm_scores = ranking.get("iptm+ptm", ranking.get("iptm", {}))
    if not isinstance(iptm_scores, dict):
        iptm_scores = {}

    # Collect per-model scores from ranking_debug.json
    iptm_values = []
    ptm_values = []
    for model_name in order:
        iptm = iptm_scores.get(model_name)
        if iptm is not None and not np.isnan(iptm):
            iptm_values.append(float(iptm))
        ptm = ranking.get("ptm", {}).get(model_name)
        if ptm is not None and not np.isnan(ptm):
            ptm_values.append(float(ptm))

    if not iptm_values:
        return {"error": "no iptm scores", "n_models_parsed": 0}

    iptm_best = max(iptm_values)
    iptm_mean = float(np.mean(iptm_values))
    ptm_best = max(ptm_values) if ptm_values else None

    # Load best model pkl for PAE and pLDDT
    best_model = order[0] if order else None
    mean_interchain_pae = None
    mean_interface_plddt = None
    n_interface_residues = 0

    if best_model:
        best_pkl = af_subdir / f"result_{best_model}.pkl"
        if best_pkl.exists():
            try:
                with open(best_pkl, "rb") as f:
                    data = pickle.load(f)

                pae = data.get("predicted_aligned_error")
                if pae is not None:
                    pae = np.array(pae)
                    n = pae.shape[0]
                    if 0 < sflt1_length < n:
                        block_ab = pae[:sflt1_length, sflt1_length:]
                        block_ba = pae[sflt1_length:, :sflt1_length]
                        mean_interchain_pae = round(
                            float(np.mean(np.concatenate([block_ab.ravel(), block_ba.ravel()]))), 2
                        )

                        # Interface residues
                        threshold = 10.0
                        min_pae_a = np.min(block_ab, axis=1)
                        min_pae_b = np.min(block_ba, axis=1)
                        interface_a = min_pae_a < threshold
                        interface_b = min_pae_b < threshold
                        n_interface_residues = int(np.sum(interface_a) + np.sum(interface_b))

                        # Interface pLDDT
                        plddt = data.get("plddt")
                        if plddt is not None and n_interface_residues > 0:
                            plddt = np.array(plddt)
                            mask = np.zeros(n, dtype=bool)
                            mask[:sflt1_length] = interface_a
                            mask[sflt1_length:] = interface_b
                            mean_interface_plddt = round(float(np.mean(plddt[mask])), 2)

            except Exception as e:
                logger.warning("Failed to extract from %s: %s", best_pkl.name, e)

    # Classify
    if mean_interchain_pae is not None:
        if iptm_best > IPTM_HIGH_CONF and mean_interchain_pae < PAE_HIGH_CONF:
            call = "high_confidence"
        elif iptm_best > IPTM_THRESHOLD and mean_interchain_pae < PAE_THRESHOLD:
            call = "predicted"
        else:
            call = "low_confidence"
    else:
        call = "no_data"

    return {
        "iptm_best": round(iptm_best, 4),
        "iptm_mean": round(iptm_mean, 4),
        "ptm_best": round(ptm_best, 4) if ptm_best is not None else None,
        "mean_interchain_pae": mean_interchain_pae,
        "mean_interface_plddt": mean_interface_plddt,
        "interaction_call": call,
        "n_interface_residues": n_interface_residues,
        "n_models_parsed": len(iptm_values),
        "sflt1_length": sflt1_length,
    }


def cleanup_pkl(af_subdir: Path) -> tuple[int, float]:
    """Delete result_model_*.pkl files. Keep features.pkl and everything else.

    Returns
    -------
    tuple[int, float]
        (n_files_deleted, bytes_freed)
    """
    pkls = sorted(af_subdir.glob("result_model_*.pkl"))
    total_bytes = 0
    for p in pkls:
        total_bytes += p.stat().st_size
        p.unlink()
    return len(pkls), total_bytes


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <af_subdir> <sflt1_length>", file=sys.stderr)
        sys.exit(1)

    af_subdir = Path(sys.argv[1])
    sflt1_length = int(sys.argv[2])

    if not af_subdir.exists():
        logger.error("Directory not found: %s", af_subdir)
        sys.exit(1)

    # Extract scores
    scores = extract_scores(af_subdir, sflt1_length)
    scores_path = af_subdir / "scores.json"
    with open(scores_path, "w") as f:
        json.dump(scores, f, indent=2)
    logger.info("Saved: %s (ipTM=%.4f, call=%s)",
                scores_path, scores.get("iptm_best", 0), scores.get("interaction_call", "?"))

    # Cleanup pkl files
    n_deleted, bytes_freed = cleanup_pkl(af_subdir)
    logger.info("Deleted %d pkl files (%.1f MB freed)", n_deleted, bytes_freed / 1048576)


if __name__ == "__main__":
    main()
