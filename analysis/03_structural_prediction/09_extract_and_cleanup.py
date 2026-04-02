#!/usr/bin/env python3
"""09_extract_and_cleanup.py -- Per-job post-AF2 score extraction and pkl cleanup.

Called automatically at the end of each AlphaFold job. Extracts scores from
ALL result pkl files into a lightweight scores.json sidecar, then deletes
result pkl files to reclaim disk space. features.pkl is retained.

Metrics extracted per prediction:
  - interchain_pae: mean PAE across inter-chain residue pairs
  - interface_plddt: mean pLDDT of interface residues (PAE < 10A to other chain)
  - n_interface_residues: count of interface residues
  - ipsae: PAE-filtered ipTM -- fraction of inter-chain residue pairs with PAE < 5A,
           weighted by pLDDT. Rescues candidates with small confident interfaces
           that raw ipTM averages away.
  - lis: Local Interaction Score -- mean pLDDT of the top-K (K=20) inter-chain
         residue pairs ranked by PAE. Captures local interface quality independent
         of complex size.

Summary also includes:
  - template_coverage_a/b: fraction of chain residues covered by at least one
    PDB template (from features.pkl). Low coverage on the partner chain flags
    template bias -- sFLT1 Ig-like domains have abundant templates, so poorly-
    templated partners may get artificially low scores.

Designed for HPC: stdlib + numpy only (no pandas dependency).

Usage (called by AF2 job wrapper, not manually):
    python3 09_extract_and_cleanup.py <af_subdir> <sflt1_length>

    af_subdir:     Path to the AF2 output subdirectory containing ranking_debug.json
    sflt1_length:  Chain A (sFLT1) residue count (304, 631, or 721)

Outputs:
    <af_subdir>/scores.json    -- extracted metrics (per-prediction + summary)
    Deletes: result_model_*.pkl (keeps PDBs, ranking_debug.json, scores.json, features.pkl)
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

# ipSAE parameters
IPSAE_PAE_CUTOFF = 5.0  # Angstroms -- pairs below this are "confident contacts"

# LIS parameters
LIS_TOP_K = 20  # number of best inter-chain pairs to average


def _extract_from_pkl(pkl_path: Path, sflt1_length: int) -> dict | None:
    """Extract all per-prediction metrics from one result pkl.

    Returns dict with interchain_pae, interface_plddt, n_interface_residues,
    ipsae, lis -- or None on failure.
    """
    try:
        with open(pkl_path, "rb") as f:
            data = pickle.load(f)
    except Exception as e:
        logger.warning("Failed to load %s: %s", pkl_path.name, e)
        return None

    pae = data.get("predicted_aligned_error")
    if pae is None:
        return None

    pae = np.array(pae)
    n = pae.shape[0]
    if not (0 < sflt1_length < n):
        return None

    plddt = data.get("plddt")
    plddt = np.array(plddt) if plddt is not None else None

    block_ab = pae[:sflt1_length, sflt1_length:]  # A->B
    block_ba = pae[sflt1_length:, :sflt1_length]  # B->A

    # -- Mean interchain PAE --
    interchain_vals = np.concatenate([block_ab.ravel(), block_ba.ravel()])
    interchain_pae = float(np.mean(interchain_vals))

    # -- Interface residues (PAE < 10A) --
    threshold = 10.0
    min_pae_a = np.min(block_ab, axis=1)
    min_pae_b = np.min(block_ba, axis=1)
    interface_a = min_pae_a < threshold
    interface_b = min_pae_b < threshold
    n_interface = int(np.sum(interface_a) + np.sum(interface_b))

    # Interface pLDDT
    interface_plddt = None
    if plddt is not None and n_interface > 0:
        mask = np.zeros(n, dtype=bool)
        mask[:sflt1_length] = interface_a
        mask[sflt1_length:] = interface_b
        interface_plddt = round(float(np.mean(plddt[mask])), 2)

    # -- ipSAE: PAE-filtered interaction score --
    # Fraction of inter-chain pairs with PAE < cutoff, weighted by mean pLDDT
    # of those pairs. Rescues small confident interfaces.
    ipsae = None
    if plddt is not None:
        chain_b_len = n - sflt1_length
        n_pairs = sflt1_length * chain_b_len * 2  # both AB and BA blocks
        confident_ab = block_ab < IPSAE_PAE_CUTOFF
        confident_ba = block_ba < IPSAE_PAE_CUTOFF
        n_confident = int(np.sum(confident_ab) + np.sum(confident_ba))
        if n_confident > 0:
            # pLDDT for residues involved in confident contacts
            plddt_a = plddt[:sflt1_length]
            plddt_b = plddt[sflt1_length:]
            # For each confident pair (i,j), score = (plddt_i + plddt_j) / 200
            # Efficient: mean pLDDT of residues with any confident contact
            has_conf_a = np.any(confident_ab, axis=1) | np.any(confident_ba.T, axis=1)
            has_conf_b = np.any(confident_ba, axis=1) | np.any(confident_ab.T, axis=1)
            conf_plddt = np.concatenate([
                plddt_a[has_conf_a], plddt_b[has_conf_b]
            ])
            mean_conf_plddt = float(np.mean(conf_plddt)) / 100.0  # normalize to 0-1
            frac_confident = n_confident / n_pairs
            ipsae = round(float(frac_confident * mean_conf_plddt), 4)
        else:
            ipsae = 0.0

    # -- LIS: Local Interaction Score --
    # Mean pLDDT of residue pairs at the top-K lowest PAE inter-chain contacts.
    # Captures local interface quality independent of complex size.
    lis = None
    if plddt is not None:
        # Build pair-level PAE + pLDDT for inter-chain contacts
        # Use block_ab (A predicting B) -- rows=A residues, cols=B residues
        plddt_a = plddt[:sflt1_length]
        plddt_b = plddt[sflt1_length:]

        # Flatten block_ab, get top-K lowest PAE pairs
        flat_pae = block_ab.ravel()
        k = min(LIS_TOP_K, len(flat_pae))
        top_k_idx = np.argpartition(flat_pae, k)[:k]
        top_k_rows = top_k_idx // block_ab.shape[1]
        top_k_cols = top_k_idx % block_ab.shape[1]

        pair_plddt = (plddt_a[top_k_rows] + plddt_b[top_k_cols]) / 2.0
        lis = round(float(np.mean(pair_plddt)), 2)

    return {
        "interchain_pae": round(interchain_pae, 2),
        "interface_plddt": interface_plddt,
        "n_interface_residues": n_interface,
        "ipsae": ipsae,
        "lis": lis,
    }


def _extract_template_coverage(af_subdir: Path, sflt1_length: int) -> dict:
    """Extract template coverage per chain from features.pkl.

    Template coverage = fraction of residues with at least one template atom
    resolved (template_all_atom_mask > 0 for any template).

    Low coverage on the partner chain flags template bias.
    """
    features_path = af_subdir / "features.pkl"
    if not features_path.exists():
        return {"template_coverage_a": None, "template_coverage_b": None,
                "n_templates": None}

    try:
        with open(features_path, "rb") as f:
            feat = pickle.load(f)
    except Exception as e:
        logger.warning("Failed to load features.pkl: %s", e)
        return {"template_coverage_a": None, "template_coverage_b": None,
                "n_templates": None}

    n_templates = int(feat.get("num_templates", 0))

    tmpl_mask = feat.get("template_all_atom_mask")  # shape [n_templates, n_res, 37]
    if tmpl_mask is None or n_templates == 0:
        return {"template_coverage_a": 0.0, "template_coverage_b": 0.0,
                "n_templates": n_templates}

    tmpl_mask = np.array(tmpl_mask)
    # Residue has template if any template has any atom resolved
    has_template = np.any(tmpl_mask > 0, axis=(0, 2))  # shape [n_res]
    n_res = len(has_template)

    if sflt1_length >= n_res or sflt1_length <= 0:
        return {"template_coverage_a": None, "template_coverage_b": None,
                "n_templates": n_templates}

    cov_a = float(np.mean(has_template[:sflt1_length]))
    cov_b = float(np.mean(has_template[sflt1_length:]))

    return {
        "template_coverage_a": round(cov_a, 3),
        "template_coverage_b": round(cov_b, 3),
        "n_templates": n_templates,
    }


def extract_scores(af_subdir: Path, sflt1_length: int) -> dict:
    """Extract scores from ranking_debug.json and ALL result pkl files.

    Returns a dict with:
    - per_prediction: list of dicts, one per prediction (ranked by ipTM+pTM)
    - summary: aggregate statistics across all predictions
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

    # Extract per-prediction scores
    per_prediction = []
    for model_name in order:
        iptm = iptm_scores.get(model_name)
        ptm = ranking.get("ptm", {}).get(model_name)

        entry = {
            "model": model_name,
            "iptm": round(float(iptm), 4) if iptm is not None and not np.isnan(iptm) else None,
            "ptm": round(float(ptm), 4) if ptm is not None and not np.isnan(ptm) else None,
        }

        null_scores = {
            "interchain_pae": None, "interface_plddt": None,
            "n_interface_residues": None, "ipsae": None, "lis": None,
        }

        pkl_path = af_subdir / f"result_{model_name}.pkl"
        if pkl_path.exists():
            pkl_scores = _extract_from_pkl(pkl_path, sflt1_length)
            entry.update(pkl_scores if pkl_scores else null_scores)
        else:
            logger.warning("Missing pkl for %s", model_name)
            entry.update(null_scores)

        per_prediction.append(entry)

    # Template coverage (from features.pkl, shared across predictions)
    template_info = _extract_template_coverage(af_subdir, sflt1_length)

    # Compute summary statistics
    def _vals(key):
        return [p[key] for p in per_prediction if p[key] is not None]

    iptm_vals = _vals("iptm")
    ptm_vals = _vals("ptm")
    pae_vals = _vals("interchain_pae")
    plddt_vals = _vals("interface_plddt")
    iface_vals = _vals("n_interface_residues")
    ipsae_vals = _vals("ipsae")
    lis_vals = _vals("lis")

    if not iptm_vals:
        return {"error": "no iptm scores", "n_models_parsed": 0}

    def _stat(vals, best_fn=max):
        """Return best, mean, std for a list of values."""
        if not vals:
            return None, None, None
        best = round(best_fn(vals), 4)
        mean = round(float(np.mean(vals)), 4)
        std = round(float(np.std(vals)), 4) if len(vals) > 1 else None
        return best, mean, std

    iptm_best, iptm_mean, iptm_std = _stat(iptm_vals, max)
    ptm_best, _, _ = _stat(ptm_vals, max)
    pae_best, pae_mean, pae_std = _stat(pae_vals, min)
    plddt_best, plddt_mean, _ = _stat(plddt_vals, max)
    ipsae_best, ipsae_mean, ipsae_std = _stat(ipsae_vals, max)
    lis_best, lis_mean, lis_std = _stat(lis_vals, max)

    # Classify using best scores
    if pae_best is not None:
        if iptm_best > IPTM_HIGH_CONF and pae_best < PAE_HIGH_CONF:
            call = "high_confidence"
        elif iptm_best > IPTM_THRESHOLD and pae_best < PAE_THRESHOLD:
            call = "predicted"
        else:
            call = "low_confidence"
    else:
        call = "no_data"

    # Template bias flag: partner chain coverage < 30% suggests scores may
    # be artifactually low due to poor template availability, not biology.
    cov_b = template_info.get("template_coverage_b")
    template_bias_flag = cov_b is not None and cov_b < 0.3

    summary = {
        "iptm_best": iptm_best,
        "iptm_mean": iptm_mean,
        "iptm_std": iptm_std,
        "ptm_best": ptm_best,
        "interchain_pae_best": pae_best,
        "interchain_pae_mean": pae_mean,
        "interchain_pae_std": pae_std,
        "interface_plddt_best": plddt_best,
        "interface_plddt_mean": plddt_mean,
        "n_interface_residues_best": max(iface_vals) if iface_vals else None,
        "ipsae_best": ipsae_best,
        "ipsae_mean": ipsae_mean,
        "ipsae_std": ipsae_std,
        "lis_best": lis_best,
        "lis_mean": lis_mean,
        "lis_std": lis_std,
        "interaction_call": call,
        "template_bias_flag": template_bias_flag,
        **template_info,
        "n_predictions": len(per_prediction),
        "n_predictions_with_pae": len(pae_vals),
        "sflt1_length": sflt1_length,
    }

    return {
        "summary": summary,
        "per_prediction": per_prediction,
    }


def cleanup_pkl(af_subdir: Path) -> tuple[int, float]:
    """Delete result_model_*.pkl files only. Retains features.pkl.

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

    # Extract scores from ALL predictions
    scores = extract_scores(af_subdir, sflt1_length)
    scores_path = af_subdir / "scores.json"
    with open(scores_path, "w") as f:
        json.dump(scores, f, indent=2)

    s = scores.get("summary", scores)
    logger.info(
        "Saved: %s (ipTM=%.4f, PAE=%.1f%s, ipSAE=%.4f, LIS=%.1f, call=%s, "
        "tmpl_cov_b=%.1f%%, n=%d)",
        scores_path,
        s.get("iptm_best", 0),
        s.get("interchain_pae_mean", 0) or 0,
        "±%.1f" % s["interchain_pae_std"] if s.get("interchain_pae_std") else "",
        s.get("ipsae_best", 0) or 0,
        s.get("lis_best", 0) or 0,
        s.get("interaction_call", "?"),
        (s.get("template_coverage_b", 0) or 0) * 100,
        s.get("n_predictions_with_pae", 0),
    )

    if s.get("template_bias_flag"):
        logger.warning(
            "TEMPLATE BIAS: partner chain template coverage %.0f%% -- "
            "low scores may reflect poor template availability, not biology",
            (s.get("template_coverage_b", 0) or 0) * 100,
        )

    # Cleanup result pkl files (retain features.pkl)
    n_deleted, bytes_freed = cleanup_pkl(af_subdir)
    logger.info("Deleted %d result pkl files (%.1f GB freed)", n_deleted, bytes_freed / 1073741824)


if __name__ == "__main__":
    main()
