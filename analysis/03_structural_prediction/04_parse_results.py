"""04_parse_results.py -- Parse AlphaFold Multimer output files.

Extracts ipTM, pTM, PAE matrices from ranking_debug.json and model outputs.
Computes inter-chain PAE and interface residues.

Scoring thresholds:
  - ipTM > 0.6 AND mean inter-chain PAE < 15A = predicted interaction
  - ipTM > 0.8 AND PAE < 10A = high-confidence predicted interaction
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Interaction call thresholds
IPTM_THRESHOLD = 0.6
PAE_THRESHOLD = 15.0
IPTM_HIGH_CONF = 0.8
PAE_HIGH_CONF = 10.0

# Interface contact distance threshold (Angstroms, CA-CA)
CONTACT_DISTANCE = 8.0


def parse_ranking_debug(ranking_path: Path) -> dict[str, dict]:
    """Parse ranking_debug.json for ipTM and pTM scores.

    Parameters
    ----------
    ranking_path : Path
        Path to ranking_debug.json.

    Returns
    -------
    dict[str, dict]
        Mapping of model name -> {"iptm": float, "ptm": float, "order": int}.
    """
    with open(ranking_path) as f:
        data = json.load(f)

    results = {}
    order = data.get("order", [])

    iptm_scores = data.get("iptm+ptm", data.get("iptm", {}))
    if not isinstance(iptm_scores, dict):
        iptm_scores = {}

    for idx, model_name in enumerate(order):
        results[model_name] = {
            "iptm": iptm_scores.get(model_name, np.nan),
            "ptm": data.get("ptm", {}).get(model_name, np.nan)
                   if "ptm" in data else np.nan,
            "order": idx,
        }

    return results


def compute_interchain_pae(
    pae_matrix: np.ndarray,
    chain_a_len: int,
) -> float:
    """Compute mean inter-chain PAE from a full PAE matrix.

    The PAE matrix is (N_total x N_total) where N_total = chain_a_len + chain_b_len.
    Inter-chain blocks are:
      - Top-right: rows [0:chain_a_len], cols [chain_a_len:]
      - Bottom-left: rows [chain_a_len:], cols [0:chain_a_len]

    Parameters
    ----------
    pae_matrix : np.ndarray
        2D PAE matrix of shape (N_total, N_total).
    chain_a_len : int
        Number of residues in chain A (sFLT1).

    Returns
    -------
    float
        Mean inter-chain PAE in Angstroms.
    """
    n = pae_matrix.shape[0]
    if chain_a_len >= n or chain_a_len <= 0:
        return np.nan

    # Extract inter-chain blocks
    block_ab = pae_matrix[:chain_a_len, chain_a_len:]
    block_ba = pae_matrix[chain_a_len:, :chain_a_len]

    interchain_values = np.concatenate([block_ab.ravel(), block_ba.ravel()])
    return float(np.mean(interchain_values))


def load_pae_matrix(pae_path: Path) -> np.ndarray | None:
    """Load a PAE matrix from a JSON or pkl file.

    AlphaFold outputs PAE as either a JSON array or a pickle.

    Parameters
    ----------
    pae_path : Path
        Path to PAE file.

    Returns
    -------
    np.ndarray or None
        2D PAE matrix, or None on failure.
    """
    if pae_path.suffix == ".json":
        try:
            with open(pae_path) as f:
                data = json.load(f)
            if isinstance(data, list) and isinstance(data[0], dict):
                # Format: [{"predicted_aligned_error": [[...]]}]
                pae = np.array(data[0]["predicted_aligned_error"])
            elif isinstance(data, list) and isinstance(data[0], list):
                pae = np.array(data)
            else:
                pae = np.array(data)
            return pae
        except Exception as e:
            logger.warning("Failed to parse PAE JSON %s: %s", pae_path, e)
            return None
    elif pae_path.suffix == ".pkl":
        import pickle
        try:
            with open(pae_path, "rb") as f:
                data = pickle.load(f)
            if isinstance(data, dict) and "predicted_aligned_error" in data:
                return np.array(data["predicted_aligned_error"])
            return None
        except Exception as e:
            logger.warning("Failed to parse PAE pkl %s: %s", pae_path, e)
            return None
    return None


def classify_interaction(iptm: float, mean_pae: float) -> str:
    """Classify predicted interaction based on ipTM and PAE thresholds.

    Parameters
    ----------
    iptm : float
        Interface pTM score.
    mean_pae : float
        Mean inter-chain predicted aligned error.

    Returns
    -------
    str
        "high_confidence", "predicted", or "low_confidence".
    """
    if np.isnan(iptm) or np.isnan(mean_pae):
        return "no_data"
    if iptm > IPTM_HIGH_CONF and mean_pae < PAE_HIGH_CONF:
        return "high_confidence"
    if iptm > IPTM_THRESHOLD and mean_pae < PAE_THRESHOLD:
        return "predicted"
    return "low_confidence"


def parse_all_results(
    candidates_path: Path,
    af_output_dir: Path,
    sflt1_length: int = 338,
) -> pd.DataFrame:
    """Parse AlphaFold results for all candidates.

    Parameters
    ----------
    candidates_path : Path
        Path to step03_candidates.csv.
    af_output_dir : Path
        Root directory containing per-candidate AlphaFold output subdirs.
    sflt1_length : int
        Length of sFLT1 chain (D1-D3 = 338).

    Returns
    -------
    pd.DataFrame
        Interaction scores with columns: target, uniprot, iptm_best,
        iptm_mean, ptm_best, mean_interchain_pae, interaction_call,
        n_models_parsed, rank_by_iptm.
    """
    candidates = pd.read_csv(candidates_path)
    rows = []

    for _, cand in candidates.iterrows():
        target = cand["target"]
        uniprot = cand["uniprot"]
        safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")

        result_dir = af_output_dir / safe_name
        if not result_dir.exists():
            logger.warning("No output directory for %s at %s", target, result_dir)
            rows.append({
                "target": target, "uniprot": uniprot,
                "iptm_best": np.nan, "iptm_mean": np.nan, "ptm_best": np.nan,
                "mean_interchain_pae": np.nan,
                "mean_interface_plddt": np.nan,
                "interaction_call": "not_run",
                "n_interface_residues": 0,
                "n_models_parsed": 0,
            })
            continue

        # Parse ranking_debug.json
        ranking_path = result_dir / "ranking_debug.json"
        if not ranking_path.exists():
            logger.warning("No ranking_debug.json for %s", target)
            rows.append({
                "target": target, "uniprot": uniprot,
                "iptm_best": np.nan, "iptm_mean": np.nan, "ptm_best": np.nan,
                "mean_interchain_pae": np.nan,
                "mean_interface_plddt": np.nan,
                "interaction_call": "failed",
                "n_interface_residues": 0,
                "n_models_parsed": 0,
            })
            continue

        model_scores = parse_ranking_debug(ranking_path)
        iptm_values = [s["iptm"] for s in model_scores.values() if not np.isnan(s["iptm"])]
        ptm_values = [s["ptm"] for s in model_scores.values() if not np.isnan(s["ptm"])]

        iptm_best = max(iptm_values) if iptm_values else np.nan
        iptm_mean = float(np.mean(iptm_values)) if iptm_values else np.nan
        ptm_best = max(ptm_values) if ptm_values else np.nan

        # Try to parse PAE for the best model
        mean_pae = np.nan
        if model_scores:
            best_model = min(model_scores, key=lambda k: model_scores[k]["order"])
            pae_candidates = [
                result_dir / f"pae_{best_model}.json",
                result_dir / f"result_{best_model}.pkl",
            ]
            for pae_path in pae_candidates:
                if pae_path.exists():
                    pae_mat = load_pae_matrix(pae_path)
                    if pae_mat is not None:
                        mean_pae = compute_interchain_pae(pae_mat, sflt1_length)
                        break

        call = classify_interaction(iptm_best, mean_pae)

        rows.append({
            "target": target,
            "uniprot": uniprot,
            "iptm_best": round(iptm_best, 4) if not np.isnan(iptm_best) else np.nan,
            "iptm_mean": round(iptm_mean, 4) if not np.isnan(iptm_mean) else np.nan,
            "ptm_best": round(ptm_best, 4) if not np.isnan(ptm_best) else np.nan,
            "mean_interchain_pae": round(mean_pae, 2) if not np.isnan(mean_pae) else np.nan,
            "mean_interface_plddt": np.nan,  # Requires structure parsing
            "interaction_call": call,
            "n_interface_residues": 0,  # Requires structure parsing
            "n_models_parsed": len(iptm_values),
        })

    df = pd.DataFrame(rows)

    # Rank by ipTM
    df["rank_by_iptm"] = df["iptm_best"].rank(ascending=False, method="min").astype("Int64")

    return df
