"""04_parse_results.py -- Parse AlphaFold Multimer output files.

Primary source: scores.json written by 09_extract_and_cleanup.py.
Fallback: direct parsing of ranking_debug.json and result pkl files.

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


def compute_interface_metrics(
    pkl_path: Path,
    chain_a_len: int,
    contact_threshold: float = 10.0,
) -> tuple[float, int]:
    """Compute interface pLDDT and interface residue count from AF2 pkl.

    Interface residues are defined as residues where the minimum inter-chain
    PAE to any residue in the other chain is below contact_threshold (Angstroms).

    Parameters
    ----------
    pkl_path : Path
        Path to result_model_*.pkl file.
    chain_a_len : int
        Number of residues in chain A (sFLT1).
    contact_threshold : float
        PAE cutoff (Angstroms) for defining interface residues. Default 10.0.

    Returns
    -------
    tuple[float, int]
        (mean_interface_plddt, n_interface_residues). Returns (nan, 0) on failure.
    """
    import pickle

    try:
        with open(pkl_path, "rb") as f:
            data = pickle.load(f)
    except Exception as e:
        logger.warning("Failed to load pkl %s: %s", pkl_path, e)
        return np.nan, 0

    plddt = data.get("plddt", None)   # shape [N_total]
    pae = data.get("predicted_aligned_error", None)  # shape [N, N]

    if plddt is None or pae is None:
        logger.debug("pkl %s missing plddt or predicted_aligned_error", pkl_path)
        return np.nan, 0

    pae = np.array(pae)
    plddt = np.array(plddt)
    n = pae.shape[0]

    if chain_a_len >= n or chain_a_len <= 0:
        logger.warning(
            "chain_a_len=%d out of range for pae shape %s in %s",
            chain_a_len, pae.shape, pkl_path,
        )
        return np.nan, 0

    # Inter-chain PAE blocks
    block_ab = pae[:chain_a_len, chain_a_len:]   # A->B: shape [chain_a_len, chain_b_len]
    block_ba = pae[chain_a_len:, :chain_a_len]   # B->A: shape [chain_b_len, chain_a_len]

    # Interface residues in chain A: min PAE to any B residue < threshold
    min_pae_a = np.min(block_ab, axis=1)   # shape [chain_a_len]
    interface_a = min_pae_a < contact_threshold

    # Interface residues in chain B: min PAE to any A residue < threshold
    min_pae_b = np.min(block_ba, axis=1)   # shape [chain_b_len]
    interface_b = min_pae_b < contact_threshold

    n_interface = int(np.sum(interface_a) + np.sum(interface_b))

    # Build full-length interface mask
    interface_mask = np.zeros(n, dtype=bool)
    interface_mask[:chain_a_len] = interface_a
    interface_mask[chain_a_len:] = interface_b

    if np.sum(interface_mask) == 0:
        return np.nan, 0

    mean_plddt = float(np.mean(plddt[interface_mask]))
    return round(mean_plddt, 2), n_interface


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


def _read_scores_json(scores_path: Path) -> dict | None:
    """Read scores.json and return the summary dict, handling both schemas."""
    try:
        with open(scores_path) as f:
            data = json.load(f)
    except Exception as e:
        logger.warning("Failed to parse %s: %s", scores_path, e)
        return None

    if "error" in data:
        return None

    if "summary" in data:
        return data["summary"]

    # Legacy flat schema
    return {
        "iptm_best": data.get("iptm_best"),
        "iptm_mean": data.get("iptm_mean"),
        "ptm_best": data.get("ptm_best"),
        "interchain_pae_best": data.get("mean_interchain_pae"),
        "interchain_pae_mean": data.get("mean_interchain_pae"),
        "interface_plddt_best": data.get("mean_interface_plddt"),
        "interface_plddt_mean": data.get("mean_interface_plddt"),
        "n_interface_residues_best": data.get("n_interface_residues"),
        "ipsae_best": None,
        "ipsae_mean": None,
        "lis_best": None,
        "lis_mean": None,
        "interaction_call": data.get("interaction_call", "no_data"),
        "template_bias_flag": None,
        "template_coverage_a": None,
        "template_coverage_b": None,
        "n_predictions": data.get("n_models_parsed", 0),
        "n_predictions_with_pae": data.get("n_models_parsed", 0),
    }


def _resolve_result_dir(af_output_dir: Path, target: str) -> Path | None:
    """Resolve the AF2 result directory for a target, handling naming variants."""
    safe = target.replace("/", "-").replace(" ", "_").replace(":", "_")
    safe_alt = target.replace("/", "_").replace(" ", "_").replace(":", "_")

    for sname in (safe, safe_alt):
        result_dir = af_output_dir / sname
        nested = result_dir / f"sflt1_vs_{sname}"
        if nested.exists():
            return nested
        if result_dir.exists():
            return result_dir

    return None


def _row_from_scores_json(summary: dict, target: str, uniprot: str) -> dict:
    """Build a DataFrame row from a scores.json summary dict."""
    def _get(key, default=np.nan):
        v = summary.get(key)
        return v if v is not None else default

    return {
        "target": target,
        "uniprot": uniprot,
        "iptm_best": _get("iptm_best"),
        "iptm_mean": _get("iptm_mean"),
        "iptm_std": _get("iptm_std"),
        "ptm_best": _get("ptm_best"),
        "mean_interchain_pae": _get("interchain_pae_mean"),
        "interchain_pae_best": _get("interchain_pae_best"),
        "interchain_pae_std": _get("interchain_pae_std"),
        "mean_interface_plddt": _get("interface_plddt_mean"),
        "interface_plddt_best": _get("interface_plddt_best"),
        "n_interface_residues": _get("n_interface_residues_best", 0),
        "ipsae_best": _get("ipsae_best"),
        "ipsae_mean": _get("ipsae_mean"),
        "lis_best": _get("lis_best"),
        "lis_mean": _get("lis_mean"),
        "interaction_call": summary.get("interaction_call", "no_data"),
        "template_bias_flag": _get("template_bias_flag", False),
        "template_coverage_a": _get("template_coverage_a"),
        "template_coverage_b": _get("template_coverage_b"),
        "n_models_parsed": _get("n_predictions", 0),
    }


def _empty_row(target: str, uniprot: str, call: str = "not_run") -> dict:
    """Build an empty DataFrame row for missing/failed targets."""
    return {
        "target": target, "uniprot": uniprot,
        "iptm_best": np.nan, "iptm_mean": np.nan, "iptm_std": np.nan,
        "ptm_best": np.nan,
        "mean_interchain_pae": np.nan, "interchain_pae_best": np.nan,
        "interchain_pae_std": np.nan,
        "mean_interface_plddt": np.nan, "interface_plddt_best": np.nan,
        "n_interface_residues": 0,
        "ipsae_best": np.nan, "ipsae_mean": np.nan,
        "lis_best": np.nan, "lis_mean": np.nan,
        "interaction_call": call,
        "template_bias_flag": False,
        "template_coverage_a": np.nan, "template_coverage_b": np.nan,
        "n_models_parsed": 0,
    }


def parse_all_results(
    candidates_path: Path,
    af_output_dir: Path,
    sflt1_length: int = 304,
) -> pd.DataFrame:
    """Parse AlphaFold results for all candidates.

    Primary source: scores.json (written by 09_extract_and_cleanup.py).
    Fallback: direct parsing of ranking_debug.json + pkl files.

    Parameters
    ----------
    candidates_path : Path
        Path to step03_candidates.csv.
    af_output_dir : Path
        Root directory containing per-candidate AlphaFold output subdirs.
    sflt1_length : int
        Length of sFLT1 chain A. Corrected constructs: D1-D3=304, D1-D6=631, D1-D7=721.

    Returns
    -------
    pd.DataFrame
        Interaction scores with all metrics from scores.json plus ranking columns.
    """
    candidates = pd.read_csv(candidates_path)
    rows = []
    n_from_json = 0
    n_from_pkl = 0

    for _, cand in candidates.iterrows():
        target = cand["target"]
        uniprot = cand["uniprot"]

        actual_result_dir = _resolve_result_dir(af_output_dir, target)

        if actual_result_dir is None:
            logger.warning("No output directory for %s", target)
            rows.append(_empty_row(target, uniprot, "not_run"))
            continue

        # --- Primary path: scores.json ---
        scores_path = actual_result_dir / "scores.json"
        if scores_path.exists():
            summary = _read_scores_json(scores_path)
            if summary is not None:
                rows.append(_row_from_scores_json(summary, target, uniprot))
                n_from_json += 1
                continue

        # --- Fallback: direct pkl parsing ---
        ranking_path = actual_result_dir / "ranking_debug.json"
        if not ranking_path.exists():
            logger.warning("No ranking_debug.json for %s", target)
            rows.append(_empty_row(target, uniprot, "failed"))
            continue

        model_scores = parse_ranking_debug(ranking_path)
        iptm_values = [s["iptm"] for s in model_scores.values() if not np.isnan(s["iptm"])]
        ptm_values = [s["ptm"] for s in model_scores.values() if not np.isnan(s["ptm"])]

        iptm_best = max(iptm_values) if iptm_values else np.nan
        iptm_mean = float(np.mean(iptm_values)) if iptm_values else np.nan
        ptm_best = max(ptm_values) if ptm_values else np.nan

        mean_pae = np.nan
        mean_interface_plddt = np.nan
        n_interface_residues = 0

        if model_scores:
            best_model = min(model_scores, key=lambda k: model_scores[k]["order"])

            pae_candidates = [
                actual_result_dir / f"pae_{best_model}.json",
                actual_result_dir / f"result_{best_model}.pkl",
            ]
            for pae_path in pae_candidates:
                if pae_path.exists():
                    pae_mat = load_pae_matrix(pae_path)
                    if pae_mat is not None:
                        mean_pae = compute_interchain_pae(pae_mat, sflt1_length)
                        break

            best_pkl = actual_result_dir / f"result_{best_model}.pkl"
            if best_pkl.exists():
                mean_interface_plddt, n_interface_residues = compute_interface_metrics(
                    best_pkl, sflt1_length
                )

        call = classify_interaction(iptm_best, mean_pae)

        row = _empty_row(target, uniprot, call)
        row.update({
            "iptm_best": round(iptm_best, 4) if not np.isnan(iptm_best) else np.nan,
            "iptm_mean": round(iptm_mean, 4) if not np.isnan(iptm_mean) else np.nan,
            "ptm_best": round(ptm_best, 4) if not np.isnan(ptm_best) else np.nan,
            "mean_interchain_pae": round(mean_pae, 2) if not np.isnan(mean_pae) else np.nan,
            "mean_interface_plddt": mean_interface_plddt,
            "n_interface_residues": n_interface_residues,
            "n_models_parsed": len(iptm_values),
        })
        rows.append(row)
        n_from_pkl += 1

    logger.info("Parsed %d targets: %d from scores.json, %d from pkl fallback",
                len(rows), n_from_json, n_from_pkl)

    df = pd.DataFrame(rows)

    # Rank by ipTM
    df["rank_by_iptm"] = df["iptm_best"].rank(ascending=False, method="min").astype("Int64")

    return df
