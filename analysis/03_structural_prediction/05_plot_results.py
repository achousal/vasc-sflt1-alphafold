"""05_plot_results.py -- Visualization for AlphaFold Multimer results.

Generates:
  - PAE heatmaps per candidate (inter-chain blocks highlighted)
  - ipTM bar chart (VEGFA positive control marked)
  - Summary table
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

# Plot style
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "figure.dpi": 150,
})


def plot_pae_heatmap(
    pae_matrix: np.ndarray,
    chain_a_len: int,
    target_name: str,
    output_path: Path,
) -> None:
    """Plot a PAE heatmap with inter-chain blocks highlighted.

    Parameters
    ----------
    pae_matrix : np.ndarray
        2D PAE matrix.
    chain_a_len : int
        Length of chain A (sFLT1).
    target_name : str
        Name of the partner protein.
    output_path : Path
        Output PDF path.
    """
    fig, ax = plt.subplots(figsize=(8, 7))

    im = ax.imshow(pae_matrix, cmap="Greens_r", vmin=0, vmax=30, aspect="equal")
    plt.colorbar(im, ax=ax, label="Predicted Aligned Error (A)", shrink=0.8)

    # Draw chain boundary lines
    n = pae_matrix.shape[0]
    ax.axhline(y=chain_a_len - 0.5, color="red", linewidth=1.5, linestyle="--")
    ax.axvline(x=chain_a_len - 0.5, color="red", linewidth=1.5, linestyle="--")

    ax.set_xlabel("Scored residue")
    ax.set_ylabel("Aligned residue")
    ax.set_title(f"PAE: sFLT1 D1-D3 vs {target_name}")

    # Label chain regions
    ax.text(chain_a_len / 2, -n * 0.03, "sFLT1", ha="center", fontsize=9, color="red")
    ax.text(
        chain_a_len + (n - chain_a_len) / 2, -n * 0.03,
        target_name, ha="center", fontsize=9, color="red",
    )

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved: %s", output_path.name)


def plot_iptm_barplot(
    scores_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot ipTM bar chart with VEGFA positive control marked.

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with columns: target, iptm_best, interaction_call.
    output_path : Path
        Output PDF path.
    """
    df = scores_df.dropna(subset=["iptm_best"]).sort_values("iptm_best", ascending=True)

    if df.empty:
        logger.warning("No ipTM scores to plot")
        return

    # Color by interaction call
    color_map = {
        "high_confidence": "#1B9E77",
        "predicted": "#D95F02",
        "low_confidence": "#7570B3",
        "not_run": "#999999",
        "failed": "#E7298A",
        "no_data": "#CCCCCC",
    }
    colors = [color_map.get(c, "#666666") for c in df["interaction_call"]]

    fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.35)))
    bars = ax.barh(range(len(df)), df["iptm_best"], color=colors, edgecolor="white")

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["target"])
    ax.set_xlabel("ipTM (best model)")
    ax.set_title("AlphaFold Multimer ipTM: sFLT1 D1-D3 vs candidates")

    # Threshold lines
    ax.axvline(x=0.6, color="gray", linestyle="--", linewidth=0.8, label="ipTM = 0.6")
    ax.axvline(x=0.8, color="gray", linestyle=":", linewidth=0.8, label="ipTM = 0.8")

    # Mark VEGFA positive control
    vegfa_idx = df.index[df["target"] == "VEGFA"].tolist()
    if vegfa_idx:
        for idx in vegfa_idx:
            pos = list(df.index).index(idx)
            ax.annotate(
                "positive control",
                xy=(df.loc[idx, "iptm_best"], pos),
                xytext=(df.loc[idx, "iptm_best"] + 0.05, pos),
                fontsize=8, color="red", fontweight="bold",
                arrowprops={"arrowstyle": "->", "color": "red"},
            )

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color_map["high_confidence"], label="High confidence"),
        Patch(facecolor=color_map["predicted"], label="Predicted interaction"),
        Patch(facecolor=color_map["low_confidence"], label="Low confidence"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved: %s", output_path.name)


def generate_summary_table(
    scores_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write a formatted summary CSV with key metrics.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Full interaction scores DataFrame.
    output_path : Path
        Output CSV path.
    """
    summary = scores_df[[
        "target", "uniprot", "iptm_best", "iptm_mean", "ptm_best",
        "mean_interchain_pae", "interaction_call", "n_models_parsed",
        "rank_by_iptm",
    ]].copy()

    summary = summary.sort_values("rank_by_iptm")
    summary.to_csv(output_path, index=False)
    logger.info("Saved: %s", output_path.name)


def generate_stats_report(
    scores_df: pd.DataFrame,
    candidates_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write a text stats report.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Interaction scores.
    candidates_df : pd.DataFrame
        Candidate list.
    output_path : Path
        Output text file path.
    """
    from datetime import datetime

    lines = [
        "Step 03: Structural Prediction -- Summary Statistics",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"== Candidates ==",
        f"  Total candidates: {len(candidates_df)}",
        f"  Positive controls: {len(candidates_df[candidates_df['tier'] == 0])}",
    ]

    # Count by priority
    for tier in sorted(candidates_df["tier"].unique()):
        n = len(candidates_df[candidates_df["tier"] == tier])
        label = "controls" if tier == 0 else f"Tier {tier}"
        lines.append(f"  {label}: {n}")

    lines.extend(["", "== AlphaFold Results =="])
    n_parsed = len(scores_df[scores_df["n_models_parsed"] > 0])
    lines.append(f"  Models parsed: {n_parsed}/{len(scores_df)}")

    for call in ["high_confidence", "predicted", "low_confidence", "not_run", "failed"]:
        n = len(scores_df[scores_df["interaction_call"] == call])
        if n > 0:
            lines.append(f"  {call}: {n}")

    # VEGFA positive control check
    vegfa = scores_df[scores_df["target"] == "VEGFA"]
    lines.extend(["", "== VEGFA Positive Control =="])
    if vegfa.empty:
        lines.append("  VEGFA not in results")
    elif vegfa.iloc[0]["n_models_parsed"] == 0:
        lines.append("  VEGFA: not yet run")
    else:
        iptm = vegfa.iloc[0]["iptm_best"]
        lines.append(f"  VEGFA ipTM: {iptm}")
        if not np.isnan(iptm) and iptm > 0.7:
            lines.append("  PASS: ipTM > 0.7 (method calibrated)")
        elif not np.isnan(iptm):
            lines.append("  WARNING: ipTM <= 0.7 (calibration concern)")

    # Top interactions
    predicted = scores_df[
        scores_df["interaction_call"].isin(["high_confidence", "predicted"])
    ].sort_values("iptm_best", ascending=False)

    if not predicted.empty:
        lines.extend(["", "== Top Predicted Interactions =="])
        for _, row in predicted.head(10).iterrows():
            lines.append(
                f"  {row['target']} ({row['uniprot']}): "
                f"ipTM={row['iptm_best']}, PAE={row['mean_interchain_pae']}, "
                f"call={row['interaction_call']}"
            )

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    logger.info("Saved: %s", output_path.name)
