"""
G1 Gate: FLT1 expression in brain perivascular macrophages (SEA-AD snRNA-seq)

Falsification gate for hyp-034. Tests whether FLT1 mRNA is expressed in
human brain PVM clusters using the SEA-AD Microglia-Immune multi-regional dataset.

Usage:
    python 01_flt1_pvm_gate.py --h5ad <path_to_h5ad> --outdir <output_directory>
"""

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, spearmanr

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Gene panel
GATE_GENES = ["FLT1"]
SCAFFOLD_GENES = ["NRP1", "SDC1", "SDC2", "SDC4", "GPC4"]
SULFATION_GENES = ["HS3ST1", "HS3ST2"]
PVM_MARKERS = ["CD163", "LYVE1", "MRC1", "F13A1"]
MICROGLIA_MARKERS = ["P2RY12", "TMEM119"]
M1_MARKERS = ["IL1B"]
M2_MARKERS = ["CD163", "MRC1"]  # overlap with PVM markers intentional
AMBIENT_CTRL = ["CLDN5", "PECAM1"]  # endothelial contamination controls

ALL_GENES = list(dict.fromkeys(
    GATE_GENES + SCAFFOLD_GENES + SULFATION_GENES +
    PVM_MARKERS + MICROGLIA_MARKERS + M1_MARKERS + M2_MARKERS + AMBIENT_CTRL
))

# Decision thresholds
DETECTION_RATE_INCONCLUSIVE = 0.05  # below this in all myeloid = inconclusive
PVM_ENRICHMENT_PVAL = 0.05  # Wilcoxon threshold for PVM > microglia


def parse_args():
    parser = argparse.ArgumentParser(description="G1 Gate: FLT1 in brain PVMs")
    parser.add_argument(
        "--h5ad", required=True,
        help="Path to SEA-AD Microglia-Immune h5ad file"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory for results and plots"
    )
    return parser.parse_args()


def load_data(h5ad_path: str) -> ad.AnnData:
    """Load h5ad and report basic stats."""
    print(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.n_obs} nuclei x {adata.n_vars} genes")
    print(f"  obs columns: {list(adata.obs.columns[:20])}...")
    return adata


def inspect_taxonomy(adata: ad.AnnData, outdir: Path) -> dict:
    """Report cell type taxonomy and identify PVM-enriched supertypes."""
    info = {}

    # Find the supertype column (may be 'Supertype', 'supertype', etc.)
    supertype_col = None
    subclass_col = None
    for col in adata.obs.columns:
        if "supertype" in col.lower():
            supertype_col = col
        if "subclass" in col.lower():
            subclass_col = col

    if supertype_col is None:
        # Fall back: print all columns for manual inspection
        print("WARNING: No 'supertype' column found. Available columns:")
        print(list(adata.obs.columns))
        info["supertype_col"] = None
        info["subclass_col"] = subclass_col
        return info

    info["supertype_col"] = supertype_col
    info["subclass_col"] = subclass_col

    # Count nuclei per supertype
    st_counts = adata.obs[supertype_col].value_counts()
    print(f"\nSupertype distribution ({supertype_col}):")
    print(st_counts.to_string())
    info["supertype_counts"] = st_counts.to_dict()

    # Save to file
    st_counts.to_csv(outdir / "supertype_counts.csv")

    return info


def identify_pvm_supertypes(adata: ad.AnnData, supertype_col: str) -> tuple:
    """
    Identify PVM-enriched supertypes by marker expression.
    PVM markers: CD163, LYVE1, MRC1, F13A1 (high)
    Microglia markers: P2RY12, TMEM119 (high in microglia, low in PVMs)

    Returns (pvm_supertypes, microglia_supertypes, marker_df)
    """
    supertypes = adata.obs[supertype_col].unique()
    records = []

    for st in supertypes:
        mask = adata.obs[supertype_col] == st
        subset = adata[mask]
        n_cells = mask.sum()

        row = {"supertype": st, "n_cells": n_cells}
        for gene in PVM_MARKERS + MICROGLIA_MARKERS:
            if gene in adata.var_names:
                expr = subset[:, gene].X
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()
                row[f"{gene}_detection"] = (expr > 0).mean()
                row[f"{gene}_mean"] = expr.mean()
            else:
                row[f"{gene}_detection"] = np.nan
                row[f"{gene}_mean"] = np.nan
        records.append(row)

    marker_df = pd.DataFrame(records).set_index("supertype")

    # Score: PVM score = mean detection of PVM markers - mean detection of microglia markers
    pvm_marker_cols = [f"{g}_detection" for g in PVM_MARKERS if f"{g}_detection" in marker_df.columns]
    mic_marker_cols = [f"{g}_detection" for g in MICROGLIA_MARKERS if f"{g}_detection" in marker_df.columns]

    marker_df["pvm_score"] = marker_df[pvm_marker_cols].mean(axis=1)
    marker_df["microglia_score"] = marker_df[mic_marker_cols].mean(axis=1)
    marker_df["pvm_minus_microglia"] = marker_df["pvm_score"] - marker_df["microglia_score"]

    marker_df = marker_df.sort_values("pvm_minus_microglia", ascending=False)

    print("\nPVM vs Microglia marker scores by supertype:")
    print(marker_df[["n_cells", "pvm_score", "microglia_score", "pvm_minus_microglia"]].to_string())

    # Classify: PVM-enriched if pvm_score > microglia_score and pvm_score > 0.1
    pvm_mask = (marker_df["pvm_minus_microglia"] > 0) & (marker_df["pvm_score"] > 0.1)
    mic_mask = (marker_df["pvm_minus_microglia"] < 0) & (marker_df["microglia_score"] > 0.1)

    pvm_supertypes = list(marker_df.index[pvm_mask])
    microglia_supertypes = list(marker_df.index[mic_mask])

    print(f"\nPVM-enriched supertypes: {pvm_supertypes}")
    print(f"Microglia supertypes: {microglia_supertypes}")

    return pvm_supertypes, microglia_supertypes, marker_df


def run_flt1_gate(
    adata: ad.AnnData,
    supertype_col: str,
    pvm_supertypes: list,
    microglia_supertypes: list,
    outdir: Path,
) -> dict:
    """
    Core gate test: FLT1 expression in PVM vs microglia supertypes.
    """
    results = {"gate": "G1_FLT1_PVM", "status": "UNKNOWN"}

    if "FLT1" not in adata.var_names:
        print("CRITICAL: FLT1 not in gene list!")
        results["status"] = "FAIL"
        results["reason"] = "FLT1 not in adata.var_names"
        return results

    # Extract FLT1 expression for PVM and microglia supertypes
    pvm_mask = adata.obs[supertype_col].isin(pvm_supertypes)
    mic_mask = adata.obs[supertype_col].isin(microglia_supertypes)

    pvm_flt1 = adata[pvm_mask, "FLT1"].X
    mic_flt1 = adata[mic_mask, "FLT1"].X
    all_flt1 = adata[:, "FLT1"].X

    if hasattr(pvm_flt1, "toarray"):
        pvm_flt1 = pvm_flt1.toarray().flatten()
    else:
        pvm_flt1 = np.asarray(pvm_flt1).flatten()
    if hasattr(mic_flt1, "toarray"):
        mic_flt1 = mic_flt1.toarray().flatten()
    else:
        mic_flt1 = np.asarray(mic_flt1).flatten()
    if hasattr(all_flt1, "toarray"):
        all_flt1 = all_flt1.toarray().flatten()
    else:
        all_flt1 = np.asarray(all_flt1).flatten()

    # Detection rates
    pvm_detection = (pvm_flt1 > 0).mean()
    mic_detection = (mic_flt1 > 0).mean()
    all_detection = (all_flt1 > 0).mean()

    # Mean expression (among expressing cells)
    pvm_mean = pvm_flt1[pvm_flt1 > 0].mean() if pvm_flt1.any() else 0.0
    mic_mean = mic_flt1[mic_flt1 > 0].mean() if mic_flt1.any() else 0.0

    results["n_pvm_nuclei"] = int(pvm_mask.sum())
    results["n_microglia_nuclei"] = int(mic_mask.sum())
    results["n_total_nuclei"] = adata.n_obs
    results["pvm_supertypes"] = pvm_supertypes
    results["microglia_supertypes"] = microglia_supertypes
    results["flt1_detection_pvm"] = float(pvm_detection)
    results["flt1_detection_microglia"] = float(mic_detection)
    results["flt1_detection_all_myeloid"] = float(all_detection)
    results["flt1_mean_expressing_pvm"] = float(pvm_mean)
    results["flt1_mean_expressing_microglia"] = float(mic_mean)

    print(f"\n{'='*60}")
    print(f"FLT1 GATE RESULTS")
    print(f"{'='*60}")
    print(f"PVM nuclei: {pvm_mask.sum()} ({len(pvm_supertypes)} supertypes)")
    print(f"Microglia nuclei: {mic_mask.sum()} ({len(microglia_supertypes)} supertypes)")
    print(f"FLT1 detection rate - PVM:       {pvm_detection:.4f} ({pvm_detection*100:.1f}%)")
    print(f"FLT1 detection rate - Microglia: {mic_detection:.4f} ({mic_detection*100:.1f}%)")
    print(f"FLT1 detection rate - All:       {all_detection:.4f} ({all_detection*100:.1f}%)")
    print(f"FLT1 mean (expressing) - PVM:       {pvm_mean:.4f}")
    print(f"FLT1 mean (expressing) - Microglia: {mic_mean:.4f}")

    # Statistical test: PVM vs microglia (Wilcoxon rank-sum)
    if pvm_flt1.any() and mic_flt1.any():
        stat, pval = mannwhitneyu(pvm_flt1, mic_flt1, alternative="greater")
        results["wilcoxon_pvm_gt_microglia_stat"] = float(stat)
        results["wilcoxon_pvm_gt_microglia_pval"] = float(pval)
        print(f"Wilcoxon PVM > Microglia: U={stat:.0f}, p={pval:.2e}")
    else:
        results["wilcoxon_pvm_gt_microglia_pval"] = np.nan
        print("Wilcoxon test: insufficient non-zero values")

    # Ambient RNA control: check CLDN5 and PECAM1 in PVM cluster
    ambient_results = {}
    for gene in AMBIENT_CTRL:
        if gene in adata.var_names:
            pvm_expr = adata[pvm_mask, gene].X
            if hasattr(pvm_expr, "toarray"):
                pvm_expr = pvm_expr.toarray().flatten()
            else:
                pvm_expr = np.asarray(pvm_expr).flatten()
            det_rate = (pvm_expr > 0).mean()
            ambient_results[gene] = float(det_rate)
            print(f"Ambient control {gene} in PVM: {det_rate:.4f} ({det_rate*100:.1f}%)")

            # Correlation between FLT1 and endothelial marker in PVM
            if pvm_flt1.any() and pvm_expr.any():
                r, p = spearmanr(pvm_flt1, pvm_expr)
                ambient_results[f"{gene}_flt1_spearman_r"] = float(r)
                ambient_results[f"{gene}_flt1_spearman_p"] = float(p)
                print(f"  FLT1-{gene} Spearman in PVM: r={r:.3f}, p={p:.2e}")

    results["ambient_control"] = ambient_results

    # Decision
    if all_detection < DETECTION_RATE_INCONCLUSIVE:
        results["status"] = "INCONCLUSIVE"
        results["reason"] = (
            f"FLT1 detection rate {all_detection:.3f} below {DETECTION_RATE_INCONCLUSIVE} "
            "in all myeloid cells. snRNA-seq dropout likely. Spatial transcriptomics fallback needed."
        )
    elif pvm_detection > mic_detection and pvm_detection > DETECTION_RATE_INCONCLUSIVE:
        pval = results.get("wilcoxon_pvm_gt_microglia_pval", 1.0)
        if pval < PVM_ENRICHMENT_PVAL:
            results["status"] = "PASS"
            results["reason"] = (
                f"FLT1 detected in PVM supertypes ({pvm_detection:.1%}) > microglia "
                f"({mic_detection:.1%}), Wilcoxon p={pval:.2e}. Gate passes."
            )
        else:
            results["status"] = "PARTIAL"
            results["reason"] = (
                f"FLT1 detected in PVM ({pvm_detection:.1%}) > microglia ({mic_detection:.1%}) "
                f"but Wilcoxon not significant (p={pval:.2e}). Trend supports but not definitive."
            )
    elif pvm_detection > DETECTION_RATE_INCONCLUSIVE:
        results["status"] = "PARTIAL"
        results["reason"] = (
            f"FLT1 detected in PVM ({pvm_detection:.1%}) but not enriched vs microglia "
            f"({mic_detection:.1%}). Local production exists but no perivascular specialization."
        )
    else:
        # Check ambient: if endothelial markers also absent, it's a true negative
        ambient_max = max(ambient_results.get(g, 0) for g in AMBIENT_CTRL)
        if ambient_max < 0.02:
            results["status"] = "FAIL"
            results["reason"] = (
                f"FLT1 undetectable in PVM ({pvm_detection:.1%}), endothelial markers also "
                f"absent (max {ambient_max:.1%}). No ambient RNA excuse. Gate fails."
            )
        else:
            results["status"] = "INCONCLUSIVE"
            results["reason"] = (
                f"FLT1 barely detected in PVM ({pvm_detection:.1%}) but endothelial markers "
                f"present ({ambient_max:.1%}). Ambient RNA may mask true signal."
            )

    print(f"\n>>> GATE STATUS: {results['status']}")
    print(f">>> {results['reason']}")

    return results


def run_scaffold_analysis(
    adata: ad.AnnData,
    supertype_col: str,
    pvm_supertypes: list,
    microglia_supertypes: list,
    outdir: Path,
) -> pd.DataFrame:
    """
    Profile binding scaffold genes (NRP1, SDC1-4, GPC4, HS3ST1-2) across
    PVM vs microglia supertypes. This is C2 piggyback.
    """
    genes = SCAFFOLD_GENES + SULFATION_GENES
    genes_present = [g for g in genes if g in adata.var_names]
    genes_missing = [g for g in genes if g not in adata.var_names]

    if genes_missing:
        print(f"WARNING: Genes not in dataset: {genes_missing}")

    pvm_mask = adata.obs[supertype_col].isin(pvm_supertypes)
    mic_mask = adata.obs[supertype_col].isin(microglia_supertypes)

    records = []
    for gene in genes_present:
        pvm_expr = adata[pvm_mask, gene].X
        mic_expr = adata[mic_mask, gene].X
        if hasattr(pvm_expr, "toarray"):
            pvm_expr = pvm_expr.toarray().flatten()
        else:
            pvm_expr = np.asarray(pvm_expr).flatten()
        if hasattr(mic_expr, "toarray"):
            mic_expr = mic_expr.toarray().flatten()
        else:
            mic_expr = np.asarray(mic_expr).flatten()

        pvm_det = (pvm_expr > 0).mean()
        mic_det = (mic_expr > 0).mean()

        stat_result = mannwhitneyu(pvm_expr, mic_expr, alternative="two-sided") if (pvm_expr.any() or mic_expr.any()) else None

        records.append({
            "gene": gene,
            "pvm_detection": pvm_det,
            "microglia_detection": mic_det,
            "pvm_mean": pvm_expr.mean(),
            "microglia_mean": mic_expr.mean(),
            "pvm_mean_expressing": pvm_expr[pvm_expr > 0].mean() if pvm_expr.any() else 0,
            "microglia_mean_expressing": mic_expr[mic_expr > 0].mean() if mic_expr.any() else 0,
            "wilcoxon_U": float(stat_result.statistic) if stat_result else np.nan,
            "wilcoxon_p": float(stat_result.pvalue) if stat_result else np.nan,
        })

    scaffold_df = pd.DataFrame(records)
    print("\nBinding scaffold gene expression (PVM vs Microglia):")
    print(scaffold_df.to_string(index=False))

    scaffold_df.to_csv(outdir / "scaffold_genes_pvm_vs_microglia.csv", index=False)
    return scaffold_df


def run_m2_correlation(
    adata: ad.AnnData,
    supertype_col: str,
    pvm_supertypes: list,
    outdir: Path,
) -> dict:
    """
    Test FLT1 correlation with M2 markers (MRC1, CD163) and inverse with M1
    markers (IL1B) in PVM clusters. Tests hyp-034 Prediction 4.
    """
    pvm_mask = adata.obs[supertype_col].isin(pvm_supertypes)
    if pvm_mask.sum() == 0:
        return {}

    flt1_expr = adata[pvm_mask, "FLT1"].X
    if hasattr(flt1_expr, "toarray"):
        flt1_expr = flt1_expr.toarray().flatten()
    else:
        flt1_expr = np.asarray(flt1_expr).flatten()

    if not flt1_expr.any():
        print("FLT1 has no expression in PVMs -- skipping correlation analysis")
        return {}

    correlations = {}
    test_genes = M2_MARKERS + M1_MARKERS
    for gene in test_genes:
        if gene not in adata.var_names:
            continue
        gene_expr = adata[pvm_mask, gene].X
        if hasattr(gene_expr, "toarray"):
            gene_expr = gene_expr.toarray().flatten()
        else:
            gene_expr = np.asarray(gene_expr).flatten()

        r, p = spearmanr(flt1_expr, gene_expr)
        category = "M2" if gene in M2_MARKERS else "M1"
        correlations[gene] = {"r": float(r), "p": float(p), "category": category}
        print(f"FLT1 ~ {gene} ({category}) in PVM: r={r:.3f}, p={p:.2e}")

    return correlations


def make_plots(
    adata: ad.AnnData,
    supertype_col: str,
    pvm_supertypes: list,
    microglia_supertypes: list,
    marker_df: pd.DataFrame,
    outdir: Path,
):
    """Generate summary plots."""
    sc.settings.figdir = str(outdir)

    # 1. Dot plot: all genes by supertype
    plot_genes = [g for g in ALL_GENES if g in adata.var_names]
    try:
        sc.pl.dotplot(
            adata, plot_genes, groupby=supertype_col,
            save="_flt1_panel_by_supertype.pdf",
            show=False,
        )
    except Exception as e:
        print(f"Dotplot failed: {e}")

    # 2. Violin plot: FLT1 by supertype
    if "FLT1" in adata.var_names:
        try:
            sc.pl.violin(
                adata, "FLT1", groupby=supertype_col,
                save="_flt1_by_supertype.pdf",
                show=False,
            )
        except Exception as e:
            print(f"Violin plot failed: {e}")

    # 3. Marker heatmap for PVM identification
    fig, ax = plt.subplots(figsize=(10, max(4, len(marker_df) * 0.4)))
    cols = [c for c in marker_df.columns if c.endswith("_detection")]
    plot_data = marker_df[cols].copy()
    plot_data.columns = [c.replace("_detection", "") for c in cols]
    im = ax.imshow(plot_data.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(plot_data.columns)))
    ax.set_xticklabels(plot_data.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(plot_data.index)))
    ax.set_yticklabels(plot_data.index)
    plt.colorbar(im, ax=ax, label="Detection rate")
    ax.set_title("PVM vs Microglia marker detection by supertype")
    plt.tight_layout()
    fig.savefig(outdir / "pvm_marker_heatmap.pdf", dpi=150)
    plt.close(fig)

    print(f"Plots saved to {outdir}")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load data
    adata = load_data(args.h5ad)

    # Inspect taxonomy
    tax_info = inspect_taxonomy(adata, outdir)
    supertype_col = tax_info["supertype_col"]

    if supertype_col is None:
        print("ERROR: Cannot find supertype column. Dumping obs columns:")
        print(list(adata.obs.columns))
        sys.exit(1)

    # Identify PVM vs microglia supertypes by marker expression
    pvm_sts, mic_sts, marker_df = identify_pvm_supertypes(adata, supertype_col)
    marker_df.to_csv(outdir / "supertype_marker_scores.csv")

    if not pvm_sts:
        print("WARNING: No PVM-enriched supertypes found. Running gate on all myeloid.")
        pvm_sts = list(adata.obs[supertype_col].unique())
        mic_sts = []

    # Run FLT1 gate
    gate_results = run_flt1_gate(adata, supertype_col, pvm_sts, mic_sts, outdir)

    # Run scaffold analysis (C2 piggyback)
    scaffold_df = run_scaffold_analysis(adata, supertype_col, pvm_sts, mic_sts, outdir)

    # Run M2 correlation (hyp-034 P4)
    m2_corr = run_m2_correlation(adata, supertype_col, pvm_sts, outdir)
    gate_results["m2_correlations"] = m2_corr

    # Generate plots
    make_plots(adata, supertype_col, pvm_sts, mic_sts, marker_df, outdir)

    # Save results
    results_path = outdir / "g1_gate_results.json"
    with open(results_path, "w") as f:
        json.dump(gate_results, f, indent=2, default=str)
    print(f"\nResults saved to {results_path}")

    # Print final verdict
    print(f"\n{'='*60}")
    print(f"FINAL VERDICT: {gate_results['status']}")
    print(f"{'='*60}")
    print(gate_results["reason"])

    # Exit code: 0 for PASS/PARTIAL, 1 for FAIL/INCONCLUSIVE
    if gate_results["status"] in ("PASS", "PARTIAL"):
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
