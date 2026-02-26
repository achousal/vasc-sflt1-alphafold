"""01_select_candidates.py -- Rank candidate proteins for AlphaFold Multimer prediction.

Priority order:
  1. Tier 1, dual-somamer, AND in enriched axon/semaphorin/neuropilin pathway
  2. Tier 1, dual-somamer
  3. Tier 2, in enriched axon/semaphorin pathway
  4. Positive controls: VEGFA, NRP1, NRP2, SEMA3A
"""

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Regex for axon guidance pathway terms in enrichment results
AXON_PATHWAY_REGEX = (
    r"axon|semaphorin|plexin|neuropilin|ephrin|guidance|slit|robo|netrin"
)

# Mandatory positive controls with UniProt accessions
POSITIVE_CONTROLS = [
    {"target": "VEGFA", "uniprot": "P15692", "entrez_gene_symbol": "VEGFA",
     "rationale": "Known sFLT1 ligand -- pipeline calibration positive control"},
    {"target": "NRP1", "uniprot": "O14786", "entrez_gene_symbol": "NRP1",
     "rationale": "Neuropilin-1, known VEGF/semaphorin co-receptor"},
    {"target": "NRP2", "uniprot": "O60462", "entrez_gene_symbol": "NRP2",
     "rationale": "Neuropilin-2, known VEGF/semaphorin co-receptor"},
    {"target": "SEMA3A", "uniprot": "Q14563", "entrez_gene_symbol": "SEMA3A",
     "rationale": "Class 3 semaphorin, canonical axon guidance ligand"},
]


def extract_pathway_genes(enrichment_dir: Path) -> set[str]:
    """Extract gene symbols from axon/semaphorin enrichment terms.

    Scans all enrichment CSVs in the directory for terms matching the
    axon pathway regex and collects gene symbols from the geneID column.

    Parameters
    ----------
    enrichment_dir : Path
        Directory containing step02_enrichment_*.csv files.

    Returns
    -------
    set[str]
        Gene symbols found in matching pathway terms.
    """
    pathway_genes = set()
    enrichment_files = list(enrichment_dir.glob("step02_enrichment_*_pos_combined.csv"))

    for fpath in enrichment_files:
        try:
            df = pd.read_csv(fpath)
        except Exception as e:
            logger.warning("Failed to read %s: %s", fpath.name, e)
            continue

        if "Description" not in df.columns or "geneID" not in df.columns:
            continue

        mask = df["Description"].str.contains(
            AXON_PATHWAY_REGEX, case=False, na=False
        )
        for genes_str in df.loc[mask, "geneID"].dropna():
            pathway_genes.update(genes_str.split("/"))

    logger.info("Extracted %d genes from axon/guidance pathway terms", len(pathway_genes))
    return pathway_genes


def select_candidates(
    consensus_path: Path,
    enrichment_dir: Path,
    n_candidates: int = 20,
) -> pd.DataFrame:
    """Select and rank candidate proteins for structural prediction.

    Parameters
    ----------
    consensus_path : Path
        Path to step01_consensus_proteins_pos.csv.
    enrichment_dir : Path
        Directory containing Step 2 enrichment CSVs.
    n_candidates : int
        Maximum number of candidates to select (excluding positive controls).

    Returns
    -------
    pd.DataFrame
        Ranked candidates with columns: rank, target, uniprot,
        entrez_gene_symbol, tier, dual_somamer, in_axon_pathway, rationale.
    """
    consensus = pd.read_csv(consensus_path)
    pathway_genes = extract_pathway_genes(enrichment_dir)

    # Deduplicate: keep the best-tier row per target
    consensus = consensus.sort_values(["tier", "dual_somamer"], ascending=[True, False])
    consensus = consensus.drop_duplicates(subset="target", keep="first")

    # Flag pathway membership
    consensus["in_axon_pathway"] = consensus["entrez_gene_symbol"].isin(pathway_genes)

    # Assign priority
    def assign_priority(row):
        if row["tier"] == 1 and row["dual_somamer"] and row["in_axon_pathway"]:
            return 1
        if row["tier"] == 1 and row["dual_somamer"]:
            return 2
        if row["tier"] <= 2 and row["in_axon_pathway"]:
            return 3
        return 4

    consensus["priority"] = consensus.apply(assign_priority, axis=1)

    # Sort by priority then by mean_beta descending
    consensus = consensus.sort_values(
        ["priority", "mean_beta"], ascending=[True, False]
    )

    # Select top candidates
    candidates = consensus.head(n_candidates).copy()

    # Build rationale
    def build_rationale(row):
        parts = []
        parts.append(f"Tier {row['tier']}")
        if row["dual_somamer"]:
            parts.append("dual-somamer")
        if row["in_axon_pathway"]:
            parts.append("in axon/semaphorin pathway")
        parts.append(f"priority {row['priority']}")
        return "; ".join(parts)

    candidates["rationale"] = candidates.apply(build_rationale, axis=1)

    # Add positive controls (if not already in candidates)
    existing_uniprots = set(candidates["uniprot"].dropna())
    controls_to_add = []
    for ctrl in POSITIVE_CONTROLS:
        if ctrl["uniprot"] not in existing_uniprots:
            controls_to_add.append({
                "target": ctrl["target"],
                "uniprot": ctrl["uniprot"],
                "entrez_gene_symbol": ctrl["entrez_gene_symbol"],
                "tier": 0,
                "dual_somamer": False,
                "in_axon_pathway": True,
                "rationale": ctrl["rationale"],
                "priority": 0,
            })

    if controls_to_add:
        ctrl_df = pd.DataFrame(controls_to_add)
        candidates = pd.concat([ctrl_df, candidates], ignore_index=True)

    # Assign rank
    candidates = candidates.reset_index(drop=True)
    candidates["rank"] = range(1, len(candidates) + 1)

    # Select output columns
    out_cols = [
        "rank", "target", "uniprot", "entrez_gene_symbol",
        "tier", "dual_somamer", "in_axon_pathway", "rationale",
    ]
    result = candidates[out_cols].copy()

    logger.info("Selected %d candidates (including positive controls)", len(result))
    return result
