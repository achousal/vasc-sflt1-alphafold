# SFLT1_VEGF

## Overview

sFLT1 (soluble VEGF receptor 1) is a decoy receptor that captures VEGFA and blocks angiogenesis. This project tests whether a similar decoy-receptor mechanism applies to axonogenesis: sFLT1 may sequester ligands needed for axonal guidance proteins (semaphorins, neuropilins, NRP1/NRP2), preventing downstream signaling.

The project identifies proteins significantly associated with sFLT1 across multiple cohorts using SomaScan proteomics, then screens candidate interactors for structural plausibility using AlphaFold and crystal structure prediction.

## Scientific Question

Which proteins are mechanistically connected to sFLT1 in a way that inhibits axonogenesis-related downstream signaling, and can structural prediction methods validate these interactions?

## Data

All data is SomaScan proteomics. Linear models regressed each protein on sFLT1 levels (two somamers: VEGFsR1 and VEGFsR1.1). Results are pre-computed and stored as significant positive/negative association tables.

### Cohorts

| Cohort | Role | Age-adjusted | Dir |
|--------|------|--------------|-----|
| MarkVCID | Discovery | No | data/MarkVCID/ |
| UCSF_AD | Discovery | Yes | data/UCSF_AD/ |
| GNPC | Validation | No | data/GNPC/ |
| WASHU | Validation | Yes | data/WASHU/ |

### File naming convention

- `values_LM_VEGFsR1_sig_pos.csv` -- significant positive associations, VEGFsR1 somamer
- `values_LM_VEGFsR1_sig_neg.csv` -- significant negative associations, VEGFsR1 somamer
- `values_LM_VEGFsR1.1_sig_pos.csv` -- significant positive associations, VEGFsR1.1 somamer
- `values_LM_VEGFsR1.1_sig_neg.csv` -- significant negative associations, VEGFsR1.1 somamer
- `discovery_*_proteins_common.csv` / `validation_*_proteins_common.csv` -- protein annotation tables

### Column schema (LM results)

| Column | Description |
|--------|-------------|
| Target | SomaScan target short name |
| beta | Effect size from linear model |
| pvalue | Raw p-value |
| padj | Adjusted p-value (BH) |
| minus_log_padj | -log10(padj) |
| sig_prot | Significance flag |
| Significant_Direction | Positive or Negative |

### Column schema (protein annotations)

| Column | Description |
|--------|-------------|
| name | SomaScan analyte name |
| Analytes / column_name | SomaScan column identifier |
| TargetFullName / target_full_name | Full protein name |
| Target / target | Short protein name |
| UniProt | UniProt accession |
| EntrezGeneID / entrez_gene_id | Entrez gene ID |
| EntrezGeneSymbol / entrez_gene_symbol | Gene symbol |

## Analysis Pipeline

### Step 1: Cross-cohort overlap (R)
- `analysis/01_cross_cohort_overlap/`
- Identify proteins replicated across discovery (MarkVCID, UCSF_AD) and validation (GNPC, WASHU) cohorts
- Stratify by somamer (VEGFsR1 vs VEGFsR1.1) and direction (pos/neg)
- Generate UpSet plots, Venn diagrams, and ranked consensus lists
- Output: consensus protein lists for downstream analysis

### Step 2: Pathway enrichment (R)
- `analysis/02_pathway_enrichment/`
- GO, KEGG, Reactome enrichment on consensus protein sets
- Separate enrichment for positive vs negative associations
- Focus on axon guidance, semaphorin signaling, neuropilin pathways
- Output: enrichment tables and dot plots

### Step 3: Structural prediction (Python)
- `analysis/03_structural_prediction/`
- AlphaFold multimer prediction for sFLT1 + top candidate proteins
- Evaluate predicted alignment error (PAE) and interface confidence (ipTM)
- Compare with existing crystal structures where available
- Output: interaction scores, PAE plots, structural models

## Languages

- **R**: Steps 1-2 (data wrangling, overlap analysis, pathway enrichment, plotting)
- **Python**: Step 3 (AlphaFold API, structural bioinformatics)

## HPC

- Cluster: Minerva (Mount Sinai)
- Scheduler: LSF
- AlphaFold jobs require GPU nodes

## How to run

```bash
# Step 1: Cross-cohort overlap
Rscript analysis/01_cross_cohort_overlap/run_overlap.R

# Step 2: Pathway enrichment
Rscript analysis/02_pathway_enrichment/run_enrichment.R

# Step 3: Structural prediction (HPC)
bsub < analysis/03_structural_prediction/submit_alphafold.lsf
```

## How to test

Tests live alongside analysis code in each step directory.

```bash
# R tests (per step)
Rscript -e "testthat::test_file('analysis/01_cross_cohort_overlap/test_overlap.R')"
Rscript -e "testthat::test_file('analysis/02_pathway_enrichment/test_enrichment.R')"

# Python tests (per step)
pytest analysis/03_structural_prediction/ -v
```
