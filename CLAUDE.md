# SFLT1_VEGF

sFLT1/VEGF interaction screen: SomaScan proteomics cross-cohort overlap -> pathway enrichment -> AlphaFold multimer structural prediction.

## Data

All data is SomaScan proteomics. Linear models regressed each protein on sFLT1 levels (two somamers: VEGFsR1 and VEGFsR1.1). Results are pre-computed.

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

| Step | Dir | Language | Output |
|------|-----|----------|--------|
| 1. Cross-cohort overlap | `analysis/01_cross_cohort_overlap/` | R | Consensus protein lists |
| 2. Pathway enrichment | `analysis/02_pathway_enrichment/` | R | Enrichment tables, dot plots |
| 3. Structural prediction | `analysis/03_structural_prediction/` | Python | Interaction scores, PAE plots, models |

## HPC

- Cluster: Minerva (Mount Sinai), Scheduler: LSF
- AlphaFold jobs require GPU nodes
- Memory: mem=32000 (128 GB) for large partners (>500 aa); mem=16000 may suffice for <500 aa

## AlphaFold Guardrails

### Protein localization filter (mandatory before modeling)

sFLT1 is soluble/extracellular. AlphaFold predicts interfaces regardless of co-localization. Filter candidates:

| Protein Type | AlphaFold Input | Action |
|-------------|-----------------|--------|
| Soluble / secreted | Full-length sequence | Model as-is |
| Transmembrane | Extracellular domain only | Truncate; extract boundaries from UniProt topology; exclude signal peptides; flag if ectodomain <50 aa |
| GPI-anchored | Full ectodomain | Model ectodomain; note surface constraint |
| Nuclear / cytoplasmic | Flag, do not prioritize | High ipTM does NOT mean biological interaction |

Annotation sources: UniProt subcellular location, GO cellular component (GO:0005576, GO:0016021, GO:0005634), Human Protein Atlas.

### PAE thresholding

No single global PAE threshold works. Calibrate using known sFLT1 partners (VEGFA, PlGF) as positive controls and known non-interactors as negatives. PAE varies by protein size, domain flexibility, and structural class.

### Template bias

High PAE confidence may reflect PDB template abundance for sFLT1 Ig-like domains, not genuine binding. Correlate PAE scores with template sequence identity. Asymmetric template coverage (well-templated sFLT1 vs poorly-templated partner) can generate artificially low PAE.

### Binary screen blind spots

VEGF-bridged co-receptor interactions (NRP1, NRP2) are invisible to binary predictions. Strong hits (ipTM >0.5) are likely direct binary interactors only.

## Upstream Constraints

- **Batch effects**: ComBat-seq harmonization required before cross-cohort testing. Can alter rankings by artifact.
- **Comorbidity confounding**: sFLT1 associations may reflect shared vascular risk, not direct mechanism.
- **Multiple testing**: BH-FDR at alpha 0.05 for all candidate screens.

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

```bash
# R tests
Rscript -e "testthat::test_file('analysis/01_cross_cohort_overlap/test_overlap.R')"
Rscript -e "testthat::test_file('analysis/02_pathway_enrichment/test_enrichment.R')"

# Python tests
pytest analysis/03_structural_prediction/ -v
```
