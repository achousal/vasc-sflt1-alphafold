# Step 4: FLT1 Expression in Brain Perivascular Macrophages (G1 Gate)

**Purpose**: Falsification gate for hyp-034. Test whether FLT1 mRNA is expressed
in human brain perivascular macrophage clusters using SEA-AD snRNA-seq data.

## Data Source

SEA-AD Microglia-and-Immune multi-regional release (AAIC pre-release, 2025-07-24)
- S3: `s3://sea-ad-single-cell-profiling/Microglia-and-Immune-for-AAIC/`
- File: `SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC-pre-release.2025-07-24.h5ad`
- Size: ~3.0 GB
- Nuclei: 240,651 (majority microglia, plus PVMs, lymphocytes)
- Donors: 84 aged, spanning AD continuum
- Regions: 10 brain regions (MTG, DLPFC, hippocampus, etc.)

## Gene Panel

| Category | Genes | Purpose |
|----------|-------|---------|
| Gate target | FLT1 | Primary gate -- is it expressed in PVMs? |
| Binding scaffold | NRP1, SDC1, SDC2, SDC4, GPC4 | sFLT1 binding machinery |
| Sulfation code | HS3ST1, HS3ST2 | 3-O sulfation for NRP1 binding specificity |
| PVM markers | CD163, LYVE1, MRC1, F13A1 | Identify PVM-enriched supertypes |
| Microglia markers | P2RY12, TMEM119 | Parenchymal microglia reference |
| M2 markers | CD163, MRC1 | Polarization correlation |
| M1 markers | IL1B | Polarization correlation |
| Ambient RNA control | CLDN5, PECAM1 | Endothelial contamination check |

## Decision Rules

| Outcome | Criteria | Consequence |
|---------|----------|-------------|
| **PASS** | FLT1 detected in PVM-enriched supertype(s) above background, PVM > microglia | hyp-034 gate passes |
| **PARTIAL** | FLT1 detected but PVM ≈ microglia | Paracrine niche specificity weakened |
| **INCONCLUSIVE** | FLT1 < 5% detection rate in all myeloid | snRNA-seq dropout; need spatial fallback |
| **FAIL** | FLT1 undetectable AND endothelial markers absent (no ambient RNA excuse) | hyp-034 self-terminates |

## Environment

- HPC: Minerva (Mount Sinai)
- Conda env: `/sc/arion/work/chousa01/envs/cellxgene_census`
- Packages: scanpy 1.11.5, anndata 0.12.10, scipy 1.17.1, matplotlib 3.10.8

## Scripts

| Script | Purpose |
|--------|---------|
| `00_download_data.sh` | Download h5ad from S3 to HPC |
| `01_flt1_pvm_gate.py` | Main analysis: FLT1 expression, PVM identification, statistics |
| `01_flt1_pvm_gate.lsf` | LSF job submission script |
