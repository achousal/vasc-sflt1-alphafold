# Upstream Methods: MarkVCID CSF Proteomics

Source: `markVCID_angiogenesis_synapticPlasticity.R` + `MarkVCID_CSF_Protein_metadata.csv`
Documented: 2026-03-19

## Purpose

Documents the upstream MarkVCID analysis that produced the LM results consumed by our Step 1 cross-cohort overlap pipeline. This is the analysis that generated the MarkVCID significance list (proteins significant for sFLT1 association).

## Data Inputs

| File | Description |
|------|-------------|
| `MarkVCID_CSF_Proteomics_data.csv` | SomaScan aptamer-based CSF proteomics (samples x SOMAmers) |
| `MarkVCID_CSF_Protein_metadata.csv` | SOMAmer annotation: 19 fields x 7,596 analytes (SeqId, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, Organism, + QC/calibration) |
| `MarkVCID with CSF samples - Somalogic.xlsx` | Clinical metadata per participant |
| `MarkVCID1_anonID_subject_map.xlsx` | Anonymized label -> ExtIdentifier mapping |

## Sample Pipeline

1. Filter `SampleType == "Sample"`, transpose, merge with protein metadata
2. Restrict to `Organism == "Human"` -> **7,335 proteins**
3. Clinical: compute `mean_age` from 6 repeated age columns (max difference 1.2 years across visits)
4. Build `DxPrim` from NACC-style binary diagnosis columns (ALZDIS, LBDIS, PARK, CVD, HYCEPH, MEDS, etc.)
   - Priority rules: Primary > Contributing > Non-contributing > Unknown
   - Excluded: brain injury (BRNINJ), drug-induced cognitive impairment (c_DYSILL), essential tremor (uk_ESSTREM)
   - No diagnosis -> `NORMAL`
5. Binary case status: `ND` (any neurological diagnosis) vs `NORMAL`
6. **N post-exclusion: ~69 participants** (69 ND + NORMAL count combined)

## Preprocessing

- **log10 transform** all protein columns
- **z-scale** per protein (column-wise)
- **Outlier check:** Grubbs test on sFLT1 soma 1 identified one extreme value (EXID40000010371278). Mean protein expression outlier is different participant. Decision: **retained** -- not a global expression artifact.

## Key Variables

| Variable | Encoding | Notes |
|----------|----------|-------|
| Age (mean_age) | Continuous | Averaged across visit timepoints |
| Gender | 1=Male, 2=Female | |
| DxPrim | Categorical | Granular diagnosis |
| ND_NORMAL | Binary: ND vs NORMAL | Primary grouping |
| Case_Sex | 4-level: ND_1, ND_2, NORMAL_1, NORMAL_2 | Sex-stratified case status |
| HYPERTEN | Binary | **Available but NOT used as covariate** |
| FZKOVRAL | Ordinal (0-3) | Fazekas overall WMH score. **Used only for sFLT1 group analysis, NOT in proteome-wide LMs** |
| FZKDWM, FZKDWMLC, FZKPVWM | Ordinal | Additional Fazekas subscores. Not used in LMs |

## sFLT1 Target

Two SOMAmers measure sFLT1: `VEGF sR1` (soma 1) and `VEGF sR1.1` (soma 2). Both show Pearson correlation. All analyses run in parallel for reproducibility; final results use intersection of significant hits.

## Proteome-Wide Linear Models (the upstream LMs we consume)

```r
# For each of 7,335 proteins:
lm(protein_i ~ VEGF_sR1 + Age)   # soma 1
lm(protein_i ~ VEGF_sR1.1 + Age) # soma 2
```

**Covariates: Age only.** No sex adjustment, no vascular comorbidities, no Fazekas.

**Significance filter:**
- Bonferroni-adjusted p <= 0.05
- |beta| >= 0.5 (hard effect-size filter)

Separate positive-beta and negative-beta significant protein lists per somamer. Final lists: intersection across both somamers (`pos_proteins_common`, `neg_proteins_common`).

## Pathway Enrichment (upstream, not ours)

### ORA (enrichGO/enrichKEGG)
- Applied to common positive-beta and negative-beta protein lists
- GO BP top hit (positive): **Axonogenesis** (98/561 genes)
- Followed by: axon guidance, synapse organization
- Semantic similarity simplification at cutoff=0.7

### GSEA on PCA loadings
- Input: all 7,335 proteins -> PCA -> PC1/PC2 loadings as ranked gene list
- `gseGO()` on PC1/PC2 loading vectors (BP ontology, BH correction)
- `gseKEGG()` on PC1 loadings
- Key terms: GO:0007411 (axon guidance), GO:0007409 (axonogenesis), GO:0007416 (synapse assembly)
- KEGG: hsa04360 (Axon guidance), hsa04130 (SNARE), hsa04141 (ER processing)

## Focused Protein Analyses

| Protein | LM with sFLT1 | Sex interaction | Case interaction |
|---------|---------------|-----------------|------------------|
| NRP1 | p < 2.2e-16 (both somamers) | NS (p=0.43) | Not reported |
| PLXNA1 | Significant | **p=0.02** (soma 1), **p=0.047** (soma 2) | Not reported |
| PLXNA4 | Significant | **p=0.02** (both somamers) | Not reported |
| PTPRJ | Significant | Trend p=0.08 (soma 1), NS (soma 2) | Not reported |
| NRP2 | Tested | Not reported | Not reported |
| Stathmin (STMN1) | Tested | Not tested | Not tested |

## sFLT1 Group Comparisons

- **ANOVA/ANCOVA:** sFLT1 ~ Case_Sex + Age (4-level factor)
- **Fazekas stratification:** sFLT1 ~ ND_NORMAL_FZK + Age (ND/NORMAL crossed with Fazekas overall)
  - Borderline: NORMAL_0 vs ND_3 (Tukey p=0.046)

## Implications for Our Pipeline

### Resolved Questions
1. **T0.2:** MarkVCID covariates = **Age only**. Hypertension available but unused. Vascular comorbidities not adjusted.
2. **T0.3:** MarkVCID N = **~69 post-exclusion** (smaller than estimated 100-300). This further confirms small-N discovery cohort gatekeeping.
3. **Processing:** log10 + z-scale. No batch correction, no ComBat-seq.

### New Concerns
1. **Vascular confounding confirmed absent:** HYPERTEN and Fazekas are in the data but NOT covariates in the proteome-wide LMs. The vault claim "vascular comorbidity burden confounds sFLT1-axonogenesis pathway associations" is directly relevant.
2. **Sex not adjusted in proteome-wide LMs:** but sex interactions exist for key targets (PLXNA1, PLXNA4). Our AF2 screen treats all targets equally regardless of sex interaction.
3. **N=69 is very small** for a 7,335-protein screen with Bonferroni correction. The |beta| > 0.5 filter does heavy lifting here.
4. **Two-somamer intersection** is a robustness check unique to SomaScan. Strengthens confidence in the target list but halves discoverable associations.

### Actions
- [x] Document MarkVCID methods and covariates
- [ ] Update T0.2 in post-af2-analysis.md with resolved status
- [ ] Flag vascular confounding as limitation in all downstream reports
- [ ] Consider sex-stratified AF2 interpretation for PLXNA1/PLXNA4 (sex-interacting targets)
- [ ] Ask PI: GNPC, WASHU, UCSF_AD upstream methods (do they match this pattern?)
