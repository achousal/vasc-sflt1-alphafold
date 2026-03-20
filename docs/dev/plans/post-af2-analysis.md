# Plan: Post-AF2 Analysis and Biological Validation

Updated: 2026-03-19

## Goal

With 18/24 AF2 Multimer jobs complete and 6 finishing by ~Mar 23, build the downstream analysis layers: template bias check, domain-resolved interface mapping, alternative scoring, biological cross-validation, and upstream refinements.

## Phase 0: Upstream Validation (CRITICAL -- before interpreting AF2 results)

Literature grounding review (2026-03-19) identified that the 24 AF2 targets derive from cross-cohort overlap of pre-computed linear models. Four upstream concerns must be resolved before AF2 results can be interpreted with confidence. See ADR-003.

### T0.1 -- Batch Correction Audit (~1h)
**Start immediately. Blocks interpretation of all downstream results.**

**Audit findings (2026-03-19):** The data files in `data/{cohort}/` are **pre-computed LM results** received from upstream collaborators. The Step 1 R code (`01_load_helpers.R`, `02_compute_overlap.R`) performs NO batch correction -- it loads CSVs and computes set overlap. Key observations:

1. **No batch correction anywhere in the local pipeline.** The code reads `values_LM_*.csv` files and intersects significant target lists. No ComBat-seq, no limma batch removal, no harmonization step exists.
2. **Heterogeneous covariate adjustment:** UCSF_AD files carry `_age` suffix (age-adjusted LMs). MarkVCID and GNPC have no suffix (unadjusted). WASHU LM files have no suffix but annotation files say `_adj_for_age`.
3. **Massively different p-value magnitudes:** GNPC/WASHU show p-values in the 1e-170 to 1e-249 range; MarkVCID shows 1e-11 to 1e-14. This implies N differs by orders of magnitude across cohorts.
4. **Asymmetric negative direction:** GNPC has only 3 negative-direction significant proteins vs MarkVCID's 318, making negative-direction overlap nearly empty.

**Implication:** The Tier 1+2 consensus list was derived from significance overlap of independently-run, heterogeneously-adjusted, unharmonized LMs across cohorts with vastly different sample sizes. The overlap approach (significant in all 4) is robust to p-value inflation from large N but vulnerable to:
- Batch artifacts that replicate across all 4 cohorts (e.g., SomaScan platform-level aptamer behavior)
- Confounders included in some LMs but not others (age in UCSF_AD but not MarkVCID)
- Vascular comorbidities not adjusted in any cohort

**Decision needed:**
- [ ] Confirm with PI: what covariates were in the upstream LMs for MarkVCID, GNPC, WASHU?
- [ ] Confirm N per cohort (see T0.3)
- [ ] Assess whether cross-cohort replication provides sufficient protection against batch artifacts, OR whether ComBat-seq re-analysis (T5.1) must precede interpretation
- [ ] Document whether the 4-cohort consensus filtering approach (rather than formal meta-analysis) was a deliberate design choice and its limitations

### T0.2 -- Vascular Comorbidity Covariate Check (~30min)
**RESOLVED for MarkVCID (2026-03-19). Partially resolved for others.**

Source: `markVCID_angiogenesis_synapticPlasticity.R` review.

| Cohort | Covariates in LMs | Vascular comorbidities? | Source |
|--------|-------------------|------------------------|--------|
| MarkVCID | Age only | **NO.** HYPERTEN + Fazekas WMH scores available in clinical data but NOT used as covariates in proteome-wide LMs | R script audit |
| UCSF_AD | Age | Unknown | File suffix `_age` |
| GNPC | Unknown | Unknown | Pre-computed LM CSVs |
| WASHU | Likely age | Unknown | Annotation `_adj_for_age` |

**Key finding:** MarkVCID has hypertension and Fazekas WMH data available (`HYPERTEN`, `FZKOVRAL`, `FZKDWM`, `FZKDWMLC`, `FZKPVWM`) but none were included in the proteome-wide `lm(protein_i ~ VEGF_sR1 + Age)` models. This directly validates the vault concern: "vascular comorbidity burden confounds sFLT1-axonogenesis pathway associations in cross-sectional SomaScan cohorts."

**Implication:** If vascular comorbidities were not adjusted in ANY of the 4 cohorts, the 4-cohort consensus list may be enriched for vascular-confounded associations that replicate precisely because the confound is universal. The |beta| > 0.5 effect-size filter does not protect against this -- a confound-driven association can have a large effect size.

- [x] MarkVCID covariate inclusion documented (Age only, vascular absent)
- [ ] Confirm GNPC, WASHU, UCSF_AD covariates with PI
- [x] Flag as limitation: vascular confounding is a systemic risk across all cohorts
- See: `docs/dev/plans/upstream-markvcid-methods.md` for full R script audit

### T0.3 -- Sample Size Documentation (~30min)
**RESOLVED for MarkVCID (2026-03-19). Updated with R script audit.**

N per cohort, with MarkVCID now confirmed from source code:

| Cohort | Role | N (confirmed) | Published N | Age-adj? | Covariates | SomaScan version |
|--------|------|---------------|-------------|----------|------------|------------------|
| MarkVCID | Discovery | **~69 post-exclusion** | 653 (Phase 1 total, multi-site) | No | **Age only** | 7,596 analytes (7,335 human) |
| UCSF_AD | Discovery | ~31-60 (estimated) | ~200-500 (ADRC subset) | Yes (age) | Age only | Unknown |
| GNPC | Validation | ~2,400 (estimated) | Up to 18,645 (V1 harmonized, 23 cohorts) | No | Unknown | 7K v4.1 (26,458 assays) |
| WASHU | Validation | ~2,200 (estimated) | ~1,100-3,300 (Knight ADRC) | Likely age | Unknown | 7K |

**MarkVCID detail from R script audit:** 7,596 total SOMAmers, 7,335 after filtering to human proteins. Clinical exclusions: brain injury, drug-induced cognitive impairment, essential tremor. ~69 participants remain (ND + NORMAL combined). Proteome-wide LMs use Bonferroni correction across 7,335 tests with |beta| > 0.5 effect-size gate.

**Key findings:**
- **GNPC is a mega-consortium** (Global Neurodegeneration Proteomics Consortium): 18,645 participants, 23 cohorts, 31,083 samples. The LM results may come from a subset, but the p-value magnitudes (1e-195 to 1e-249) confirm N >> 1,000.
- **WASHU (Knight ADRC)**: ~1,100-3,300 participants with SomaScan, consistent with p-value magnitudes (1e-163).
- **MarkVCID Phase 1**: 653 total across 6 sites. SomaScan subset could be smaller (~100-300). P-values (1e-9 to 1e-14) consistent with N ~100-300.
- **UCSF_AD**: Smallest. P-values (1e-6 to 1e-7) suggest N ~50-200.

**Power imbalance is ~10-50x between validation and discovery cohorts.** This means:
1. GNPC/WASHU detect much weaker effects as significant, inflating their significant protein lists
2. The |beta| > 0.5 filter partially mitigates this (enforces minimum effect size)
3. But the Tier 1 requirement (significant in ALL 4) means MarkVCID/UCSF_AD are the bottleneck -- proteins must survive the smallest-N filter
4. This is actually a **conservative** design: the small discovery cohorts gate the target list, not the large validation cohorts

**Assessment:** The 4-cohort overlap approach is more robust than initially feared. The small-N discovery cohorts act as a stringency filter. However, the |beta| > 0.5 effect-size filter was applied upstream and is the primary safeguard -- without it, the large GNPC/WASHU cohorts would contribute thousands of trivially-significant proteins. The main remaining risk is confounding (batch artifacts or comorbidity effects that replicate across all 4 cohorts).

- [x] N per cohort estimated (see table above)
- [ ] Confirm exact N with PI (especially MarkVCID SomaScan subset and UCSF_AD subset)
- [x] Power imbalance assessed: conservative design, small cohorts gate

### T0.4 -- VEGF-Depletion Inversion Test (~2h)
**Independent. Zero new data needed. High impact.**

- Script: `analysis/02_pathway_enrichment/run_vegf_adjustment.R` (new)
- Test vault inversion: "sflt1-associated reductions in axonogenesis pathway scores reflect VEGF pathway suppression rather than direct sflt1 activity"
- Counter-evidence exists: "plasma sFLT1 predicts axonal injury markers independently of VEGF across multi-cohort SomaScan proteomics"
- Method: mutual VEGF/sFLT1 adjustment on pathway enrichment scores
- If VEGF adjustment eliminates signal: the VEGF-independent mechanism (hyp-004/005/016/028/031) collapses
- If signal persists: strongest evidence yet for direct sFLT1 activity
- Acceptance:
  - [ ] Pathway enrichment with and without VEGF adjustment
  - [ ] Effect attenuation quantified (% reduction in enrichment scores)
  - [ ] Inversion claim updated in vault with result

## Phase 1: Immediate post-AF2 analysis (no new GPU jobs)

### T1.2 -- Template Bias Quantification (~2h)
**Start now on 18 targets. Rerun trivially when 24 complete.**
**Why first:** if r > 0.3, template bias confounds all downstream interpretation.

- Script: `analysis/03_structural_prediction/07_template_bias_analysis.py`
- Output: `results/03_structural_prediction/step03_template_bias_report.txt`, scatter plot
- Parse `msas/` dirs from each AF2 job for PDB template hits
- Spearman correlation: ipTM rank vs template count
- Acceptance:
  - [ ] Spearman rho + p-value reported
  - [ ] If r > 0.3: bias flagged in report + downstream outputs carry caveat
  - [ ] Scatter: ipTM vs template count with regression line + CI, VEGFA labeled
- Hypothesis link: hyp-015 (falsification-gated dual-filter)

### T1.1 -- Domain-Resolved Interface Mapping (~4h)
**Start now on 18 targets. Final version needs all 24.**

- Script: `analysis/03_structural_prediction/06_domain_interface_analysis.py`
- Output: `results/03_structural_prediction/step03_domain_interface_scores.csv`, domain PAE heatmap
- sFLT1 domain boundaries: D1 (1-110), D2 (111-220), D3 (221-338)
- Per target: inter-chain PAE matrix, per-domain mean PAE + interface residue count (PAE < 12A)
- Classify: D2-engaging (VEGF-competitive) vs D1/D3-engaging (non-competitive)
- Cross-tabulate against partner subcellular localization (UniProt annotation)
- Acceptance:
  - [ ] VEGFA maps predominantly to D2 (positive control)
  - [ ] CSV: target, d1_mean_pae, d2_mean_pae, d3_mean_pae, d1/d2/d3_interface_res, domain_class, localization
  - [ ] Heatmap: rows=targets (sorted by ipTM), cols=domains, fill=mean inter-chain PAE
- Hypothesis link: hyp-031 (domain-resolved binding mode classification)

### T1.3 -- Alternative Scoring Metrics (~3h)
**Depends on T1.1 (uses interface residue definitions).**

- Script: `analysis/03_structural_prediction/08_alternative_scores.py`
- Output: `results/03_structural_prediction/step03_score_comparison.csv`, rank concordance plot
- Implement ipSAE (Dunbrack): PAE-filtered ipTM for disordered partners
- Implement LIS (PAE <= 12A interface-focused): better for syndecans, NRP1
- Rank all targets by ipTM, ipSAE, LIS; compute Kendall tau
- Acceptance:
  - [ ] VEGFA remains top-3 under all three metrics
  - [ ] CSV: target, ipTM, ipSAE, LIS, rank_ipTM, rank_ipSAE, rank_LIS
  - [ ] Kendall tau values reported

## Phase 2: Trimer and expanded structural predictions (~2 weeks GPU)

### T2.1 -- VEGF-bridged trimer predictions
NRP1 scored ipTM 0.22 in binary screen (expected false negative). Run 3-chain AF2:
- sFLT1-D1D3 + VEGF-A + NRP1 (b1/b2 domains, residues 275-586)
- sFLT1-D1D3 + VEGF-A + NRP2, sFLT1-D1D3 + VEGF-A + SEMA3A
- Acceptance: NRP1 trimer ipTM > 0.5 (recovers known interaction)
- Hypothesis link: hyp-028 (NRP1-dependent co-receptor sequestration)

### T2.2 -- sFLT1 D1-D7 full ectodomain screen
Current screen uses D1-D3 only (338 aa). Re-run top 5 hits + top 5 misses with full ectodomain (residues 1-750).
- Caveat: 750 + partner = >1000 residues, substantial GPU cost increase

## Phase 3: Biological validation integration

### T3.1 -- M2 Surfaceome Dual-Filter (~3h)
**Independent. Can start now.**

- Script: `analysis/03_structural_prediction/09_m2_surfaceome_filter.py`
- Output: `results/03_structural_prediction/step03_dual_filter_candidates.csv`, Venn diagram
- Download TMT surfaceome data: PRIDE PXD032801, PXD032823, PXD032967
- Cross-reference: AF2 structural hits (ipTM > 0.4) AND M2 macrophage surface expression
- Map against SomaScan panel for three-way overlap
- Acceptance:
  - [ ] Dual-filter table: target, ipTM, m2_surface_expressed, somascan_panel
  - [ ] Venn: AF2 hits / M2 surface / SomaScan overlap
- Hypothesis link: hyp-015, hyp-006

### T3.2 -- sFLT1-macrophage mediation analysis
**Blocked on cad-simoa-assembly completion.**
- SEM: sFLT1 -> sCD163 -> NfL path decomposition
- Control for NfL-myeloid reverse causation (2025 Cell Reports)
- Sex-stratified (sCD163 female-specific in late-stage PD)
- Hypothesis link: hyp-004, hyp-016

### T3.3 -- Brain PVM/Microglia Receptor Validation (~2h)
**Independent. Can start now.**

- Script: `analysis/03_structural_prediction/10_scrna_receptor_expression.py`
- Query SEA-AD + Allen Brain Cell Atlas for FLT1, NRP1, NRP2, SDC1, GPC1
- Cell types: PVMs, microglia, endothelial, pericytes
- Acceptance:
  - [ ] Expression matrix: gene x cell_type
  - [ ] FLT1 expression in PVMs confirmed or refuted

## Phase 4: Experimental validation planning

### T4.1 -- Crystallization propensity screening
XtalPred/PPCpred on top AF2 hits to triage: X-ray vs cryo-EM vs XL-MS/HDX-MS.

### T4.2 -- sFLT1 concentration-effect curve
Brain-specific dose-response characterization. Compare renal thresholds (preeclampsia) vs cerebrovascular context.

### T4.3 -- Heparanase-mediated local sFLT1 release
Mine CADASIL/CSVD proteomics for HPSE and sFLT1 co-elevation.

## Phase 5: SomaScan upstream refinements

### T5.1 -- ComBat-seq batch harmonization sensitivity
Re-run Step 1 with/without batch correction. Parameter sweep on covariates.

### T5.2 -- Negative-direction overlap analysis (~1h)
**Independent. Can start now.**
- Report 3-9 negative-direction consensus proteins
- If significant replication: flag as AF2 candidates (competitive displacement)

## Execution Order

| Priority | Task | Parallel? | Start condition | Est. effort |
|----------|------|-----------|-----------------|-------------|
| **0a** | **T0.1 Batch Correction Audit** | -- | **Now** | **1h** |
| **0b** | **T0.2 Comorbidity Check** | Parallel w/ T0.1 | **Now** | **30min** |
| **0c** | **T0.3 Sample Size Doc** | Parallel w/ T0.1 | **Now** | **30min** |
| **0d** | **T0.4 VEGF-Depletion Test** | Parallel w/ above | **Now** | **2h** |
| 1 | T1.2 Template Bias | -- | Now (18 targets) | 2h |
| 2 | T1.1 Domain Interface | Parallel w/ T1.2 | Now (18 targets) | 4h |
| 3 | T3.1 M2 Surfaceome | Parallel w/ above | Now | 3h |
| 4 | T3.3 scRNA Receptors | Parallel w/ above | Now | 2h |
| 5 | T5.2 Negative Overlap | Parallel w/ above | Now | 1h |
| 6 | T1.3 Alt Scores | After T1.1 | T1.1 complete | 3h |
| -- | Rerun T1.1-T1.3 | -- | All 24 complete (~Mar 23) | 1h |

**GATE:** If T0.1 reveals no batch correction in upstream LMs, promote T5.1 (ComBat-seq sensitivity) before final AF2 interpretation. Phase 1 analysis can proceed in parallel as preliminary, but results carry a "batch-unvalidated" caveat until T5.1 resolves.

**Total Phase 0: ~4h. Total Phase 1: ~16h dev + 1h rerun.**

## Dependency Graph

```
Phase 0 (upstream validation) ───────────────────────────┐
  0.1 Batch audit ─── GATE ──→ T5.1 if unharmonized     │
  0.2 Comorbidity check ──→ limitations flag             │
  0.3 Sample size doc  ──→ power assessment              │
  0.4 VEGF-depletion test ──→ mechanism viability gate   │
         │                                                │
         ▼ (if target list valid)                         │
Phase 1 (post-AF2, no new GPU) ──────────────────────────┤
  1.2 Template bias  ──→ Phase 3.1 (dual-filter)        │
  1.1 Domain interface ──→ Phase 2.2 (full ectodomain)   │
  1.3 Alternative scores ──→ (informs all downstream)    │
                                                          │
Phase 2 (new GPU jobs, ~2 weeks) ────────────────────────┤
  2.1 Trimer predictions ──→ Phase 4.1 (xtal screening)  │
  2.2 Full ectodomain    ──→ Phase 4.1                   │
                                                          │
Phase 3 (biological validation) ─────────────────────────┤
  3.1 M2 surfaceome   (unblocked)                        │
  3.2 Mediation        (blocked: cad-simoa-assembly)      │
  3.3 Brain scRNA-seq  (unblocked)                        │
                                                          │
Phase 4 (experimental, needs Phase 1-3) ─────────────────┤
Phase 5 (upstream refinement, independent) ───────────────┘
```

## GPU Budget

| Item | N jobs | Est. GPU-hours | Timeline |
|------|--------|---------------|----------|
| Current batch (6 remaining) | 6 | ~400 | Mar 19-23 |
| Trimers (Phase 2.1) | 3 | ~360 | 1 week |
| Full ectodomain (Phase 2.2) | 10 | ~500 | 1-2 weeks |
| **Total remaining** | **19** | **~1260** | **~3 weeks** |

All constrained to A100 (ADR-001).

## Risks and Assumptions

- **MSA directory structure:** verify AF2 `msas/` format on Minerva before T1.2
- **PRIDE data access:** TMT surfaceome may need manual download. Cache locally.
- **SEA-AD API stability:** fallback = pre-downloaded expression matrix
- **ipSAE/LIS:** no reference implementation; validate against VEGFA positive control
- **Domain boundaries:** D1-D3 from UniProt P17948 (FLT1_HUMAN). Verify against AF2 numbering (signal peptide offset?)

## Unresolved Questions

1. **[RESOLVED: NO]** Were upstream LMs batch-corrected? **No batch correction exists in the local pipeline. LMs were pre-computed upstream. Cross-cohort overlap uses raw significance lists.**
2. **[PARTIALLY RESOLVED]** MarkVCID: Age only (R script audit). GNPC, WASHU: **still needs PI**. UCSF_AD: age-adjusted.
3. **[PARTIALLY RESOLVED]** MarkVCID: **NO** -- HYPERTEN and Fazekas available but not used. Others: **still needs PI**.
4. **[PARTIALLY RESOLVED]** MarkVCID: **~69 post-exclusion** (confirmed from R script). UCSF_AD ~50-200, GNPC ~2,400+, WASHU ~2,200+ (still estimated). Confirm exact N with PI.
9. **[NEW]** Sex interactions exist for PLXNA1 and PLXNA4 (p=0.02) but not NRP1 (p=0.43) in MarkVCID. Should AF2 structural interpretation be sex-stratified for these targets?
10. **[NEW]** MarkVCID N=69 with Bonferroni across 7,335 proteins is underpowered. How many of our 24 AF2 targets survived only because of the large-N validation cohorts (GNPC/WASHU)?
5. **[RESOLVED: conservative]** 4-cohort overlap is conservative because small discovery cohorts (MarkVCID, UCSF_AD) gate the target list. Combined with |beta| > 0.5 effect-size filter, the design is more robust than initially feared. Remaining risk: confounders that replicate across all 4 cohorts.
6. AF2 MSA output format on Minerva -- exact dir structure for template hits in `msas/`? (Blocks T1.2)
7. sFLT1 numbering -- did FASTA include signal peptide (residues 1-26)? If yes, domain boundaries shift +26. (Blocks T1.1 correctness)
8. Negative overlap threshold for AF2 forwarding (T5.2) -- how many cohorts must agree? (User decision)
