# ADR-003: Upstream validation required before AF2 interpretation

Date: 2026-03-19
Status: accepted

## Context

Literature grounding review (2026-03-19) compared project methods against 40+ vault claims spanning SomaScan methodology, sFLT1 biology, AlphaFold scoring, and confounders. Code audit of `analysis/01_cross_cohort_overlap/` confirmed:

1. **No batch correction in the local pipeline.** The 24 AF2 targets were derived from cross-cohort overlap of pre-computed LM results received as CSVs from upstream collaborators. The R code (`01_load_helpers.R`, `02_compute_overlap.R`) loads these CSVs and computes set intersection. No ComBat-seq, no limma batch removal, no harmonization exists anywhere in the pipeline.

2. **Heterogeneous covariate adjustment across cohorts.** UCSF_AD LM files carry `_age` suffix (age-adjusted). MarkVCID and GNPC files have no suffix (unadjusted or unknown). WASHU LM files have no suffix but annotation files say `_adj_for_age`. The covariate specifications for MarkVCID, GNPC, and WASHU LMs are not documented in this repo.

3. **Massively different p-value magnitudes suggest huge sample size differences.** GNPC/WASHU show p-values in the 1e-170 to 1e-249 range; MarkVCID shows 1e-11 to 1e-14; UCSF_AD is intermediate. This implies N differs by orders of magnitude. The overlap approach (significant in ≥3 cohorts) is driven by the largest cohorts.

4. **Asymmetric negative direction.** GNPC has 3 negative-direction significant proteins vs MarkVCID's 318, making negative-direction overlap nearly empty (only 3 Tier 2 neg total).

5. **VEGF-depletion inversion untested.** The central mechanistic claim -- sFLT1 acts independently of VEGF -- has a direct counter-claim in the vault. Resolvable with mutual adjustment on existing data.

Vault evidence documents that batch harmonization alters sFLT1 candidate rankings (orientation claim) and that this alteration could be artifact rather than biology (untested inversion).

## Decision

Accept Phase 0 (upstream validation) in the post-AF2 analysis plan. Three items require PI input (covariate specs, sample sizes, batch correction rationale). Phase 1 analysis proceeds in parallel as preliminary with "upstream-unvalidated" caveat.

The 4-cohort significance overlap design may provide adequate protection against batch artifacts IF the cohorts were processed on different platforms or at different times (replication = implicit batch protection). But this assumption requires PI confirmation.

## Consequence

- AF2 results carry "upstream-unvalidated" caveat until PI clarifies covariate specs and sample sizes
- ComBat-seq sensitivity analysis (T5.1) conditionally promoted depending on PI response
- Does not block ongoing AF2 GPU jobs or Phase 1 preliminary analyses
- T0.4 (VEGF-depletion test) can proceed independently

## Evidence

- [[cross-cohort-somascan-batch-effect-harmonization-alters-the-ranking-and-significance-of-sflt1-protein-protein-interaction-candidates]]
- [[somascan-batch-harmonization-introduces-systematic-signal-distortion-that-changes-sflt1-association-rankings-by-artifact-rather-than-biology]]
- [[vascular-comorbidity-burden-confounds-sflt1-axonogenesis-pathway-associations-in-cross-sectional-somascan-cohorts]]
- [[somascan-multi-cohort-sample-size-is-unknown-pending-data-access-confirmation-across-all-four-cohorts]]
- [[sflt1-associated-reductions-in-axonogenesis-pathway-scores-reflect-vegf-pathway-suppression-rather-than-direct-sflt1-activity]]
- [[plasma sflt1 predicts axonal injury markers independently of vegf across multi-cohort somascan proteomics]]
- [[quantitative-assessment-of-pre-analytical-variability-sources-is-required-for-valid-batch-harmonization-in-multi-site-somascan-proteomics-studies]]
