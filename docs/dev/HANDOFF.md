---
updated: "2026-03-20T10:30"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Two things this session: (1) added ectodomain truncation for transmembrane AF2 partners, (2) ran the T0.4 VEGF-depletion inversion test. Also verified all protein sequences against UniProt and submitted all 44 AF2 jobs.

## Current State

- **All 44 AF2 jobs COMPLETED** on Minerva across 3 batches (d1d3_corrected x24, d1d6 x11, d1d7 x9). All using ectodomain-only sequences for TM proteins.
- **Phase 0 upstream validation COMPLETE** (T0.1-T0.4 all resolved)
- **T0.4 result:** 37/38 axon/guidance enrichment terms retained after removing 7 VEGF-pathway genes. Signal is VEGF-independent. Inversion hypothesis not supported.
- **Ectodomain filter** added to `02_fetch_sequences.py`. 16 TM proteins truncated, 5,447 aa total reduction (23.7%).
- **All sequences verified** against live UniProt -- 100% match for all 24 candidates.
- **28/28 tests passing.**

## Next Steps

1. **Phase 1 analysis (UNBLOCKED):**
   - T1.2: Template bias quantification -- parse AF2 `msas/` for PDB template hits, correlate with ipTM
   - T1.1: Domain-resolved interface mapping -- per-domain (D1/D2/D3) PAE scores, classify binding modes
   - T1.3: Alternative scoring (ipSAE, LIS) -- depends on T1.1
2. **Compare old construct (1-338) vs corrected (27-330) scores** -- quantify signal peptide effect
3. **Compare D1-D3 vs D1-D6 vs D1-D7 scores** -- construct length effect
4. **Independent tasks:**
   - T3.1: M2 surfaceome dual-filter (PRIDE data)
   - T3.3: Brain PVM/microglia receptor validation (SEA-AD)
   - T5.2: Negative-direction overlap analysis
5. **Ask PI:** covariates in GNPC, WASHU, UCSF_AD upstream LMs

## Key Decisions

- Ectodomain-only modeling for TM proteins: `fetch_all_sequences(apply_ectodomain_filter=True)` is now default. Min ectodomain length: 50 aa.
- T0.4 VEGF-depletion: used broad definition (VEGF ligands + receptors + NRPs + PDGF + FGF + angiopoietins). Only 7 of 25 defined genes were in the consensus list.
- Phase 0 gate passed: target list is valid for interpretation. Remaining caveat: vascular comorbidity confounding not adjusted in any cohort (documented limitation).

## Open Questions

- Sex-stratified AF2 interpretation for PLXNA1/PLXNA4 (sFLT1 × Sex interaction p=0.02)?
- Negative-direction overlap threshold: relax to 3/4 cohorts? (GNPC has only 3 neg proteins)
- Confirm exact N and covariates with PI for GNPC, WASHU, UCSF_AD

## Files Modified This Session

- `analysis/03_structural_prediction/02_fetch_sequences.py` -- ProteinTopology, ectodomain truncation
- `analysis/03_structural_prediction/test_structural.py` -- TestProteinTopology
- `analysis/02_pathway_enrichment/run_vegf_depletion.R` -- NEW: T0.4 inversion test
- `docs/dev/STATUS.md` -- T0.4 marked complete
- `docs/dev/HANDOFF.md` -- this file
