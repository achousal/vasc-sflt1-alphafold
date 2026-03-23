---
updated: "2026-03-22T21:30"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Material balance audit of the full pipeline (LR → overlap → enrichment → AF2). Discovered and fixed signal peptide inclusion in 5 target FASTA files. Regenerated all 44 LSF jobs across 3 constructs.

## Current State

- **material-balance.md written** (`docs/dev/material-balance.md`) -- documents the full data flow from upstream LR equation (`protein_i ~ VEGFsR1 + Age`, log10+z-scored, Bonferroni, |beta|>=0.5) through cross-cohort overlap (678 Tier 1+2 proteins) → pathway enrichment (axonogenesis padj=8.21e-52) → candidate selection (24 targets) → AF2 structural prediction on Minerva.
- **Signal peptide bug found and fixed** in `02_fetch_sequences.py`. The ectodomain filter only activated for transmembrane proteins. Secreted proteins with signal peptides (but no TM domain) got full-length sequences including their signal peptides. Same class of error previously corrected for sFLT1.
- **5 targets affected:** SEMA3A (-20 aa), Contactin-5 (-46 aa, signal + GPI tail), FSTL4 (-22 aa), NOE1 (-16 aa), SLIT2 (-30 aa).
- **Code fix:** Added `chain`, `gpi_anchor` fields to `ProteinTopology`. GPI-anchored proteins now use UniProt Chain boundaries. Soluble proteins with signal peptides get signal peptide trimmed.
- **All FASTA files regenerated** across d1d3 (24), d1d6 (11), d1d7 (9).
- **All 44 LSF jobs regenerated** via `06_generate_ec_batches.py` with correct per-target walltimes. Walltimes stayed in same bins (deltas too small to cross thresholds).
- **19 previously-correct targets verified unchanged** (14 ectodomain-extracted + 3 legitimately full-length: VEGFA, STMN3, Calcineurin B a).

## Next Steps

1. **Commit and push** the code fix + regenerated files
2. **`git pull` on Minerva** and resubmit AF2 jobs
   - If VEGFA test job (LSF 235867492) succeeded: submit all 3 batches
   - If not: check result first, then batch submit
3. **Update material-balance.md** FASTA inventory table and audit section to reflect corrected values (the doc was written before the fix)
4. **Phase 1 analysis modules** (blocked on AF2 results):
   - `07_template_bias.py`, `08_domain_resolved.py`, `09_alternative_scoring.py`

## Key Decisions

- Signal peptide removal applies to ALL secreted/surface proteins, not just sFLT1. The localization filter now handles 3 cases: transmembrane (ectodomain extraction), GPI-anchored (Chain boundaries), soluble with signal peptide (trim signal).
- The enrichment step (Step 2) is not a filtering dependency -- all 20 data-driven candidates are Tier 1 + dual-somamer + in axon pathway regardless. It serves scientific justification, not analytical necessity.
- d1d6/d1d7 walltimes were fixed from flat 24h to proper per-target scaling (72-144h).

## Open Questions

- Age adjustment asymmetry across cohorts (MarkVCID/GNPC: no age adj in data filenames; UCSF_AD/WASHU: `_age` suffix). The LR equation in the MarkVCID script already includes Age -- unclear what "age-adjusted" means for the other cohorts.
- Annotation CSV schemas differ across cohorts (column names, column counts). Need to verify the overlap code handles this correctly.
- Sex-stratified AF2 interpretation for PLXNA1/PLXNA4 (carried forward)
- Confirm exact N and covariates with PI for GNPC, WASHU, UCSF_AD (carried forward)

## Files Modified This Session

- `analysis/03_structural_prediction/02_fetch_sequences.py` -- signal peptide fix for soluble + GPI-anchored proteins
- `docs/dev/material-balance.md` -- new file, full pipeline material balance
- `results/03_structural_prediction/d1d3/` -- regenerated FASTA + LSF jobs
- `results/03_structural_prediction/d1d6/` -- regenerated FASTA + LSF jobs
- `results/03_structural_prediction/d1d7/` -- regenerated FASTA + LSF jobs
