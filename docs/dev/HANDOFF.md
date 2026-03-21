---
updated: "2026-03-21T13:00"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Diagnosing and fixing systematic AF2 job failures. All 44 jobs failed due to HHblits BFD bug, not /tmp exhaustion. Applied reduced_dbs fix, regenerated all scripts, submitted test job.

## Current State

- **Root cause confirmed:** HHblits crashes on BFD titin sequences (`hhalignment.cpp:1244: Compress: sequences in merged A3M file do not all have the same number of columns`). NOT /tmp exhaustion (node had 710G free).
- **Fix applied:** `--db_preset=reduced_dbs` (skips full BFD, uses UniRef30 + small_bfd). Committed and pushed (`153ec48`).
- **All 44 LSF scripts regenerated** on Minerva for all 3 batches (d1d3_corrected: 24, d1d6: 11, d1d7: 9).
- **Stale results cleared:** d1d3_corrected old MSAs (wrong construct 1-338) moved to `results_old_construct_338/`. Sentinels cleared.
- **Test job submitted:** af2_VEGFA (LSF 235867492), d1d3_corrected batch, reduced_dbs. Expect result ~13:45.
- **Wrapper bug fixed** (previous session): EXIT trap routes to `failed.log` (rc!=0) vs `completed.log` (rc=0).
- **Construct lengths confirmed:** D1-D3=304aa, D1-D6=631aa, D1-D7=721aa.
- **Domain boundaries confirmed** for PAE slicing (0-indexed): D1=5-103, D2=105-198, D3=199-303.

## Next Steps

1. **Check test job** -- `ssh minerva "bjobs 235867492"` then inspect output
   - **If SUCCESS:** Submit remaining 23 d1d3_corrected jobs, then d1d6 (11) and d1d7 (9) batches. Build Phase 1 modules.
   - **If FAIL:** Escalate to HPC support or try AF2.3.2 with `--db_preset=reduced_dbs --use_precomputed_msas` after manually running jackhmmer for chain B.
2. **Build Phase 1 modules** (blocked on AF2 results):
   - `07_template_bias.py` (T1.2) -- PDB template count vs ipTM correlation
   - `08_domain_resolved.py` (T1.1) -- per-domain PAE submatrices, binding mode classification
   - `09_alternative_scoring.py` (T1.3) -- ipSAE, LIS metrics
3. **Tests for Phase 1** -- mock PAE matrices, verify VEGFA->D2 mapping

## Key Decisions

- Build Phase 1 analysis code only after real AF2 outputs exist (avoid rework from guessing output structure).
- Test single job before batch submission to avoid burning GPU-hours on systematic failure.
- `reduced_dbs` is acceptable for human protein multimer PPI -- quality difference vs full_dbs is marginal for well-represented sequences.
- Stale d1d3_corrected results (wrong construct MSAs) backed up rather than deleted.

## Open Questions

- Sex-stratified AF2 interpretation for PLXNA1/PLXNA4 (carried forward)
- Confirm exact N and covariates with PI for GNPC, WASHU, UCSF_AD (carried forward)

## Files Modified This Session

- `analysis/03_structural_prediction/03_generate_lsf_jobs.py` -- `full_dbs` -> `reduced_dbs`
