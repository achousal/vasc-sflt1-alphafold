---
updated: "2026-03-21T12:05"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Starting Phase 1 analysis. Discovered all 44 AF2 jobs failed silently -- HHblits crashed on chain B MSA search but the wrapper's EXIT trap logged them as "completed". Zero model outputs exist across all 3 batches.

## Current State

- **All 44 AF2 jobs FAILED** -- zero `ranking_debug.json`, zero `result_model_*.pkl`. Only chain A (sFLT1) MSAs exist. Chain B HHblits crashed after ~28 min in every job.
- **Root cause hypothesis:** `/tmp` exhaustion on shared GPU nodes. Multiple AF jobs writing large HHblits temp files to the same node `/tmp`.
- **Diagnostic job submitted:** `diag_VEGFA` (LSF 235865567) on `lg07c01`, started 11:59. Fix: `TMPDIR` set to project scratch + `--use_precomputed_msas`. Expect result by ~12:50.
- **Wrapper bug fixed:** EXIT trap now routes to `failed.log` (rc!=0) vs `completed.log` (rc=0).
- **sflt1_length corrected:** D1-D3=304aa (was 338), D1-D6=631, D1-D7=721. `--construct` flag added to `run_structural.py`.
- **Domain boundaries confirmed** for PAE slicing (0-indexed in 304-aa FASTA): D1=5-103, D2=105-198, D3=199-303.

## Next Steps

1. **Check diagnostic job** -- `ssh minerva "bjobs 235865567"` then inspect output log
   - **If SUCCESS:** Regenerate all 44 LSF scripts with TMPDIR fix, resubmit in 3 batches. Build Phase 1 modules against real output structure.
   - **If FAIL (same HHblits error):** Try `--db_preset=reduced_dbs` (skip BFD, use only UniRef30) as fallback. Smaller search = less /tmp pressure. Or request dedicated GPU node.
2. **Build Phase 1 modules** (blocked on AF2 results):
   - `07_template_bias.py` (T1.2) -- PDB template count vs ipTM correlation
   - `08_domain_resolved.py` (T1.1) -- per-domain PAE submatrices, binding mode classification
   - `09_alternative_scoring.py` (T1.3) -- ipSAE, LIS metrics
3. **Tests for Phase 1** -- mock PAE matrices, verify VEGFA->D2 mapping

## Key Decisions

- Build Phase 1 analysis code only after real AF2 outputs exist (avoid rework from guessing output structure).
- Diagnostic before batch resubmit (option b) to avoid burning 44 more GPU-hours on systematic failure.
- TMPDIR to project scratch as primary fix hypothesis for HHblits /tmp exhaustion.

## Open Questions

- Is /tmp exhaustion the actual cause, or is it a BFD database issue?
- If TMPDIR fix works, should we regenerate all 44 LSF scripts from `03_generate_lsf_jobs.py` (adds TMPDIR globally) or patch the existing scripts?
- Sex-stratified AF2 interpretation for PLXNA1/PLXNA4 (carried forward)
- Confirm exact N and covariates with PI for GNPC, WASHU, UCSF_AD (carried forward)

## Files Modified This Session

- `analysis/03_structural_prediction/03_generate_lsf_jobs.py` -- wrapper EXIT trap fix (failed.log routing)
- `analysis/03_structural_prediction/04_parse_results.py` -- sflt1_length default 338->304
- `analysis/03_structural_prediction/run_structural.py` -- `--construct` flag, construct-aware sflt1_length
- `results/.../d1d3_corrected/jobs/diag_VEGFA.lsf` -- NEW: diagnostic job (gitignored)
