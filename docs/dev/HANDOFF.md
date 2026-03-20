---
updated: "2026-03-20T10:05"
---

## What I Was Doing

Added ectodomain truncation for transmembrane proteins in AlphaFold Multimer FASTA generation. 16 of 24 candidates are transmembrane -- previously modeled full-length (including cytoplasmic domains that cannot interact with extracellular sFLT1). Now uses UniProt topology annotations to extract the full extracellular domain only.

## Current State

- **Ectodomain filter implemented** in `02_fetch_sequences.py`: fetches UniProt JSON topology, extracts signal peptide and TM helix positions, truncates to ectodomain for TM proteins
- **All 44 AF2 jobs regenerated** with ectodomain-filtered partner sequences (3 batches: d1d3_corrected x24, d1d6 x11, d1d7 x9)
- **5,447 aa total reduction** (23.7%) across all candidates -- reduces GPU cost and eliminates biologically impossible interface predictions
- **All sequences verified** against live UniProt -- 100% match
- **28/28 tests passing** including new topology tests
- **44 AF2 jobs NOT YET SUBMITTED** -- need git push + minerva pull first
- **T0.4 (VEGF-depletion inversion test): NOT STARTED**

## Next Steps

1. **Commit, push, and deploy to Minerva** (`git push` + `ssh minerva "cd ... && git pull"`)
2. **Submit 44 AF2 jobs** on Minerva (`bash results/03_structural_prediction/d1d3_corrected/jobs/submit_all.sh`)
3. **Ask PI:** covariates in GNPC, WASHU, UCSF_AD upstream LMs
4. **T0.4:** Run VEGF-depletion inversion test
5. Start T1.2 (template bias) and T1.1 (domain interface) once D1-D3 corrected batch completes

## Key Decisions

- Ectodomain-only modeling for TM proteins: biologically required because sFLT1 is extracellular. Guardrail was already in CLAUDE.md but not implemented in code.
- FASTA headers now include `ectodomain_X-Y` provenance for truncated proteins
- `fetch_all_sequences()` gains `apply_ectodomain_filter` parameter (default True) for backwards compat
- Minimum ectodomain length threshold: 50 aa (below this, fall back to full-length)

## Files Modified This Session

- `analysis/03_structural_prediction/02_fetch_sequences.py` -- ProteinTopology dataclass, `fetch_topology_from_uniprot()`, ectodomain truncation in `fetch_all_sequences()`
- `analysis/03_structural_prediction/test_structural.py` -- TestProteinTopology class, partner_region FASTA test
- `results/03_structural_prediction/d1d3_corrected/` -- regenerated FASTAs and LSF jobs
- `results/03_structural_prediction/d1d6/` -- regenerated FASTAs and LSF jobs
- `results/03_structural_prediction/d1d7/` -- regenerated FASTAs and LSF jobs
