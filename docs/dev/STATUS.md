# Status

Updated: 2026-03-21

## Now

- [ ] **AF2 test job running** -- af2_VEGFA (LSF 235867492), reduced_dbs fix, expect ~13:45
  - [x] Root cause: HHblits BFD titin column mismatch (not /tmp)
  - [x] Fix: `--db_preset=reduced_dbs` committed (153ec48)
  - [x] All 44 LSF scripts regenerated on Minerva
  - [x] Stale d1d3_corrected results backed up
  - [ ] If test passes: submit remaining 43 jobs (23 d1d3 + 11 d1d6 + 9 d1d7)

## Next

- Phase 1 analysis (blocked on AF2 results):
  - Template bias quantification (T1.2)
  - Domain-resolved interface mapping (T1.1)
  - Alternative scoring metrics (T1.3) -- ipSAE, LIS
- M2 surfaceome dual-filter (T3.1) -- independent, PRIDE data public
- Brain PVM/microglia receptor validation (T3.3) -- independent, SEA-AD public
- Negative-direction overlap analysis (T5.2) -- independent, upstream R re-analysis

## Later

- Trimer predictions: sFLT1 + VEGFA + NRP1/NRP2/SEMA3A (Phase 2, ~360 GPU-hrs)
- sFLT1-macrophage mediation analysis (blocked on cad-simoa-assembly)
- Crystallization propensity screening (needs Phase 1-3 results)
- See [post-af2-analysis.md](plans/post-af2-analysis.md) for full dependency graph and GPU budget

## Done

- [x] Step 1: Cross-cohort overlap -- 4 SomaScan cohorts, Tier 1+2 consensus (2026-02)
- [x] Step 2: Pathway enrichment -- ORA (GO/KEGG/Reactome), axon/semaphorin highlight (2026-02)
- [x] Step 3a: Candidate selection -- 24 targets (4 controls + 20 Tier 1) (2026-02)
- [x] Step 3b partial: 18/24 old D1-D3 runs (wrong construct, 1-338) (2026-03-19)
- [x] LSF path migration: SFLT1_VEGF -> vasc-sflt1-alphafold (2026-03-19)
- [x] **Construct correction**: signal peptide removed, 5T89 structural boundaries (2026-03-19)
- [x] **3-batch job generation**: D1-D3 corrected (24) + D1-D6 (11) + D1-D7 (9) = 44 jobs (2026-03-19)
- [x] Wrapper path bug fix + regeneration on Minerva (2026-03-19)
- [x] Phase 0 upstream validation COMPLETE (T0.1-T0.4) (2026-03-20)
- [x] **AF2 failure diagnosis**: HHblits BFD bug, reduced_dbs fix applied (2026-03-21)

## Key Results So Far (OLD CONSTRUCT -- will be superseded)

| Category | Target | ipTM | Notes |
|----------|--------|------|-------|
| Positive control | VEGFA | 0.816 | Pipeline calibrated (>0.7 gate) |
| High-confidence novel | PCDH9 | 0.806 | Strongest novel hit |
| High-confidence novel | NOE1 | 0.708 | |
| High-confidence novel | SLIK4 | 0.602 | |
| Expected false negative | NRP1 | 0.22 | VEGF-bridged ternary; binary screen blind spot |
| Expected false negative | NRP2 | 0.21 | Same mechanism as NRP1 |

*Scores above used construct 1-338 (includes signal peptide + D4 bleed). Corrected construct (27-330) scores pending -- blocked on AF2 test job.*
