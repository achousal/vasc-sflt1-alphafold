---
updated: "2026-04-01T13:15"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Batch score extraction and analysis of the full D1-D3 AF2 screen. Ran `09_extract_and_cleanup.py` on all completed targets on Minerva, merged into a unified CSV via `10_merge_scores.py`, and updated the progress report with 148-target results.

## Current State

### D1-D3 screen: 148/365 scored

- Orchestrator job 236397016 still managing remaining 217 targets
- All 148 completed targets now have `scores.json` sidecars on Minerva
- Unified scores CSV: `d1d3/step03_interaction_scores_all.csv` (Minerva + local)
- 62 result directories exist but have no `ranking_debug.json` (in-progress or failed)

### Key finding: VEGFA is the only consistent hit

- 30/148 targets cross ipTM > 0.6, but **only VEGFA has spread < 0.10** (all 25 models agree)
- 4 targets get "predicted" call (DDAH1, CAH11, KBRS1, CNBP1) — all fail consistency filter (spread > 0.20)
- 20% apparent FPR at ipTM > 0.6, consistent with published AF2 multimer benchmarks
- No novel sFLT1 binary interaction with axon guidance proteins detected in D1-D3

### Cross-construct signal (unchanged from prior session)

- NRP1: 0.25 (D1-D3) → 0.55 (D1-D6) → 0.60 (D1-D7) — gains with longer construct
- NGL1: 0.26 → 0.25 → 0.61 — gains with D1-D7
- Consistency analysis not yet applied to D1-D6/D1-D7

### Code changes this session

- `10_merge_scores.py` — added `safe_name_alt()` fallback for Minerva `/` -> `_` convention
- `11_batch_extract_and_merge.sh` — new script for batch extraction + merge on Minerva (login node, no GPU)
- Both scripts synced to Minerva

### Artifacts

- `d1d3/step03_interaction_scores_all.csv` — 365 rows (148 scored, 217 `not_run`)
- `docs/dev/progress-report-2026-04-01.md` — full progress report with 148-target analysis
- ~1.3 GB pkl files freed on Minerva during extraction

### Caveat: pilot 24 targets have `no_data` PAE

The original 24 pilot targets were extracted this session but their pkl files only contained `ranking_debug.json` scores (best-model pkl had already been cleaned up in a prior run, leaving only 1 pkl per target). ipTM scores are correct from `ranking_debug.json` but inter-chain PAE and interface pLDDT are missing (`no_data`). Not blocking — consistency filter works on ipTM alone.

## Next Steps

1. **Re-run `11_batch_extract_and_merge.sh` as orchestrator completes batches** — script is idempotent, skips already-extracted targets
2. **Run FLT1 PVM gate** — download SEA-AD data, submit LSF job. Highest-leverage experiment
3. **Complete D1-D6/D1-D7** and apply consistency filter to cross-construct results
4. **Apply ipSAE/LIS alternative scoring** once full screen completes — may rescue candidates with template bias
5. **Decide on 3 complex targets and 8 OOM targets** (carried forward)

## Key Decisions

- **Minerva HPC path is `d1d3/` not `d1d3_corrected/`.** The local `batch_summary.json` references `d1d3_corrected` but Minerva filesystem uses `d1d3`. The batch script uses the Minerva path.
- **Consistency filter (ipTM spread) confirmed as primary hit classifier.** At 148 targets, only VEGFA passes. This is the methodological contribution of the screen.
- **Progress report framed as "148 of 365 scored, here's what we see"** — not as pilot-then-expansion narrative.

## Open Questions

- 8 OOM kills at 64 GB (megalin, DSCAM, etc.) — need higher memory or construct truncation
- 3 multi-subunit complex targets: model as 3-chain or drop?
- Should Tier 2 candidates be added given 0% consistent hit rate in Tier 1?
