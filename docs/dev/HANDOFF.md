---
updated: "2026-04-02T16:00"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Full audit and correction of the AF2 target preparation pipeline before resubmitting the fullscreen batch. Upgraded scoring extraction (per-prediction metrics, ipSAE, LIS, template coverage), rewrote sequence fetch with UniProt-grounded topology, and resolved 12 decision points for resubmission.

## Current State

### Pipeline code updated (6 scripts)

- **02_fetch_sequences.py** — rewritten. Uses annotated UniProt `Topological domain: Extracellular` as primary ECD source, inferred SP->TM as fallback. Drops lumenal-only TM proteins. Takes largest single extracellular segment for multi-pass TM. Removes signal peptides from soluble/secreted.
- **03_generate_lsf_jobs.py** — memory tier (32 GB/core for > 1600 aa total), GPU constraint updated to V100 || A100 (exclude H100).
- **09_extract_and_cleanup.py** — extracts ipSAE, LIS, template coverage, per-prediction PAE/pLDDT from ALL predictions. Retains features.pkl.
- **10_merge_scores.py** — reads new nested scores.json schema (summary + per_prediction).
- **04_parse_results.py** — scores.json as primary source, pkl fallback.
- **05_plot_results.py** — ipSAE rescue candidates, template bias warnings in reports.
- **00_uniprot_api.R** — self-contained UniProt explorer for interactive grounding in RStudio.

### Target disposition (from topology audit)

| Category | Count |
|---|---|
| Total candidates | 365 |
| Dropped (19 lumenal-only + 6 tiny ECD no annotation) | 25 |
| To model | 340 |
| Keep existing results (FASTA unchanged) | 68 |
| Re-run (FASTA changed by corrections) | 95 |
| Never ran | 166 |
| **GPU jobs needed** | **261 (d1d3)** |

d1d6 and d1d7 also need regeneration with same corrections -- all three constructs submitting simultaneously.

### Decisions documented (ADR-004)

See `docs/dev/decisions/004-target-topology-filters.md` for full rationale:
1. Keep intracellular proteins (negative controls + dead-cell-leaking hypothesis)
2. Drop lumenal-only TM (Golgi/ER interior != extracellular)
3. Largest single extracellular segment for multi-pass TM
4. 32 GB/core memory for complexes > 1600 aa total
5. V100 + A100 allowed (ADR-001 updated)

### 25 running d1d3 jobs

Still running from prior batch. NOMO2 will timeout (only 12/25 predictions done at 24h wall). NEUM and NPTN should finish. These used old FASTAs -- results will be superseded by the corrected re-run for any target whose sequence changed.

### Scoring pipeline ready

New `scores.json` schema extracts per-prediction: ipTM, pTM, interchain PAE, interface pLDDT, ipSAE, LIS. Summary includes mean+/-std for all metrics + template coverage from features.pkl + template_bias_flag. `11_batch_extract_and_merge.sh` has `--force` flag to re-extract with new schema.

## Next Steps

1. **Regenerate FASTAs** -- run `02_fetch_sequences.py` for all 3 constructs (d1d3, d1d6, d1d7) with corrected topology. ~2h for 365 UniProt API calls per construct.
2. **Regenerate LSF scripts** -- run `03_generate_lsf_jobs.py` with new walltimes + memory tiers.
3. **Clear stale sentinel** -- `echo -n > .../d1d3/logs/sentinels/completed.log`
4. **Preflight check** -- run `analysis/checks/run_preflight.sh` to validate all FASTAs, LSF scripts, manifest.
5. **Submit** -- orchestrator handles chunked submission. 261 d1d3 + d1d6 + d1d7 jobs.
6. **After completion** -- run `11_batch_extract_and_merge.sh --force` to extract with new schema for all targets (including the 68 kept from prior run that have old-schema scores.json).

## Key Decisions

- **No intracellular filter.** SomaScan detected these in CSF -- modeling them is informative regardless of canonical localization.
- **Lumenal != extracellular.** Golgi/ER lumenal domains cannot interact with extracellular sFLT1. 19 targets dropped.
- **Annotated ECD over inferred.** UniProt `Topological domain: Extracellular` is ground truth. Fixes type II TM proteins (TWEAK, AT1B1, AT1B2) that the old SP->TM inference got wrong.
- **Largest single segment for multi-pass TM.** Union span includes TM helices and cytoplasmic loops -- AF2 would try to fold them.
- **V100 validated.** 181 successful runs on V100 nodes. Pool is larger than A100-only.

## Open Questions

- sLRP1 and megalin (~4700 aa total) may still OOM at 128 GB. Will find out on submission.
- 18 no-keyword proteins have unknown localization. Kept for now; review after results.
- d1d6/d1d7 target lists: apply same 340-target filter, or subset to d1d3 hits only? (Decided: same full set)
- Age adjustment asymmetry across upstream cohorts (carried forward from prior sessions)
