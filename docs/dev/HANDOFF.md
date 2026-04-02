---
updated: "2026-04-02T15:40"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Fullscreen batch submission for all 3 sFLT1 constructs (D1-D3, D1-D6, D1-D7) with corrected topology pipeline. Added `--construct` flag to `07_generate_fullscreen_batch.py` so all constructs use one workflow. Regenerated FASTAs, LSF scripts, ran preflight, archived old results, and submitted to Minerva.

## Current State

### Submitted to Minerva (in progress)

| Construct | Targets | LSF Scripts | GPU-hours | Output dir |
|-----------|---------|-------------|-----------|------------|
| D1-D3 | 280 | 280 | 17,616 | `d1d3_fullscreen/` |
| D1-D6 | 280 | 280 | 22,704 | `d1d6_fullscreen/` |
| D1-D7 | 280 | 280 | 23,928 | `d1d7_fullscreen/` |
| **Total** | | **840** | **64,248** | |

- Submission was in progress at session end (`submit_all.sh` running for each construct)
- **Verify all 840 jobs landed**: `ssh minerva "bjobs 2>&1 | grep -c af2_"` — should be 840

### Topology filtering applied

- 364 Tier 1 consensus targets → 280 modeled (84 filtered)
- 30 dropped: lumenal-only TM (Golgi/ER interior, not extracellular; ADR-004)
- 7 dropped: ECD < 50 aa
- 47 dropped: intracellular-only, fetch failures, multi-accession complexes

### Archive

- Old results moved to `_archive_pre_fullscreen/` on Minerva
- Slimmed from 1.1 TB → 216 MB (kept scores.json, ranking_debug.json, ranked_0.pdb, timings)
- 12 legacy AF2 jobs killed before submission

### Pipeline code changes (committed + pushed)

- `07_generate_fullscreen_batch.py`: `--construct d1d3|d1d6|d1d7` flag, removed broken `AF_WALLTIME` monkey-patching
- `03_generate_lsf_jobs.py`: `WALLTIME_TIERS` and `WALLTIME_MAX` capped at 144h (gpu queue hard limit is 8640 min)
- `00_uniprot_api.R`: added `show_all()` helper
- Preflight report updated (4/4 passed)

### Walltime fix (applied on Minerva, not yet committed locally)

- 9 jobs (3 per construct) had 192h walltimes exceeding the 144h gpu queue limit
- Patched in-place on Minerva via sed: `192:00 → 144:00`
- Source fix in `03_generate_lsf_jobs.py` done locally but not yet pushed

## Next Steps

1. **Verify submission** — confirm 840 jobs on queue. If `submit_all.sh` stalled, resubmit the remaining construct(s)
2. **Commit walltime fix** — `git add analysis/03_structural_prediction/03_generate_lsf_jobs.py && git commit && git push`
3. **Monitor batch** — `ssh minerva "bjobs 2>&1 | grep af2_ | awk '{print \$3}' | sort | uniq -c"` for PEND/RUN/DONE counts
4. **After completion** — `bash 11_batch_extract_and_merge.sh --force` on each construct to extract scores with new schema
5. **Positive controls missing** — VEGFA/PlGF are not in the consensus target list (they're the ligand, not SomaScan-associated). Prior d1d3 batch had them as manual additions. Consider whether to add them back as controls for the fullscreen batches
6. **Multi-accession targets** — `UBE2N/UB2V1 Complex.1` (P61088|Q13404) failed because UniProt API doesn't accept pipe-delimited accessions. Needs a complex-aware fetch if these matter
7. **sLRP1/megalin OOM risk** — capped at 144h walltime, but may still OOM at 128 GB. Check these targets first after results come in

## Key Decisions

- **One script, all constructs.** `07_generate_fullscreen_batch.py --construct X --include-existing` is the canonical workflow going forward
- **144h walltime cap.** Minerva gpu queue hard limit. Targets needing >144h (megalin-class, >4000 aa total) will need checkpoint-restart or splitting
- **Archive, don't delete.** Old results slimmed to lightweight archive preserving scores + best model
- **No intracellular filter.** Kept from prior session — SomaScan detected these in CSF, modeling is informative

## Open Questions

- Submission may not have completed for all 3 constructs — verify first thing
- Age adjustment asymmetry across upstream cohorts (carried forward)
- 18 no-keyword proteins with unknown localization — kept, review after results
