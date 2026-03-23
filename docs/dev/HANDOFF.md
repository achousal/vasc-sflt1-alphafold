---
updated: "2026-03-23T15:15"
project: "vasc-sflt1-alphafold"
---

## What I Was Doing

Full pipeline audit, signal peptide fix, preflight check system, clean slate production submission.

## Current State

- **24 d1d3 jobs running on Minerva** (21 RUN, 3 PEND). Job IDs 236035782 (VEGFA) + 236035796-236035819. Expected completion: 48-144h depending on target size.
- **VEGFA test job validated** -- MSA phase running, chain_id_map.json confirms correct sequences (sFLT1 starts Ser27, 304 aa; VEGFA 395 aa full-length).
- **d1d6 and d1d7 not yet submitted** -- ready to go, all FASTA/LSF regenerated.
- **Stale sentinel in `d1d3/logs/sentinels/completed.log`** has one premature `af2_VEGFA` entry from a prior failed run. Clear before checking final completion counts.

## Session Summary

1. **Material balance doc** (`docs/dev/material-balance.md`) -- full pipeline data flow from upstream LR equation through AF2 submission.
2. **Signal peptide bug found and fixed** in `02_fetch_sequences.py` -- 5 secreted/GPI-anchored targets had signal peptides included. Fixed with `chain`, `gpi_anchor` fields in ProteinTopology.
3. **Preflight check system** (`analysis/checks/`) -- 3 scripts + runner:
   - `check_fasta_integrity.py` -- sequence lengths, headers, signal peptides, UniProt verification (default on), bait-target inventory table
   - `check_candidates_complete.py` -- FASTA/LSF/manifest completeness, orphan detection
   - `check_af2_inputs.py` -- decodes LSF base64, verifies AF2 flags, databases, paths
   - `run_preflight.sh` -- runs all 3, writes `preflight_report.txt`
4. **Renamed `d1d3_corrected` → `d1d3`** everywhere (code, docs, Minerva).
5. **Removed `--run_relax`** (not supported by Minerva container).
6. **Added `--use_gpu_relax=false`** (required by Minerva AF2.3.2).
7. **Archived old logs** to `logs_archive_pre_20260323/` on Minerva.
8. **Deleted all stale artifacts** -- root fasta/, root jobs/, orphan FASTA in d1d6/d1d7, old docs.

## Next Steps

1. **Monitor d1d3 jobs** -- `ssh minerva "bjobs -w | grep af2_"`
2. **Clear stale sentinel** -- `ssh minerva "echo -n > .../d1d3/logs/sentinels/completed.log"` before checking final counts
3. **Submit d1d6** (11 jobs) and **d1d7** (9 jobs) when ready -- no dependency on d1d3 finishing
4. **Phase 1 analysis** (blocked on AF2 results): `04_parse_results.py` to extract ipTM/PAE scores

## Key Decisions

- Single test job before batch submission catches container flag issues early.
- `--use_gpu_relax=false` is the correct flag for Minerva's AF2.3.2 singularity build. `--run_relax` is not recognized.
- Preflight checks run UniProt sequence verification by default (`--skip-uniprot` for offline).
- Batch name is `d1d3` not `d1d3_corrected` -- it's the canonical batch.

## Open Questions

- Age adjustment asymmetry across cohorts (carried forward)
- Confirm exact N and covariates with PI for GNPC, WASHU, UCSF_AD (carried forward)
- STMN3 and Calcineurin B a are cytoplasmic -- expect low ipTM as implicit negative controls
