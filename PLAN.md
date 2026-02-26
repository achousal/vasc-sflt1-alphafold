# SFLT1_VEGF Analysis Plan

Executable plan for the 3-step analysis pipeline. Read CLAUDE.md for project context.

---

## Data Schema Notes (verified from files)

Row-index column presence (leading `""` column):
- MarkVCID, GNPC: **have** row-index column
- UCSF_AD, WASHU: **no** row-index column

UCSF_AD files have `_age` suffix (e.g., `values_LM_VEGFsR1_sig_pos_age.csv`).

UCSF_AD annotation files have extra columns (SeqId, SomaId) -- drop during normalization.

GNPC annotation uses snake_case headers (`target_full_name`, `uni_prot`, `entrez_gene_id`); others use CamelCase (`TargetFullName`, `UniProt`, `EntrezGeneID`).

Neg annotation files for UCSF_AD and WASHU are single-column (header `"x"`) with just protein names. GNPC has no neg annotation file. MarkVCID neg annotation is full multi-column.

Negative-direction overlap will likely be near-zero: GNPC neg has 3 proteins (VEGFsR1) and 5 (VEGFsR1.1); WASHU neg has 5 and 9. Handle empty consensus lists gracefully.

---

## Step 1: Cross-Cohort Overlap (R, local)

### Objective
Identify proteins with replicated sFLT1 associations across cohorts, stratified by somamer and direction.

### Input
All 16 LM result CSVs under `data/{cohort}/` plus annotation CSVs.

### Scripts (in `analysis/01_cross_cohort_overlap/`)

**`01_load_helpers.R`** -- Shared functions:
- `load_lm_results(path, cohort)`: read CSV, detect/drop row-index column (check if first column header is `""`), add `cohort` column. Return: `Target, beta, pvalue, padj, minus_log_padj, direction, cohort`.
- `load_annotation(path, cohort)`: read CSV, normalize to snake_case (`Target` -> `target`, `UniProt` -> `uniprot`, `EntrezGeneSymbol` -> `entrez_gene_symbol`). Handle single-column neg files (header `"x"`) -> return tibble with `target` column only. Drop SeqId/SomaId extras.
- `build_cohort_registry()`: return tibble mapping `cohort, role, somamer, direction, file_path` for all 16 LM files. Handle UCSF_AD `_age` suffix programmatically.

**`02_compute_overlap.R`** -- Core overlap logic:

Replication rule per (somamer, direction):
- **Tier 1:** significant in all 4 cohorts
- **Tier 2:** both discovery (MarkVCID AND UCSF_AD) + at least 1 validation (GNPC or WASHU)
- **Tier 3:** both discovery only (no validation) -- reported but NOT forwarded to Steps 2-3

Cross-somamer flag: proteins replicating in BOTH VEGFsR1 and VEGFsR1.1 marked `dual_somamer = TRUE`.

Consensus list schema:
```
target, uniprot, entrez_gene_symbol, entrez_gene_id, direction, somamer,
tier, n_cohorts, cohorts_sig, mean_beta, min_padj, dual_somamer
```

Annotation join priority: MarkVCID -> UCSF_AD -> GNPC -> WASHU. Fuzzy match on `.N` suffix stripping if exact match fails.

**`03_plot_overlap.R`** -- Visualizations:
- 4 UpSet plots (one per somamer x direction) via `UpSetR` or `ComplexUpset`, replicated intersections highlighted
- Beta heatmap: rows = Tier 1+2 proteins, columns = cohorts, fill = beta, faceted by somamer
- Tier summary bar chart: counts per tier per somamer-direction

**`run_overlap.R`** -- Entrypoint. Sources 01-03. CLI args via `optparse`: `--data-dir` (default `data/`), `--results-dir` (default `results/01_cross_cohort_overlap/`).

### R packages
`tidyverse`, `optparse`, `UpSetR` or `ComplexUpset`, `ggplot2`

### Output (`results/01_cross_cohort_overlap/`)
```
step01_consensus_proteins_pos.csv      -- Tier 1+2 positive, both somamers
step01_consensus_proteins_neg.csv      -- Tier 1+2 negative, both somamers
step01_all_overlaps.csv                -- Full table all tiers
step01_upset_VEGFsR1_pos.pdf
step01_upset_VEGFsR1_neg.pdf
step01_upset_VEGFsR1.1_pos.pdf
step01_upset_VEGFsR1.1_neg.pdf
step01_beta_heatmap.pdf
step01_tier_summary.pdf
step01_stats_report.txt                -- Counts per tier/direction/somamer, annotation join rate
README.md
```

### Test: `test_overlap.R`
- `load_lm_results()` handles files with and without row-index columns
- `load_annotation()` normalizes across all 4 cohorts, handles single-column neg files
- Overlap logic: toy 4-cohort scenario -> verify tier assignment
- No duplicate targets within (somamer, direction, cohort)
- Output CSV has all required columns, no NA in target/direction
- Warn if >10% of targets fail annotation join

### Acceptance criteria
- [ ] All 16 LM files load without error
- [ ] Consensus lists contain only Tier 1+2 (both discovery + >=1 validation)
- [ ] Every consensus protein has uniprot and entrez_gene_symbol (or flagged missing)
- [ ] UpSet plots render with replicated intersections highlighted
- [ ] `test_overlap.R` passes
- [ ] Stats report includes totals per cohort, overlap counts per tier, annotation join rate

---

## Step 2: Pathway Enrichment (R, local)

### Objective
ORA enrichment on replicated consensus sets. Focus on axon guidance, semaphorin, neuropilin pathways.

### Input
```
results/01_cross_cohort_overlap/step01_consensus_proteins_pos.csv
results/01_cross_cohort_overlap/step01_consensus_proteins_neg.csv
```

### Method
**ORA (Over-Representation Analysis)** -- input is a discrete significant set, not a ranked list. Use `entrez_gene_id` as primary key.

Databases: GO:BP, GO:MF, GO:CC, KEGG, Reactome.

Background: SomaScan v4.1 panel size (N ~7289). Accept `--background` CSV to override with full unsorted LM results if available.

Parameters: padj < 0.05, qvalue < 0.1, min gene set 10, max gene set 500.

### Scripts (in `analysis/02_pathway_enrichment/`)

**`01_prepare_gene_lists.R`** -- Read consensus CSVs, extract Entrez IDs. Produce 6 lists: pos/neg x VEGFsR1/VEGFsR1.1/combined. Rescue missing Entrez IDs via UniProt-to-Entrez mapping (`AnnotationDbi`).

**`02_run_enrichment.R`** -- Run `clusterProfiler::enrichGO()` (BP, MF, CC), `clusterProfiler::enrichKEGG()`, `ReactomePA::enrichPathway()` per gene list. Save `.rds` and `.csv`.

**`03_plot_enrichment.R`** -- Plots:
- Dot plots: top 20 terms per database per direction
- Combined dot plot faceted by database for positive consensus
- Highlight plot: filter for `axon|semaphorin|neuropilin|plexin|VEGF|angiogenesis` terms, focused bar chart
- Concept network plot (`cnetplot`) for top GO:BP terms

**`run_enrichment.R`** -- Entrypoint. CLI args: `--consensus-dir`, `--results-dir`, `--background`.

### R packages
`clusterProfiler`, `ReactomePA`, `enrichplot`, `org.Hs.eg.db`, `AnnotationDbi`, `tidyverse`, `optparse`

### Output (`results/02_pathway_enrichment/`)
```
step02_enrichment_GO_BP_pos.csv
step02_enrichment_GO_MF_pos.csv
step02_enrichment_GO_CC_pos.csv
step02_enrichment_KEGG_pos.csv
step02_enrichment_Reactome_pos.csv
step02_enrichment_GO_BP_neg.csv   (and MF, CC, KEGG, Reactome)
step02_dotplot_GO_BP_pos.pdf
step02_dotplot_KEGG_pos.pdf
step02_dotplot_Reactome_pos.pdf
step02_dotplot_combined_pos.pdf
step02_highlight_axon_semaphorin.pdf
step02_cnetplot_GO_BP_pos.pdf
step02_enrichment_objects.rds
step02_stats_report.txt
README.md
```

### Test: `test_enrichment.R`
- No NA Entrez IDs passed to enrichGO
- Mapping rate >= 80% of consensus proteins
- Enrichment runs on mock gene list (SEMA3A, NRP1, NRP2, PLXNA1, FLT1)
- Output CSV schema correct
- Highlight regex captures "Semaphorin-plexin signaling"

### Acceptance criteria
- [ ] >=80% Entrez ID mapping rate
- [ ] Enrichment runs for all 5 databases x 2 directions
- [ ] Stats report: N genes tested, N significant terms per database, top 5 per database
- [ ] Axon/semaphorin terms highlighted if significant; explicit null reported if not
- [ ] `test_enrichment.R` passes

---

## Step 3: Structural Prediction (Python, local + HPC)

### Objective
AlphaFold Multimer prediction for sFLT1 (D1-D3, residues 1-338 of FLT1 P17948) vs top candidate proteins.

### Input
```
results/01_cross_cohort_overlap/step01_consensus_proteins_pos.csv
results/02_pathway_enrichment/step02_enrichment_GO_BP_pos.csv
```

### Candidate selection (N = 10-20)
Priority order:
1. Tier 1, dual-somamer, AND in enriched axon/semaphorin/neuropilin pathway
2. Tier 1, dual-somamer
3. Tier 2, in enriched axon/semaphorin pathway
4. Positive controls: VEGFA (P15692, known sFLT1 ligand -- pipeline calibration), NRP1, NRP2, SEMA3A

VEGFA is mandatory positive control. If AlphaFold does not produce ipTM > 0.7 for sFLT1-VEGFA, the method has a calibration problem.

### sFLT1 sequence
FLT1 P17948, Ig-like domains 1-3 (residues 1-338). Soluble form. Document exact sequence used.

### Scoring metrics
- **ipTM > 0.6 AND mean inter-chain PAE < 15A** = predicted interaction
- **ipTM > 0.8 AND PAE < 10A** = high-confidence predicted interaction
- Also report: pTM, mean pLDDT at interface residues (contact < 8A), N interface residues

### Scripts (in `analysis/03_structural_prediction/`)

**`01_select_candidates.py`** -- Rank candidates by priority. Output `step03_candidates.csv`: `rank, target, uniprot, entrez_gene_symbol, tier, dual_somamer, in_axon_pathway, rationale`.

**`02_fetch_sequences.py`** -- Fetch FASTA from UniProt REST API. Write two-chain FASTA files to `fasta/sflt1_vs_{target}.fasta`. Trim sFLT1 to D1-D3. Retry + rate limiting.

**`03_generate_lsf_jobs.py`** -- Generate LSF submission scripts per candidate:
```
#BSUB -J af2_{target}
#BSUB -q gpu
#BSUB -n 4
#BSUB -R "rusage[mem=16000:ngpus_excl_p=1] span[hosts=1]"
#BSUB -W 24:00
```
AlphaFold multimer, 5 predictions per model (25 total). Generate `submit_all.sh`.

**`04_parse_results.py`** -- Parse `ranking_debug.json` for ipTM/pTM. Load PAE matrices, compute mean inter-chain PAE. Identify interface residues (CA < 8A). Output `step03_interaction_scores.csv`: `target, uniprot, iptm_best, iptm_mean, ptm_best, mean_interchain_pae, mean_interface_plddt, interaction_call, n_interface_residues, rank_by_iptm`.

**`05_plot_results.py`** -- PAE heatmaps per candidate (inter-chain blocks highlighted), ipTM bar chart (VEGFA positive control marked), summary table.

**`run_structural.py`** -- Entrypoint. CLI: `--step` (select/fetch/generate/parse/plot/all), `--consensus-dir`, `--results-dir`, `--n-candidates`.

### Python packages
`pandas`, `numpy`, `requests`, `biopython`, `matplotlib`, `seaborn`, `argparse`, `logging`

### Output (`results/03_structural_prediction/`)
```
step03_candidates.csv
step03_interaction_scores.csv
step03_pae_heatmap_{target}.pdf
step03_iptm_barplot.pdf
step03_summary_table.csv
step03_stats_report.txt
README.md
```

### Test: `test_structural.py`
- Candidate list includes VEGFA positive control
- Two-chain FASTA format correct (mock with cached local FASTA)
- LSF scripts have required #BSUB directives
- ipTM extraction from mock ranking_debug.json
- PAE inter-chain computation on mock array
- Interaction call logic: ipTM=0.85/PAE=8 -> high-confidence; ipTM=0.4/PAE=20 -> low

### Acceptance criteria
- [ ] VEGFA positive control included
- [ ] All FASTA files have exactly 2 chains, sFLT1 = D1-D3
- [ ] LSF scripts syntactically valid
- [ ] Positive control achieves ipTM > 0.7
- [ ] Parsing handles failed runs gracefully
- [ ] `pytest` passes on mock data

---

## Execution Sequence

```
Step 1 (R, local)
  |
  v
Step 2 (R, local) -- needs Step 1 consensus CSVs
  |
  v
Step 3a (Python, local) -- candidate selection + seq fetch + job generation
  |                        needs Step 1 consensus + Step 2 enrichment
  v
Step 3b (HPC, GPU) -- AlphaFold runs (~24h per candidate, 10-20 GPU-days)
  |
  v
Step 3c (Python, local) -- result parsing + plotting
```

## Known Risks

1. **Sparse negative overlap:** validation cohorts have 3-9 neg proteins. Tier 2 neg may be empty. Report as finding.
2. **Annotation join misses:** SomaScan Target names may not match annotations (`.N` suffixes). Fuzzy matching needed.
3. **Background set approximate:** only significant LM results available, not full panel. Use SomaScan v4.1 size (~7289) as conservative universe.
4. **AlphaFold compute:** 12-24h/candidate on GPU. Budget 10-20 GPU-days for 10-20 candidates.
5. **sFLT1 isoform:** document whether using canonical D1-D3 (338 residues) or including sFLT1-specific C-terminal tail (31 residues from intron 13 read-through).
