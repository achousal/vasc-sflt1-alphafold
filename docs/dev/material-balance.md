# Material Balance: vasc-sflt1-alphafold

Two stages: (1) linear regression analysis on SomaScan proteomics, (2) AlphaFold multimer structural prediction on Minerva HPC. Stage 1 ran upstream (Elahi Lab R script). Stage 2 runs locally + HPC.

---

## Stage 1: Linear Regression Analysis (upstream)

Ran by Elahi Lab. Source script: `data/markVCID_angiogenesis_synapticPlasticity.R`. Outputs are frozen CSVs in `data/`.

### Input data

SomaScan 11k proteomics (CSF). Four cohorts, two of which include age adjustment:

| Cohort | Role | Age-adjusted | Platform |
|--------|------|-------------|----------|
| MarkVCID | Discovery | No | SomaScan |
| UCSF_AD | Discovery | Yes | SomaScan |
| GNPC | Validation | No | SomaScan |
| WASHU | Validation | Yes | SomaScan |

Raw inputs to the R script (not in this repo):
- `MarkVCID_CSF_Proteomics_data.csv` -- SomaScan RFU matrix, samples x analytes
- `MarkVCID_CSF_Protein_metadata.csv` -- analyte annotation (Target, UniProt, EntrezGeneSymbol)
- `MarkVCID with CSF samples - Somalogic.xlsx` -- clinical data
- `MarkVCID1_anonID_subject_map.xlsx` -- sample-to-subject mapping

### Transformation

**Preprocessing:**
1. Filter to `SampleType == "Sample"` (drop calibrators/buffers)
2. Merge proteomics with clinical data via ExtIdentifier
3. Log10-transform all protein columns
4. Z-score scale each column (`scale()`)

**Linear model (per protein, per somamer):**

```
protein_i ~ VEGFsR1 + Age
```

and separately:

```
protein_i ~ VEGFsR1.1 + Age
```

- `VEGFsR1` = SomaScan somamer 1 for sFLT1 (VEGF sR1)
- `VEGFsR1.1` = SomaScan somamer 2 for sFLT1 (VEGF sR1.1)
- `beta` = coefficient on the sFLT1 somamer term (row 2, col 1 of `summary(lm)$coefficients`)
- `pvalue` = p-value for that coefficient (row 2, col 4)
- `padj` = Bonferroni correction across all ~7,335 proteins
- Significance gate: `padj <= 0.05` AND `|beta| >= 0.5`
- Direction: sign of beta

For UCSF_AD and WASHU, the same model was run but with age adjustment applied upstream in the data (details unconfirmed -- the `_age` suffix on their filenames is the only indicator).

### Output files (frozen in `data/`)

Per cohort, 4 LM result CSVs + annotation CSVs:

**MarkVCID** (`data/MarkVCID/`):

| File | Rows | Content |
|------|------|---------|
| `values_LM_VEGFsR1_sig_pos.csv` | 1,117 | Significant positive associations, somamer 1 |
| `values_LM_VEGFsR1_sig_neg.csv` | 318 | Significant negative associations, somamer 1 |
| `values_LM_VEGFsR1.1_sig_pos.csv` | 1,205 | Significant positive associations, somamer 2 |
| `values_LM_VEGFsR1.1_sig_neg.csv` | 360 | Significant negative associations, somamer 2 |
| `discovery_vasc_markvcid_pos_proteins_common.csv` | -- | Annotation: name, Analytes, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol |
| `discovery_vasc_markvcid_neg_proteins_common.csv` | -- | Same schema, negative direction |

**UCSF_AD** (`data/UCSF_AD/`): same 4 LM CSVs (filenames have `_age` suffix), plus `discovery_ucsfad_{pos,neg}_proteins_common_adj_age.csv`.

**GNPC** (`data/GNPC/`): same 4 LM CSVs, plus `validation_gnpc_pos_proteins_common.csv`.

**WASHU** (`data/WASHU/`): same 4 LM CSVs, plus `validation_washu_{pos,neg}_proteins_common_between_the_2_somamers_adj_for_age.csv`.

**LM CSV schema:**
```
"","Target","beta","pvalue","padj","minus_log_padj","sig_prot","Significant_Direction"
```

**Annotation CSV schema (MarkVCID):**
```
"","name","Analytes","TargetFullName","Target","UniProt","EntrezGeneID","EntrezGeneSymbol","Organism","ID"
```

**Annotation CSV schema (UCSF_AD):**
```
"name","Analytes","SeqId","SomaId","TargetFullName","Target","UniProt","EntrezGeneID","EntrezGeneSymbol","Organism","ID"
```

Note: annotation schemas differ across cohorts. Column names and column counts are not uniform.

### Audit points for Stage 1

- [ ] UCSF_AD and WASHU have `_age` filenames but the LM equation in the MarkVCID script already includes Age. Confirm whether "age-adjusted" means a different model or different preprocessing.
- [ ] The MarkVCID script self-excludes `VEGF sR1` and `VEGF sR1.1` from results (line 626). Confirm the other cohorts do the same.
- [ ] Bonferroni correction is conservative. Confirm all 4 cohorts use the same multiple testing method.
- [ ] The `|beta| >= 0.5` threshold is on Z-scored log10 data. This is a standardized effect size filter, not a fold-change filter.
- [ ] Annotation schemas differ across cohorts. The overlap code must handle column name mismatches.

---

## Stage 1 → Stage 2 bridge: Cross-Cohort Overlap + Pathway Enrichment

Ran locally. Code in `analysis/01_cross_cohort_overlap/` (R) and `analysis/02_pathway_enrichment/` (R).

### Step 1: Cross-Cohort Overlap

**Input:** 16 LM result CSVs (4 cohorts x 2 somamers x 2 directions) + annotation CSVs.

**Process:**
1. Extract significant protein sets from each CSV
2. Compute set intersections across 4 cohorts per somamer-direction stratum
3. Assign tier: Tier 1 = all 4 cohorts, Tier 2 = 3 cohorts, Tier 3 = 2, Tier 4 = 1
4. Flag dual-somamer proteins (significant under both VEGFsR1 and VEGFsR1.1)
5. Join annotations (99.1% join rate for Tier 1+2)

**Counts from stats report:**

| Somamer | Direction | Tier 1 | Tier 2 | Tier 3 |
|---------|-----------|--------|--------|--------|
| VEGFsR1 | pos | 328 | 51 | 36 |
| VEGFsR1 | neg | 0 | 3 | 24 |
| VEGFsR1.1 | pos | 218 | 447 | 168 |
| VEGFsR1.1 | neg | 0 | 1 | 33 |

Dual-somamer Tier 1+2: **678 proteins**

**Output files** (`results/01_cross_cohort_overlap/`):

| File | Content |
|------|---------|
| `step01_consensus_proteins_pos.csv` | Tier 1+2 positive proteins. **Feeds Step 2 and Step 3a.** |
| `step01_consensus_proteins_neg.csv` | Tier 1+2 negative proteins. Feeds Step 2 only. |
| `step01_all_overlaps.csv` | All proteins, all tiers. Audit reference. |
| `step01_beta_heatmap.pdf` | Beta consistency across cohorts |
| `step01_upset_VEGFsR1_pos.pdf` | UpSet plot, somamer 1, positive |
| `step01_upset_VEGFsR1_neg.pdf` | UpSet plot, somamer 1, negative |
| `step01_upset_VEGFsR1_1_pos.pdf` | UpSet plot, somamer 2, positive |
| `step01_upset_VEGFsR1_1_neg.pdf` | UpSet plot, somamer 2, negative |
| `step01_tier_summary.pdf` | Tier distribution |
| `step01_stats_report.txt` | Full counts and protein list |

**Consensus CSV schema:**
```
target,uniprot,entrez_gene_symbol,entrez_gene_id,direction,somamer,tier,n_cohorts,cohorts_sig,mean_beta,min_padj,dual_somamer
```

### Step 2: Pathway Enrichment

**Input:** `step01_consensus_proteins_pos.csv`, `step01_consensus_proteins_neg.csv`.

**Process:**
1. Map proteins to Entrez gene IDs (97.9% mapping rate for positive, 100% for negative)
2. Run Over-Representation Analysis (ORA) via clusterProfiler against: GO BP, GO CC, GO MF, KEGG, Reactome
3. BH-FDR correction at alpha = 0.05
4. Run per gene list variant: `pos_combined` (625 genes), `pos_VEGFsR1` (348), `pos_VEGFsR1.1` (587), `neg_combined` (4), `neg_VEGFsR1` (3), `neg_VEGFsR1.1` (1)

**Top pathway hits (pos_combined):**
- Axonogenesis: padj = 8.21e-52 (99/581 genes)
- Synapse assembly: padj = 2.71e-43
- Axon guidance (KEGG): padj = 6.47e-32 (54/355 genes)

49 enriched terms matched axon/semaphorin regex filter.

**Output files** (`results/02_pathway_enrichment/`):

| File pattern | Count | Content |
|-------------|-------|---------|
| `step02_enrichment_{DB}_{DIR}_{LIST}.csv` | 30 | Per-database, per-direction, per-list enrichment results |
| `step02_enrichment_objects.rds` | 1 | R enrichResult objects (21 MB) |
| `step02_dotplot_*.pdf` | 12 | Dotplots per database |
| `step02_cnetplot_*.pdf` | 1 | Gene-term network (GO BP positive) |
| `step02_highlight_axon_semaphorin.pdf` | 1 | Axon/semaphorin focused view |
| `step02_stats_report.txt` | 1 | Counts, top terms, axon regex matches |
| `vegf_depletion/` | dir | Subset analysis with VEGF genes removed |

### Step 3a: Candidate Selection

**Input:** `step01_consensus_proteins_pos.csv` + axon-pathway gene list from Step 2.

**Process:**
1. Intersect consensus proteins with axon-pathway gene list
2. Rank: Tier 1 + dual-somamer + in axon pathway = highest priority
3. Add 4 positive controls: VEGFA, NRP1, NRP2, SEMA3A (Tier 0)
4. Select top 20 data-driven + 4 controls = **24 targets**

**Result:** All 20 data-driven candidates are Tier 1, dual-somamer, in axon pathway. No Tier 2 candidates were needed.

**Output:** `results/03_structural_prediction/step03_candidates.csv`

```
rank,target,uniprot,entrez_gene_symbol,tier,dual_somamer,in_axon_pathway,rationale
1,VEGFA,P15692,VEGFA,0,False,True,Known sFLT1 ligand -- pipeline calibration positive control
2,NRP1,O14786,NRP1,0,False,True,"Neuropilin-1, known VEGF/semaphorin co-receptor"
3,NRP2,O60462,NRP2,0,False,True,"Neuropilin-2, known VEGF/semaphorin co-receptor"
4,SEMA3A,Q14563,SEMA3A,0,False,True,"Class 3 semaphorin, canonical axon guidance ligand"
5-24: BAI1, PLXA1, PLXA4, STMN3, Contactin-5, MER, Dtk, NGL1, NOE1, FSTL4,
      Amyloid-like protein 1, BASI, DSCAM, PCDH9, NRX1B.1, SLIT2, SLIK4,
      APLP2.1, Calcineurin B a, ROBO2
```

---

## Stage 2: AlphaFold Multimer Structural Prediction (HPC)

Code in `analysis/03_structural_prediction/` (Python). Runs on Minerva HPC.

### Step 3b: Sequence Fetch

**Input:** 24 targets from `step03_candidates.csv` with UniProt accessions.

**Process:**
1. Query UniProt REST API for each target: full sequence + topology annotations
2. Apply protein localization filter:
   - Soluble/secreted: full-length
   - Transmembrane: extracellular domain only
   - GPI-anchored: ectodomain only
3. Pair each target with sFLT1 chain in two-chain FASTA format

**sFLT1 construct definitions:**

| Construct | Residues | Length | Description |
|-----------|----------|--------|-------------|
| D1-D3 | 27-330 | 304 aa | Ig-like domains 1-3 (VEGF binding region). Signal peptide (1-26) removed. |
| D1-D6 | 27-657 | 631 aa | Ig-like domains 1-6 |
| D1-D7 | 27-747 | 721 aa | Full extracellular domain |

**Output:** `results/03_structural_prediction/fasta/` -- 25 FASTA files (24 pairs + 1 reference)

FASTA header format:
```
>sFLT1_D1-D3|P17948|residues_27-330
[sequence]
>TARGET|UNIPROT_ID
[sequence]
```

### FASTA inventory (d1d3)

Cross-reference: candidates CSV UniProt → FASTA header UniProt → sequence lengths → batch_summary.json totals. All 24 match.

| Target | UniProt | sFLT1 (aa) | Target (aa) | Total (aa) | Domain | Walltime |
|--------|---------|-----------|-------------|-----------|--------|----------|
| VEGFA | P15692 | 304 | 395 | 699 | full-length | 48h |
| NRP1 | O14786 | 304 | 835 | 1,139 | ectodomain 22-856 | 72h |
| NRP2 | O60462 | 304 | 844 | 1,148 | ectodomain 21-864 | 72h |
| SEMA3A | Q14563 | 304 | 771 | 1,075 | full-length | 72h |
| BAI1 | O14514 | 304 | 918 | 1,222 | ectodomain 31-948 | 96h |
| PLXA1 | Q9UIW2 | 304 | 1,218 | 1,522 | ectodomain 27-1244 | 96h |
| PLXA4 | Q9HCM2 | 304 | 1,214 | 1,518 | ectodomain 24-1237 | 96h |
| STMN3 | Q9NZ72 | 304 | 180 | 484 | full-length | 48h |
| Contactin-5 | O94779 | 304 | 1,100 | 1,404 | full-length | 96h |
| MER | Q12866 | 304 | 485 | 789 | ectodomain 21-505 | 48h |
| Dtk | Q06418 | 304 | 389 | 693 | ectodomain 41-429 | 48h |
| NGL1 | Q9HCJ2 | 304 | 483 | 787 | ectodomain 45-527 | 48h |
| NOE1 | Q99784 | 304 | 485 | 789 | full-length | 48h |
| FSTL4 | Q6MZW2 | 304 | 842 | 1,146 | full-length | 72h |
| Amyloid-like protein 1 | P51693 | 304 | 542 | 846 | ectodomain 39-580 | 72h |
| BASI | P35613 | 304 | 302 | 606 | ectodomain 22-323 | 48h |
| DSCAM | O60469 | 304 | 1,578 | 1,882 | ectodomain 18-1595 | 144h |
| PCDH9 | Q9HC56 | 304 | 791 | 1,095 | ectodomain 24-814 | 72h |
| NRX1B.1 | P58400 | 304 | 346 | 650 | ectodomain 51-396 | 48h |
| SLIT2 | O94813 | 304 | 1,529 | 1,833 | full-length | 144h |
| SLIK4 | Q8IW52 | 304 | 600 | 904 | ectodomain 19-618 | 72h |
| APLP2.1 | Q06481 | 304 | 661 | 965 | ectodomain 32-692 | 72h |
| Calcineurin B a | P63098 | 304 | 170 | 474 | full-length | 48h |
| ROBO2 | Q9HCK4 | 304 | 838 | 1,142 | ectodomain 22-859 | 72h |

**Verification status:**
- All 24 FASTA files have sFLT1 chain = 304 aa (correct, signal peptide removed)
- All total residue counts match `batch_summary.json` exactly
- 14/24 targets use ectodomain extraction (transmembrane proteins)
- 10/24 targets use full-length sequence (soluble/secreted)
- Uncorrected `fasta/` root directory has sFLT1 = 338 aa (wrong, includes signal peptide)

### Step 3c: LSF Job Generation

**Input:** FASTA files + HPC parameters.

**Process:** Generate per-target LSF batch scripts for AlphaFold 2.3.2.

**AF2 command (decoded from LSF base64):**
```bash
cd "/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold/results/03_structural_prediction/d1d3"
module purge
module load alphafold/2.3.2-singularity

export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

singularity run --nv \
    --bind "${AF2DATA}":/data \
    --bind /sc/arion:/sc/arion \
    "$AF2IMAGE" \
    --fasta_paths="$FASTA_ABS" \
    --output_dir="$OUTPUT_ABS" \
    --data_dir=/data \
    --uniref90_database_path=/data/uniref90/uniref90.fasta \
    --mgnify_database_path=/data/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=/data/pdb_mmcif/mmcif_files \
    --obsolete_pdbs_path=/data/pdb_mmcif/obsolete.dat \
    --pdb_seqres_database_path=/data/pdb_seqres/pdb_seqres.txt \
    --uniprot_database_path=/data/uniprot/uniprot.fasta \
    --small_bfd_database_path=/data/small_bfd/bfd-first_non_consensus_sequences.fasta \
    --model_preset=multimer \
    --db_preset=reduced_dbs \
    --max_template_date=2024-01-01 \
    --num_multimer_predictions_per_model=5 \
    --use_precomputed_msas
```

**HPC reference databases (bind-mounted via `$AF2DATA` → `/data`):**

| Container path | Database | Description |
|---------------|----------|-------------|
| `/data/uniref90/uniref90.fasta` | UniRef90 | MSA search (JackHMMER) |
| `/data/mgnify/mgy_clusters_2022_05.fa` | MGnify 2022-05 | MSA search (JackHMMER) |
| `/data/pdb_mmcif/mmcif_files` | PDB mmCIF | Template structures |
| `/data/pdb_mmcif/obsolete.dat` | PDB obsolete list | Filters retired entries |
| `/data/pdb_seqres/pdb_seqres.txt` | PDB seqres | Template sequence search (HMMSearch) |
| `/data/uniprot/uniprot.fasta` | UniProt | Pairing for multimer MSAs |
| `/data/small_bfd/bfd-first_non_consensus_sequences.fasta` | Small BFD | Reduced DB replacement for full BFD (avoids HHblits titin crash) |

`$AF2DATA` and `$AF2IMAGE` are set by `module load alphafold/2.3.2-singularity` on Minerva. The `--bind /sc/arion:/sc/arion` mount gives the container access to project FASTA and output directories.

`--db_preset=reduced_dbs` uses Small BFD + JackHMMER instead of full BFD + HHblits. This avoids the HHblits titin-length sequence crash (see ADR 001).

**HPC resources per job:**

| Parameter | Value |
|-----------|-------|
| Queue | gpu |
| Cores | 4 |
| Memory | 16 GB/core |
| GPU | 1x A100 (40 GB) |
| Account | acc_vascbrain |

**Three batches generated:**

| Batch | sFLT1 construct | Targets | Job scripts | Walltime range | HPC root |
|-------|----------------|---------|-------------|----------------|----------|
| `d1d3` | D1-D3 (27-330) | 24 | 24 | 48-144h | `/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold/results/03_structural_prediction/d1d3` |
| `d1d6` | D1-D6 (27-657) | 11 | 11 | 24h (all) | same pattern |
| `d1d7` | D1-D7 (27-747) | 9 | 9 | 24h (all) | same pattern |

**d1d6 FASTA inventory (11 targets, sFLT1 = 631 aa):**

| Target | sFLT1 (aa) | Target (aa) | Total (aa) |
|--------|-----------|-------------|-----------|
| VEGFA | 631 | 395 | 1,026 |
| NRP1 | 631 | 835 | 1,466 |
| SEMA3A | 631 | 771 | 1,402 |
| Contactin-5 | 631 | 1,100 | 1,731 |
| NGL1 | 631 | 483 | 1,114 |
| NOE1 | 631 | 485 | 1,116 |
| Amyloid-like protein 1 | 631 | 542 | 1,173 |
| BASI | 631 | 302 | 933 |
| PCDH9 | 631 | 791 | 1,422 |
| SLIT2 | 631 | 1,529 | 2,160 |
| SLIK4 | 631 | 600 | 1,231 |

**d1d7 FASTA inventory (9 targets, sFLT1 = 721 aa):**

| Target | sFLT1 (aa) | Target (aa) | Total (aa) |
|--------|-----------|-------------|-----------|
| VEGFA | 721 | 395 | 1,116 |
| NRP1 | 721 | 835 | 1,556 |
| SEMA3A | 721 | 771 | 1,492 |
| Contactin-5 | 721 | 1,100 | 1,821 |
| NGL1 | 721 | 483 | 1,204 |
| NOE1 | 721 | 485 | 1,206 |
| Amyloid-like protein 1 | 721 | 542 | 1,263 |
| BASI | 721 | 302 | 1,023 |
| SLIK4 | 721 | 600 | 1,321 |

Note: d1d6 and d1d7 are subsets of the d1d3 targets. d1d6 drops 13 targets (all large ectodomain proteins). d1d7 drops PCDH9 and SLIT2 further. Target sequences are identical across batches -- only the sFLT1 chain changes. All d1d6/d1d7 walltimes are set to 24h in the manifest (may need adjustment for SLIT2 at 2,160 total residues).

**Walltime scaling (d1d3, total residues = sFLT1 304aa + target):**

| Residues | Walltime | Example targets |
|----------|----------|----------------|
| <700 | 48h | VEGFA (699), STMN3 (484), Calcineurin B a (474) |
| 700-1200 | 72h | NRP1 (1139), SEMA3A (1075), PCDH9 (1095) |
| 1200-1600 | 96h | BAI1 (1222), PLXA1 (1522), Contactin-5 (1404) |
| >1600 | 144h | DSCAM (1882), SLIT2 (1833) |

**Output files per batch** (`results/03_structural_prediction/{batch}/jobs/`):

| File | Content |
|------|---------|
| `af2_{TARGET}.lsf` | Per-target LSF script |
| `manifest.json` | Job metadata for all targets |
| `wrapper.sh` | Sentinel-based exit trap (writes completed.log or failed.log) |
| `submit_all.sh` | Sequential submission script |
| `submit_orchestrator.sh` | Parallel submission with monitoring |
| `batch_summary.json` | Construct info, walltime per target |

### HPC directory layout (Minerva)

```
/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold/
└── results/03_structural_prediction/
    ├── d1d3/
    │   ├── fasta/                    # 24 corrected two-chain FASTA files
    │   │   └── sflt1_vs_{TARGET}.fasta
    │   ├── jobs/                     # 24 LSF scripts + orchestration
    │   │   ├── af2_{TARGET}.lsf
    │   │   ├── manifest.json
    │   │   ├── wrapper.sh
    │   │   ├── submit_all.sh
    │   │   └── submit_orchestrator.sh
    │   └── results/                  # AF2 output (one subdir per target)
    │       └── {TARGET}/
    │           └── sflt1_vs_{TARGET}/
    │               ├── ranking_debug.json
    │               ├── result_model_*.pkl
    │               ├── relaxed_model_*.pdb  (if relax=true)
    │               ├── unrelaxed_model_*.pdb
    │               └── msas/
    ├── d1d6/                         # same structure, 11 targets
    └── d1d7/                         # same structure, 9 targets
```

### Step 3d: Parse AF2 Results (pending)

**Input:** AF2 output directories on Minerva.

**Process:**
1. Read `ranking_debug.json` per target: extract best model, ipTM, pTM
2. Load PAE matrix from best model
3. Compute inter-chain PAE: mean of off-diagonal blocks (chain A→B and B→A)
4. Classify:
   - `high_confidence`: ipTM > 0.8 AND mean inter-chain PAE < 10 A
   - `predicted`: ipTM > 0.6 AND mean inter-chain PAE < 15 A
   - `low_confidence`: below thresholds

**Output:** `results/03_structural_prediction/step03_interaction_scores.csv`

```
target,uniprot,entrez_gene_symbol,iptm_best,ptm_best,mean_inter_chain_pae,
n_interface_residues,interface_plddt,interaction_call,n_models_parsed
```

### Step 3e: Plot Results (pending)

**Output:**
- `step03_iptm_barplot_{construct}.pdf` -- ranked bar chart
- `step03_pae_heatmap_{TARGET}.pdf` -- per-target PAE matrix
- `step03_summary_table.csv` -- ranked summary
- `step03_stats_report.txt` -- positive control validation

---

## Files to audit

### Process-critical (errors invalidate downstream)

| File | What to verify |
|------|---------------|
| `data/MarkVCID/values_LM_VEGFsR1_sig_pos.csv` | Correct LM equation, correct significance gate, self-exclusion of sFLT1 somamers |
| `data/UCSF_AD/values_LM_VEGFsR1_sig_pos_age.csv` | Same model as MarkVCID? Or different age adjustment? |
| `data/GNPC/`, `data/WASHU/` | Same model? Same significance thresholds? |
| Annotation CSVs (all 4 cohorts) | Column name differences between cohorts -- does the overlap code handle them? |
| `step01_consensus_proteins_pos.csv` | Tier assignment logic correct; annotation join rate (99.1% -- what about the 0.9%?) |
| `d1d3/fasta/sflt1_vs_*.fasta` | All headers say `residues_27-330` (not `residues_1-338`) |
| `step03_candidates.csv` | All 24 targets present with correct UniProt accessions |

### Interpretation-critical (errors bias results)

| File | What to verify |
|------|---------------|
| `step02_enrichment_*.csv` | Gene universe definition. What background set was used for ORA? |
| AF2 scoring thresholds | ipTM/PAE cutoffs should be calibrated against VEGFA positive control, not hardcoded |
| `--db_preset=reduced_dbs` | Accepted trade-off for HHblits BFD bug. Quality margin documented? |

### Reproducibility (errors block replication)

| File | What to verify |
|------|---------------|
| UniProt sequences | Fetched at what date? Pinned or live API? |
| `step02_enrichment_objects.rds` | GO/KEGG/Reactome database versions embedded in R session |
| `--max_template_date=2024-01-01` | Appropriate PDB template cutoff? |
| AF2 container version | `alphafold/2.3.2-singularity` on Minerva |

### Target signal peptide audit

The sFLT1 signal peptide was corrected (338→304 aa) in `d1d3/`. But 5 of the 10 "full-length" targets still include their signal peptides. The 14 ectodomain-extracted targets are correct -- their start residues match UniProt signal peptide boundaries (e.g., APLP1 ectodomain starts at 39, signal ends at 38).

| Target | UniProt | Signal peptide | In FASTA? | Extra residues | Notes |
|--------|---------|---------------|-----------|---------------|-------|
| VEGFA | P15692 | None | N/A | 0 | OK -- no cleavable signal peptide |
| STMN3 | Q9NZ72 | None | N/A | 0 | OK -- cytoplasmic protein |
| Calcineurin B a | P63098 | None | N/A | 0 | OK -- cytoplasmic protein |
| SEMA3A | Q14563 | 1-20 | **Yes** | +20 | Secreted protein, signal peptide included |
| FSTL4 | Q6MZW2 | 1-22 | **Yes** | +22 | Secreted protein, signal peptide included |
| NOE1 | Q99784 | 1-16 | **Yes** | +16 | Secreted protein, signal peptide included |
| SLIT2 | O94813 | 1-30 | **Yes** | +30 | Secreted protein, signal peptide included |
| Contactin-5 | O94779 | 1-18 | **Yes** | +46 | Signal peptide included + GPI-anchor tail (chain ends at 1072, FASTA has 1100) |

**Impact:** Signal peptides are cleaved during maturation and never present on the functional protein. AF2 will likely predict them as disordered, adding noise to the PAE matrix (inflating mean inter-chain PAE) without strongly affecting interface predictions. Contactin-5 has the largest effect (+46 residues, 4.2% of its sequence).

**Affected batches:** Same target sequences are used across d1d3, d1d6, and d1d7 -- the issue propagates to all three constructs.

**Corrective action:** Regenerate FASTA for SEMA3A (21-771), FSTL4 (23-842), NOE1 (17-485), SLIT2 (31-1529), and Contactin-5 (19-1072). Requires re-running AF2 for these 5 targets per construct.

### Known hazards

| Hazard | Status | Location |
|--------|--------|----------|
| Age adjustment asymmetry | Unresolved | MarkVCID/GNPC: no age adj; UCSF_AD/WASHU: age adj |
| HHblits BFD titin crash | Fixed | `--db_preset=reduced_dbs` |
| NRP1/NRP2 binary screen blind spot | Documented | VEGF-bridged co-receptor interactions invisible to AF2 binary multimer |
