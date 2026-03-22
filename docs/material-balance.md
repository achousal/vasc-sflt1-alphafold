# Project Material Balance

End-to-end data flow for the sFLT1 interaction screening pipeline. Tracks every file this project creates, transforms, and consumes -- from SomaScan linear model results through AlphaFold structural predictions to scored interaction candidates.

---

## Pipeline Overview

```
 data/                       results/01_             results/02_             results/03_
 4 cohorts                   cross_cohort            pathway_                structural_prediction/
 × 2 somamers                _overlap/               enrichment/             d1d3_corrected/  d1d6/  d1d7/
 × 2 directions
 ═══════════                 ═══════════             ═══════════             ═══════════════════════════════

 16 LM result CSVs ──────>   Consensus list  ──────> Pathway flags  ──────>  Candidate CSVs
 + 8 annotation CSVs         (Tier 1-3)              (axon/sema)             │
                                                                             ├──> FASTA files (UniProt fetch)
 ~3800 significant           1044 proteins           pathway labels          │
 protein associations        across tiers             enriched               ├──> LSF job scripts
                                                                             │
                                                                             ├──> AF2 outputs (Minerva GPU)
                                                                             │    msas/ features/ models/
                                                                             │
                                                                             ├──> Interaction scores CSV
                                                                             │
                                                                             └──> JPEG renders

 FUNNEL:   ~3800 assoc ──> 1044 consensus ──> 24 candidates ──> 44 AF2 jobs ──> ranked hits
```

---

## Stage 1: Input Data

### What we start with

Pre-computed SomaScan linear model results. Each cohort ran a linear model regressing every SomaScan protein on sFLT1 plasma levels, producing per-protein effect sizes and significance. We did not run these models -- they are upstream deliverables.

```
data/
├── MarkVCID/                              Discovery cohort, no age adjustment
│   ├── values_LM_VEGFsR1_sig_pos.csv     1117 proteins positively associated with sFLT1
│   ├── values_LM_VEGFsR1_sig_neg.csv      318 proteins negatively associated
│   ├── values_LM_VEGFsR1.1_sig_pos.csv   VEGFsR1.1 somamer (same protein, different aptamer)
│   ├── values_LM_VEGFsR1.1_sig_neg.csv
│   ├── discovery_vasc_markvcid_pos_proteins_common.csv   Protein annotations (pos)
│   └── discovery_vasc_markvcid_neg_proteins_common.csv   Protein annotations (neg)
│
├── UCSF_AD/                               Discovery cohort, age-adjusted
│   └── (same 6-file structure, suffix "_age")
│
├── GNPC/                                  Validation cohort, no age adjustment
│   └── (same 6-file structure)
│
└── WASHU/                                 Validation cohort, age-adjusted
    └── (same 6-file structure)
```

**Per-file schema (LM results):**

| Column | Type | Description |
|--------|------|-------------|
| Target | string | SomaScan protein short name (e.g., "c-Raf") |
| beta | float | Linear model effect size (sFLT1 → protein association) |
| pvalue | float | Raw p-value |
| padj | float | BH-adjusted p-value |
| minus_log_padj | float | -log10(padj) for visualization |
| sig_prot | string | "Significant" flag |
| Significant_Direction | string | "Positive" or "Negative" |

**Per-file schema (protein annotations):**

| Column | Type | Description |
|--------|------|-------------|
| Target | string | SomaScan protein name |
| UniProt | string | UniProt accession (e.g., P15692) |
| EntrezGeneSymbol | string | HUGO gene symbol |
| EntrezGeneID | int | NCBI gene ID |

**Key properties of input data:**
- Two sFLT1 somamers (VEGFsR1, VEGFsR1.1) measure the same protein with different aptamers. Concordance between somamers increases confidence.
- UCSF_AD and WASHU include age as a covariate; MarkVCID and GNPC do not. This creates a built-in sensitivity analysis for age confounding.
- "Significant" means BH-adjusted p < 0.05 within that cohort's linear model.

---

## Stage 2: Cross-Cohort Overlap (Step 1)

### Transformation

Each cohort independently identified hundreds of sFLT1-associated proteins. Step 1 intersects these lists to find proteins that replicate across cohorts and somamers.

```
 16 LM CSVs (sig proteins per cohort × somamer × direction)
      │
      │  01_select_candidates.py / run_overlap.R
      │  - Union all targets across cohort-somamer combinations
      │  - Count: in how many cohorts is each protein significant?
      │  - Tier assignment by replication depth
      │
      ▼
 step01_consensus_proteins_pos.csv        1044 proteins (positive direction)
 step01_consensus_proteins_neg.csv        analogous (negative direction)
 step01_all_overlaps.csv                  1309 proteins (all tiers, both directions)
```

**Tier assignment logic:**

| Tier | Criteria | Count (pos) | Interpretation |
|------|----------|-------------|----------------|
| 1 | Significant in all 4 cohorts for ≥1 somamer | ~50-80 | Highest confidence: replicates everywhere |
| 2 | Significant in 3 cohorts | ~100-150 | Strong: one cohort dissents |
| 3 | Significant in 2 cohorts | ~200-300 | Moderate: discovery + 1 validation |

**Output schema (consensus):**

| Column | Type | Source |
|--------|------|--------|
| target | string | SomaScan Target name |
| uniprot | string | From annotation CSVs |
| entrez_gene_symbol | string | From annotation CSVs |
| direction | string | pos or neg |
| somamer | string | VEGFsR1 or VEGFsR1.1 |
| tier | int | Computed from cross-cohort overlap |
| n_cohorts | int | Number of cohorts where significant |
| cohorts_sig | string | Semicolon-separated cohort names |
| mean_beta | float | Mean effect size across significant cohorts |
| min_padj | float | Best (smallest) adjusted p-value |
| dual_somamer | bool | Significant for both VEGFsR1 and VEGFsR1.1 |

**What gets dropped:** ~2750 associations that did not replicate in ≥2 cohorts (single-cohort hits, likely false positives or cohort-specific effects).

---

## Stage 3: Pathway Enrichment (Step 2)

### Transformation

The consensus protein list is tested for overrepresentation in biological pathways. This does not filter proteins -- it adds pathway annotations used downstream for candidate prioritization.

```
 step01_consensus_proteins_pos.csv (1044 proteins)
      │
      │  02_pathway_enrichment.R (clusterProfiler ORA)
      │  - GO Biological Process, Cellular Component, Molecular Function
      │  - KEGG pathways
      │  - Reactome pathways
      │  - Per somamer and combined
      │
      ▼
 step02_enrichment_GO_BP_pos_combined.csv     Enriched GO-BP terms
 step02_enrichment_KEGG_pos_combined.csv      Enriched KEGG pathways
 step02_enrichment_Reactome_pos_combined.csv  Enriched Reactome pathways
 (... 15 enrichment CSVs total across ontologies × somamers × directions)
 step02_highlight_axon_semaphorin.pdf         Visual flag: axonogenesis pathway
 step02_enrichment_objects.rds                R objects for downstream reuse
```

**Key pathways flagged:**
- Axon guidance / axonogenesis (GO:0007411) -- primary biological hypothesis
- Semaphorin-plexin signaling -- known sFLT1/NRP1 co-receptor context
- VEGF signaling pathway (KEGG hsa04370) -- positive control pathway

**What this stage produces:** Pathway membership flags (`in_axon_pathway`) that inform candidate ranking in Stage 4. No proteins are filtered out.

---

## Stage 4: Candidate Selection (Step 3a)

### Transformation

Selects the top 20 interaction candidates from the Tier 1 consensus list, plus 4 controls, based on ranking criteria.

```
 step01_consensus_proteins_pos.csv (1044 proteins)
 + step02_enrichment_*.csv (pathway annotations)
      │
      │  01_select_candidates.py
      │  - Filter to Tier 1 (all 4 cohorts significant)
      │  - Rank by: dual_somamer > in_axon_pathway > mean_beta
      │  - Add 4 controls: VEGFA (positive), NRP1, NRP2, SEMA3A (known biology)
      │  - Cap at 24 total
      │
      ▼
 step03_candidates.csv                        24 rows
```

**Output schema:**

| Column | Type | Source |
|--------|------|--------|
| rank | int | Priority order (1 = VEGFA positive control) |
| target | string | SomaScan protein name |
| uniprot | string | UniProt accession for sequence fetch |
| entrez_gene_symbol | string | Gene symbol |
| tier | int | 0 = control, 1 = Tier 1 candidate |
| dual_somamer | bool | From cross-cohort overlap |
| in_axon_pathway | bool | From pathway enrichment |
| rationale | string | Why this candidate was selected |

**Funnel:** 1044 consensus → 24 candidates (97.7% filtered). Selection criteria: Tier 1 replication + pathway relevance + effect size.

---

## Stage 5: Sequence Fetch (Step 3b)

### Transformation

Each candidate's UniProt accession is used to fetch its full amino acid sequence. This is combined with the sFLT1 construct sequence into a two-chain FASTA for AlphaFold.

```
 step03_candidates.csv (24 rows with UniProt accessions)
      │
      │  02_fetch_sequences.py
      │  - HTTP GET to UniProt REST API per accession
      │  - Prepend sFLT1 construct as chain A:
      │      D1-D3: residues 27-330 (304 aa, from PDB 5T89)
      │      D1-D6: 631 aa
      │      D1-D7: 721 aa
      │  - Write two-chain FASTA
      │
      ▼
 fasta/
 ├── sflt1_vs_VEGFA.fasta           >sFLT1_D1-D3|P17948|residues_27-330
 ├── sflt1_vs_NRP1.fasta             SKLKDPELSLK... (304 aa)
 ├── sflt1_vs_PCDH9.fasta           >VEGFA|P15692
 └── ... (24 files)                  MNFLLSWVHW... (variable aa)
```

**Construct variants and their batches:**

| Batch | sFLT1 Construct | Chain A Length | Candidates | Rationale |
|-------|----------------|---------------|------------|-----------|
| d1d3_corrected | D1-D3 (27-330) | 304 aa | 24 (all) | Ligand-binding domains only |
| d1d6 | D1-D6 | 631 aa | 11 (subset) | + Ig-like domains 4-6 |
| d1d7 | D1-D7 | 721 aa | 9 (subset) | Full extracellular domain |

**Total FASTA files:** 24 + 11 + 9 = 44 across 3 batches.

**What the FASTA encodes:** The FASTA is the sole unique input to AlphaFold. Everything AF2 does -- MSA search, template matching, structure prediction -- derives from these two sequences. Chain A is constant (sFLT1 construct); chain B varies per candidate.

---

## Stage 6: AlphaFold Structural Prediction (Step 3c-d, HPC)

### Transformation

Each two-chain FASTA is submitted to AlphaFold Multimer on Minerva GPU nodes. AF2 searches sequence databases for evolutionary homologs, finds structural templates, and predicts 25 3D structures of the complex with confidence scores.

```
 fasta/sflt1_vs_{TARGET}.fasta (~1 KB)
      │
      │  AlphaFold 2.3.2 Multimer (Singularity container, A100/V100 GPU)
      │  --db_preset=reduced_dbs --run_relax=false
      │  ~4-5 hrs per job
      │
      │  Internal stages (not directly visible to us):
      │  1. MSA search: find evolutionary homologs per chain
      │     jackhmmer vs UniRef90, MGnify, Small BFD, UniProt
      │     → 10 .sto alignment files (reusable across reruns)
      │  2. Template search: find solved structures of homologs
      │     hmmsearch vs PDB → structural coordinate priors
      │  3. Feature assembly: merge MSAs + templates into input tensor
      │     → features.pkl (50-500 MB)
      │  4. Model inference: 5 models × 5 seeds = 25 predictions
      │     Evoformer attention → 3D coordinates + confidence scores
      │     → 25 result_*.pkl + 25 unrelaxed_*.pdb
      │  5. Ranking: sort by 0.8×ipTM + 0.2×pTM
      │     → ranking_debug.json
      │
      ▼
 results/{TARGET}/sflt1_vs_{TARGET}/
 ├── ranking_debug.json              Ranked model list + ipTM/pTM scores
 ├── result_model_*_pred_*.pkl ×25   PAE matrices, ipTM, pTM, coordinates
 ├── unrelaxed_model_*_pred_*.pdb ×25  Viewable 3D structures
 ├── features.pkl                    Merged input features
 ├── msas/{A,B}/*.sto               Evolutionary alignments (reusable)
 └── timings.json                    Per-stage wall clock
```

**What AF2 produces that we actually use:**

| Output | Size | What We Extract | Used In |
|--------|------|----------------|---------|
| ranking_debug.json | ~2 KB | Model ranking order, ipTM, pTM per model | Parse step: identifies best model |
| result_*.pkl (best model only) | 100-300 MB | PAE matrix [N×N], ipTM (scalar), pTM (scalar), pLDDT [N] | Parse: domain-resolved scoring. Render: PAE heatmaps |
| unrelaxed_*.pdb (best model only) | ~400 KB | 3D atom coordinates, pLDDT in B-factor | Render: structure JPEGs |
| result_*.pkl (remaining 24) | 2.5-7 GB | Ensemble variance (optional) | Archivable after score extraction |
| msas/*.sto | 50-600 MB | Not used directly | Reusable if job needs resubmission |
| features.pkl | 50-500 MB | Not used directly | Intermediate, archivable |

**Storage per job:** ~3-8 GB total, of which ~90% is redundant `.pkl` files from non-top-ranked models.

**Compute per job:**

| Resource | D1-D3 (304aa) | D1-D6 (631aa) | D1-D7 (721aa) |
|----------|--------------|--------------|--------------|
| Wall clock | ~4-5 hrs | ~6-8 hrs | ~8-12 hrs |
| GPU memory | ~8-10 GB | ~12-16 GB | ~16-24 GB |
| CPU (MSA) | ~45 min | ~45 min | ~45 min |
| GPU (inference) | ~2-4 hrs | ~4-6 hrs | ~6-10 hrs |

---

## Stage 7: Score Parsing (Step 3e, local)

### Transformation

Extracts the biologically meaningful scores from AF2's raw outputs and classifies each candidate by binding mode.

```
 results/{TARGET}/sflt1_vs_{TARGET}/ranking_debug.json
 results/{TARGET}/sflt1_vs_{TARGET}/result_{best_model}.pkl
      │
      │  04_parse_results.py
      │  - Read ranking_debug.json → identify best model
      │  - Load best model's .pkl → extract PAE matrix, ipTM, pTM
      │  - Slice PAE by sFLT1 domain boundaries:
      │      D1 (5-103), D2 (105-198), D3 (199-303)
      │  - Classify binding mode: which domain has lowest inter-chain PAE?
      │  - Compute alternative metrics: ipSAE, LIS (PAE-filtered scores)
      │
      ▼
 step03_interaction_scores.csv              24 rows × scoring columns
```

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| target | string | Protein name |
| uniprot | string | UniProt accession |
| iptm | float | Inter-chain predicted TM-score (0-1, primary metric) |
| ptm | float | Predicted TM-score of full complex |
| iptm_ptm | float | 0.8×ipTM + 0.2×pTM (AF2 ranking score) |
| pae_inter_mean | float | Mean inter-chain PAE (Å) |
| pae_D1_mean | float | Mean PAE for D1 interface residues |
| pae_D2_mean | float | Mean PAE for D2 interface residues |
| pae_D3_mean | float | Mean PAE for D3 interface residues |
| binding_domain | string | D1, D2, D3, or multi (lowest PAE domain) |
| tier | int | From candidate selection |
| category | string | control, novel_high, novel_moderate, novel_low |

**Interpretation thresholds:**

| ipTM | Category | Downstream Action |
|------|----------|------------------|
| > 0.7 | High confidence | Proceed to Phase 2 (trimer, extended constructs) |
| 0.5-0.7 | Moderate | Consider alternative metrics, biological priors |
| < 0.5 | Low / no interaction | Deprioritize |

---

## Stage 8: Visualization (Steps 3f-g, local)

### Transformation

Converts numeric scores and 3D coordinates into publication-ready figures and lightweight images.

```
 step03_interaction_scores.csv
      │
      │  05_plot_results.py
      │
      ├──> step03_iptm_barplot_{construct}.pdf     Bar chart: ipTM by target
      ├──> step03_summary_table.csv                 Formatted results table
      └──> step03_stats_report.txt                  Text summary

 result_{best}.pkl + unrelaxed_{best}.pdb
      │
      │  06_render_structures.py (PyMOL headless)
      │
      ├──> {target}_structure.jpg  (~80 KB)         Cartoon + transparent surface
      │    Chain A (sFLT1) = blue                    Chain B (partner) = orange
      │
      └──> {target}_pae.jpg  (~80 KB)               PAE heatmap with domain lines
           Green = low PAE = high confidence          Red dashes = chain/domain bounds
```

---

## Full Funnel Summary

```
 STAGE          FILES                      ROWS/COUNT          SIZE         OPERATION
 ─────────────  ─────────────────────────  ──────────────────  ───────────  ──────────────────
 Input data     16 LM CSVs + 8 annot      ~3800 sig assoc     ~5 MB        Upstream delivery
                                                │
 Step 1         step01_consensus_pos.csv   1044 proteins       ~200 KB      Cross-cohort
 Overlap                                        │                            intersection
                                           ▼ 97.7% filtered
 Step 2         15 enrichment CSVs         pathway labels      ~2 MB        ORA annotation
 Enrichment     (no filtering, adds flags)      │                            (additive)
                                                │
 Step 3a        step03_candidates.csv      24 candidates       ~3 KB        Tier + pathway
 Selection                                      │                            + effect rank
                                                │
 Step 3b        44 FASTA files             44 sequence pairs   ~50 KB       UniProt fetch +
 Fetch          (3 construct batches)           │                            construct join
                                                │
 Step 3c-d      44 AF2 output dirs         25 models × 44      ~220 GB      GPU inference
 AF2            ranking + pkl + pdb             │               (bulk)       on Minerva
                                                │
 Step 3e        interaction_scores.csv     24-44 scored rows   ~5 KB        PAE extraction
 Parse                                          │                            + domain classify
                                                │
 Step 3f-g      48-88 JPEGs               structure + PAE      ~7 MB        PyMOL render
 Render         per target per batch            │
                                                ▼
                                           RANKED HIT LIST
                                           ipTM > 0.7 → Phase 2
```

---

## Storage Budget

| Category | Per Job | 44 Jobs Total | Archivable? |
|----------|---------|---------------|-------------|
| Input (LM CSVs) | -- | ~5 MB | No (source data) |
| Step 1-2 outputs | -- | ~2 MB | No (small) |
| FASTA files | ~1 KB | ~50 KB | No (small) |
| MSAs (.sto) | 50-600 MB | ~10 GB | Yes, after rerun confidence |
| features.pkl | 50-500 MB | ~10 GB | Yes, after parsing |
| result*.pkl (×25) | 2.5-7.5 GB | **~220 GB** | 24/25 per job after parsing |
| unrelaxed*.pdb (×25) | ~10 MB | ~440 MB | 24/25 per job after rendering |
| ranking_debug.json | ~2 KB | ~90 KB | No (key metadata) |
| interaction_scores.csv | ~5 KB | ~15 KB | No (final output) |
| JPEG renders | ~160 KB | ~7 MB | No (deliverable) |

**Reduction opportunity:** After parsing and rendering, retain only the top-ranked model's `.pkl` and `.pdb` per job. Archive or delete the other 24. **~200 GB → ~15 GB** (93% reduction).

---

## Configuration (Current)

| Parameter | Value | Why |
|-----------|-------|-----|
| `--db_preset` | `reduced_dbs` | Full BFD causes HHblits crash on titin sequences |
| `--run_relax` | `false` | Amber relax fails on sFLT1 complexes; unnecessary for PPI scoring |
| `--use_precomputed_msas` | `true` | Reuse MSAs from prior runs on resubmission (~45 min saved) |
| `--num_multimer_predictions_per_model` | `5` | 5 models × 5 seeds = 25 predictions for ensemble |
| `--max_template_date` | `2024-01-01` | Exclude structures deposited after this date |
| GPU constraint | `select[a100 \|\| v100]` | AF2 2.3.2 (CUDA 11.1) incompatible with CC > 8.0 GPUs |
| `TF_FORCE_UNIFIED_MEMORY` | `1` | Allow GPU memory overflow to system RAM for D1-D6/D1-D7 |
