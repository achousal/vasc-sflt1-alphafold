# AlphaFold 2.3.2 Multimer -- Process Material Balance

## Purpose

Document every data object that enters, transforms within, and exits the AlphaFold multimer prediction pipeline. This enables auditing of information provenance, storage planning, and identification of reducible intermediates.

Our configuration: `--model_preset=multimer`, `--db_preset=reduced_dbs`, `--run_relax=false`, `--num_multimer_predictions_per_model=5`.

---

## Input Databases (Shared, Static)

These databases are maintained by the HPC facility and mounted read-only into the Singularity container at `/data`. They are not consumed or modified by the pipeline.

| Database | Container Path | Size | Format | Source | Role in Pipeline |
|----------|---------------|------|--------|--------|-----------------|
| UniRef90 | `/data/uniref90/uniref90.fasta` | ~67 GB | FASTA | UniProt | Primary MSA search target for jackhmmer. Provides deep evolutionary sequence diversity. Each chain's query sequence is searched against this database to find homologous sequences that inform co-evolutionary contact prediction. |
| MGnify | `/data/mgnify/mgy_clusters.fa` | ~67 GB | FASTA | EBI MGnify | Metagenomic MSA source. Searched by jackhmmer per chain. Adds environmental sequence diversity not captured in curated databases -- metagenomics samples from soil, gut, ocean, etc. Particularly valuable for proteins with limited representation in UniProt. |
| Small BFD | `/data/small_bfd/bfd-first_non_consensus_sequences.fasta` | ~17 GB | FASTA | BFD (reduced) | Compact version of the Big Fantastic Database. Searched by jackhmmer (not HHblits) under `reduced_dbs` preset. Provides additional sequence diversity from clustered UniProt+metagenomic sources without the full BFD's titin-related HHblits crash risk. |
| PDB mmCIF | `/data/pdb_mmcif/mmcif_files/` | ~250 GB | mmCIF | RCSB PDB | Structural template library. After hmmsearch identifies homologous PDB entries, their 3D coordinates are extracted from these files and used as structural templates to guide folding. Templates provide direct 3D priors where experimental structures exist. |
| PDB obsolete | `/data/pdb_mmcif/obsolete.dat` | ~1 KB | TXT | RCSB PDB | Obsolete PDB entry list. Entries in this file are skipped during template extraction to avoid using superseded or retracted structures. |
| PDB seqres | `/data/pdb_seqres/pdb_seqres.txt` | ~250 MB | FASTA | RCSB PDB | PDB sequence database for template search. hmmsearch queries each chain against this to find structurally characterized homologs. Matches are then cross-referenced against mmCIF files for coordinate extraction. |
| UniProt | `/data/uniprot/uniprot.fasta` | ~100 GB | FASTA | UniProt | Multimer-specific paired MSA construction. jackhmmer searches each chain independently against UniProt, then cross-chain species pairing identifies co-occurring orthologs. This paired MSA provides the inter-chain co-evolutionary signal that distinguishes multimer from monomer prediction. |
| Model parameters | Container-internal | ~12 GB | NumPy (.npz) | DeepMind | Five independently trained neural network weight sets (multimer_v3). Each model has identical architecture but different random initialization during training, providing ensemble diversity. Running all 5 models × 5 random seeds = 25 predictions samples the prediction uncertainty. |

**Total shared input: ~512 GB** (mounted read-only, zero per-job cost)

---

## Per-Job Input

| File | Size | Format | Source | Description |
|------|------|--------|--------|-------------|
| `sflt1_vs_{TARGET}.fasta` | ~1 KB | Two-chain FASTA | `02_fetch_sequences.py` | Contains exactly two sequences separated by `>` headers. Chain A is always the sFLT1 construct (D1-D3: residues 27-330, 304 aa; D1-D6: 631 aa; D1-D7: 721 aa). Chain B is the candidate interaction partner fetched from UniProt. The FASTA stem (`sflt1_vs_VEGFA`) becomes the output subdirectory name. |

---

## Processing Stages

### Stage 1: MSA Search (~30-45 min, CPU-bound)

**What happens:** Each chain's amino acid sequence is independently searched against sequence databases to build a multiple sequence alignment (MSA). The MSA captures evolutionary conservation patterns -- positions that co-vary across species reveal residue-residue contacts.

**Why it matters:** MSA depth and diversity are the primary determinant of prediction quality. A shallow MSA (few homologs) produces uncertain predictions. The multimer-specific UniProt search enables cross-chain pairing: if chain A and chain B orthologs co-occur in the same organism, they likely interact, providing direct inter-chain evolutionary signal.

**Tools and transformations:**

```
Per chain (A, B):

  Query FASTA ──> jackhmmer (1 iter, E=1e-4) vs UniRef90
                  Iterative profile HMM search. Builds an HMM from the query,
                  searches the database, adds hits to the profile, repeats.
                  Output: uniref90_hits.sto (Stockholm alignment, 1-50 MB)

  Query FASTA ──> jackhmmer (1 iter, E=1e-4) vs MGnify
                  Same algorithm, metagenomic database.
                  Output: mgnify_hits.sto (1-50 MB)

  Query FASTA ──> jackhmmer (3 iter, E=1e-3) vs Small BFD
                  More iterations + relaxed E-value for the smaller database.
                  Output: small_bfd_hits.sto (5-200 MB)

  Query FASTA ──> jackhmmer (1 iter) vs UniProt
                  For paired MSA construction (multimer only).
                  Output: uniprot_hits.sto (1-50 MB)
```

**Intermediate products:** 8 Stockholm alignment files (4 per chain) containing aligned homologous sequences. These are the pipeline's most expensive intermediate -- ~45 min of CPU time to generate but reusable across reruns via `--use_precomputed_msas`.

### Stage 2: Template Search (~1-5 min, CPU-bound)

**What happens:** Each chain's sequence is searched against PDB sequences to find experimentally solved structures of homologous proteins. Matching structures provide direct 3D coordinate priors.

**Why it matters:** Templates provide "shortcuts" -- if a close homolog has been crystallized, the model can use its backbone as a starting point rather than folding from scratch. For well-studied proteins like sFLT1 (PDB: 5T89) and VEGFA (multiple structures), templates substantially improve accuracy. For novel proteins without PDB homologs, the model relies entirely on MSA-derived co-evolution.

```
Per chain (A, B):

  Query FASTA ──> hmmsearch vs pdb_seqres.txt
                  Profile HMM search against PDB sequences.
                  Hits filtered by max_template_date (2024-01-01).
                  Output: pdb_hits.sto (<1 MB)

  pdb_hits.sto ──> Template extraction
                   Match PDB IDs, load mmCIF coordinates.
                   Skip entries in obsolete.dat.
                   Extract: backbone coords, sequence, resolution.
```

**Intermediate products:** 2 pdb_hits.sto files + extracted template coordinate arrays. Templates are embedded into `features.pkl` in Stage 3.

### Stage 3: Feature Pipeline (~1-2 min, CPU-bound)

**What happens:** All MSA alignments and template structures are merged into a single feature dictionary -- a set of NumPy arrays that serve as the neural network's input tensor.

**Why it matters:** This is the serialization boundary between data preparation (Stages 1-2, CPU) and neural network inference (Stage 4, GPU). The feature dict standardizes heterogeneous inputs (variable-length MSAs, different template counts) into fixed-shape tensors.

```
  12 .sto files ──┐
                  ├──> DataPipeline.process()
  Template coords ┘
                        │
                        ▼
                  features.pkl (50-500 MB)
                  ├── aatype           [N_res]           residue identity (0-20)
                  ├── residue_index    [N_res]           position in sequence
                  ├── msa              [N_msa, N_res]    aligned sequences
                  ├── msa_mask         [N_msa, N_res]    gap mask
                  ├── deletion_matrix  [N_msa, N_res]    insertion counts
                  ├── template_aatype  [N_tmpl, N_res]   template sequences
                  ├── template_all_atom_positions [N_tmpl, N_res, 37, 3]
                  ├── asym_id          [N_total]         chain assignment
                  └── entity_id        [N_total]         unique sequence ID
```

**Key transformation:** Variable-depth MSAs (thousands to millions of sequences) are cropped/sampled to a fixed `N_msa` for memory. The multimer pipeline concatenates both chains' features and marks chain boundaries via `asym_id`.

### Stage 4: Model Inference (~2-4 hrs, GPU-bound)

**What happens:** Five independently trained neural networks each run the feature dict through the Evoformer attention architecture, producing 3D coordinates and confidence metrics. Each model runs 5 times with different random seeds, yielding 25 total predictions.

**Why it matters:** The Evoformer processes MSA rows and residue pairs simultaneously, learning which co-evolutionary patterns correspond to spatial contacts. The Structure Module then converts these learned representations into 3D atom coordinates through iterative refinement (3 recycling rounds).

```
  For each of 25 (model_i, pred_j) combinations:

  features.pkl ──> Evoformer (48 blocks)
                   ├── MSA row attention:    which aligned sequences are informative?
                   ├── MSA column attention:  which positions co-evolve?
                   └── Pair representation:   residue-residue relationship matrix
                          │
                          ▼
                   Structure Module (8 layers)
                   ├── Invariant Point Attention: 3D-aware attention
                   ├── Backbone update:           refine atom positions
                   └── Side-chain packing:        chi-angle prediction
                          │
                          ▼
                   Confidence Heads
                   ├── pLDDT  [N_res]         per-residue confidence (0-100)
                   ├── PAE    [N_res, N_res]   predicted aligned error (Angstroms)
                   ├── pTM    scalar           predicted TM-score of full complex
                   └── ipTM   scalar           inter-chain TM-score (our key metric)
                          │
                          ▼
                   Outputs per prediction:
                   ├── result_model_{i}_multimer_v3_pred_{j}.pkl    (100-300 MB)
                   └── unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb (~400 KB)
```

**Critical confidence metrics:**
- **ipTM** (inter-chain predicted TM-score): Measures predicted structural accuracy of the interface between chains A and B. Range 0-1, higher = more confident. This is our primary screening metric for sFLT1 interaction candidates.
- **PAE** (predicted aligned error): Matrix of predicted positional errors between all residue pairs. Low inter-chain PAE values indicate confident interface predictions. We extract per-domain PAE submatrices (D1, D2, D3) to classify binding mode.
- **pLDDT**: Per-residue confidence. Stored in the B-factor column of PDB files for visualization.

**GPU memory:** ~8-16 GB for D1-D3 (304 + variable aa). D1-D7 (721 aa) may require unified memory (`TF_FORCE_UNIFIED_MEMORY=1`).

### Stage 5: Ranking (~seconds, CPU)

**What happens:** All 25 predictions are ranked by a composite score: `0.8 × ipTM + 0.2 × pTM`. The best-scoring model becomes `ranked_0`.

```
  25 result*.pkl ──> Extract ipTM, pTM per model
                     Sort by (0.8 × ipTM + 0.2 × pTM)
                     ──> ranking_debug.json
```

---

## Output Inventory (Per Job)

```
{OUTPUT_DIR}/{FASTA_STEM}/
├── msas/                                        MSA intermediates
│   ├── chain_id_map.json          1 KB          Chain→sequence mapping
│   ├── A/                                       sFLT1 chain MSAs
│   │   ├── uniref90_hits.sto      1-50 MB       Evolutionary homologs
│   │   ├── mgnify_hits.sto        1-50 MB       Metagenomic homologs
│   │   ├── small_bfd_hits.sto     5-200 MB      BFD cluster homologs
│   │   ├── uniprot_hits.sto       1-50 MB       For cross-chain pairing
│   │   └── pdb_hits.sto           <1 MB         Structural template hits
│   └── B/                                       Partner chain MSAs
│       └── (same 5 files)
│
├── features.pkl                   50-500 MB     Merged input tensor
│
├── result_model_*_pred_*.pkl      ×25           Model outputs (PAE, ipTM, coords)
│   (100-300 MB each)              2.5-7.5 GB    BULK OF STORAGE
│
├── unrelaxed_model_*_pred_*.pdb   ×25           3D structures (pLDDT in B-factor)
│   (~400 KB each)                 ~10 MB
│
├── ranking_debug.json             ~2 KB         Ranked model list + scores
│
└── timings.json                   ~1 KB         Wall-clock per stage
```

---

## Material Balance Summary

### Per Job (D1-D3 construct, 304 + variable aa)

| Category | In | Out | Ratio |
|----------|-----|-----|-------|
| Unique input | ~1 KB (FASTA) | 2.6-8.6 GB | ~10⁶× amplification |
| Shared input | ~512 GB (databases) | -- (read-only) | -- |
| CPU time | -- | ~45 min (MSA) + ~5 min (features/templates) | -- |
| GPU time | -- | ~2-4 hrs (inference, A100) | -- |
| Peak RAM | -- | ~10 GB | -- |
| Peak VRAM | -- | ~8-16 GB | -- |

### Full Campaign (44 jobs: 24 D1-D3 + 11 D1-D6 + 9 D1-D7)

| Metric | Estimate |
|--------|----------|
| Total storage | ~220 GB |
| Total GPU-hours | ~120-176 hrs (A100) |
| Total wall-clock | ~4-5 hrs per job (parallelizable) |
| MSA reuse savings | ~33 hrs saved on resubmits |

### Storage Dominance

```
  result*.pkl (25 per job)    ████████████████████████████████████  87%
  MSAs (.sto files)           █████                                 10%
  features.pkl                ██                                     3%
  PDB files (25 per job)      ▏                                    <1%
  JSON metadata               ▏                                    <1%
```

The 25 `.pkl` files dominate storage. For downstream analysis, only the top-ranked model's `.pkl` is needed (PAE matrix, ipTM, pTM). The remaining 24 could be archived or deleted after score extraction, yielding **~90% storage reduction**.

---

## Post-Processing (Local)

```
  result*.pkl ──┬──> 04_parse_results.py ──> step03_interaction_scores.csv
                │    Extracts: ipTM, pTM, per-domain PAE submatrices
                │    Classifies: binding mode (D1/D2/D3 dominant interface)
                │    Derives: ipSAE, LIS alternative scoring metrics
                │
                ├──> 05_plot_results.py  ──> step03_iptm_barplot.pdf
                │    Ranks all targets by ipTM with control annotations
                │
                └──> 06_render_structures.py (PyMOL headless)
                     ├── {target}_structure.jpg  (~80 KB, cartoon+surface)
                     └── {target}_pae.jpg        (~80 KB, PAE heatmap)
```

### Downstream Decision Thresholds

| ipTM Range | Interpretation | Action |
|------------|---------------|--------|
| > 0.7 | High-confidence interaction | Proceed to Phase 2 (trimer, extended constructs) |
| 0.5 - 0.7 | Moderate confidence | Consider alternative metrics (ipSAE, LIS), trimer context |
| < 0.5 | Low confidence / no interaction | Deprioritize unless biological prior is strong |

---

## Known Failure Modes

| Failure | Root Cause | Detection | Mitigation |
|---------|-----------|-----------|------------|
| HHblits BFD crash | Titin sequences cause A3M column mismatch | `hhalignment.cpp:1244` in stderr | Use `reduced_dbs` (skip full BFD) |
| Singularity rc=0 on Python error | Container exec layer swallows internal exceptions | Output validation guard checks for `ranking_debug.json` | Post-singularity file existence check in wrapper |
| Amber relax exhaustion | Minimization fails after 100 attempts on some complexes | `ValueError: Minimization failed` in stderr | `--run_relax=false` (relax unnecessary for PPI scoring) |
| GPU CC incompatibility | AF2 2.3.2 (CUDA 11.1) fails on CC > 8.0 GPUs | `ptxas does not support CC 8.9` in stderr | LSF constraint: `select[a100 \|\| v100]` |
| Sentinel race condition | `set -e` in wrapper aborts before `_af_rc` assignment | `completed.log` written for failed jobs | `set +e` around eval in wrapper |
