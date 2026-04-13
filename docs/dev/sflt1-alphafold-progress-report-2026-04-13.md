# sFLT1 AlphaFold Multimer Interaction Screen — Progress Report

**Date:** 2026-04-13
**Project:** vasc-sflt1-alphafold
**PI:** Fanny Elahi
**Analyst:** Andres Chousal

---

## Objective

Test whether sFLT1 physically interacts with axon guidance proteins beyond its known VEGF ligand, using AlphaFold Multimer predictions on candidates surfaced by cross-cohort SomaScan proteomics.

---

## Screen Design

Candidates derive from SomaScan cross-cohort overlap (MarkVCID, UCSF_AD, GNPC, WASHU), filtered for dual-somamer replication and axon/semaphorin pathway membership. Each candidate is screened as a binary pair with sFLT1 across three construct lengths:

| Construct | sFLT1 residues | Length | Purpose |
|-----------|---------------|--------|---------|
| D1-D3 | 27-330 | 304 aa | VEGF-binding domains — known interaction surface |
| D1-D6 | 27-654 | 631 aa | Adds dimerization/co-receptor domains |
| D1-D7 | 27-750 | 721 aa | Full ectodomain |

Each target gets 25 models (5 predictions × 5 model seeds). VEGFA is included as a positive control in the fullscreen run.

**Compute:** AlphaFold 2.3.2, multimer preset, reduced_dbs, NVIDIA V100/A100 GPUs (Minerva HPC). Pool orchestrator manages submission with a 60-job cap across 10 GPU nodes.

---

## Current Completion

**Pool orchestrator status (2026-04-13 11:17):**

| Construct | Pool attempted | Scorable targets | Fully scored (≥5 models) | % of scorable |
|-----------|----------------|------------------|--------------------------|---------------|
| D1-D3 | 220/283 (ok 217, fail 3) | 196 | **108** (+ VEGFA ctrl) | **55%** |
| D1-D6 | 73/283 (ok 72, fail 1) | 196 | **33** (+ VEGFA ctrl) | **17%** |
| D1-D7 | 75/283 (ok 74, fail 1) | 196 | **32** (+ VEGFA ctrl) | **16%** |

**Total:** 368/849 pool attempts (43%). HPC pool is active.

**Change since 2026-04-10:** D1-D3 scored 87 → 108 (+21). D1-D6 27 → 33 (+6). D1-D7 25 → 32 (+7). Cross-construct N (all 3 constructs) grew 24 → 33 targets.

---

## Positive Control

**VEGFA (all constructs):** D1-D3 mean 0.790 | D1-D6 mean 0.605 | D1-D7 mean 0.561. The control fires strongly in D1-D3 (the known binding interface) and attenuates modestly in longer constructs, consistent with template coverage bias. **Screen is functioning** — a true binary partner produces a clearly-separated signal (~9 SD above the D1-D3 background median).

---

## Results — D1-D3 (108 targets scored + VEGFA control)

Mean ranking confidence (iptm+ptm averaged across all 25 models) is the reliable metric.

**Top targets by mean confidence (VEGFA shown as positive control):**

| Target | mean | best | Δ | n |
|--------|------|------|---|---|
| **VEGFA (ctrl)** | **0.790** | 0.800 | 0.01 | 25 |
| **TESC** | **0.464** | 0.665 | 0.20 | 25 |
| AP2A2 | 0.445 | 0.840 | 0.39 | 25 |
| UNC5H4.2 | 0.418 | 0.675 | 0.26 | 25 |
| BASI | 0.412 | 0.653 | 0.24 | 25 |
| CLC2L | 0.400 | 0.594 | 0.19 | 25 |
| AGRB2 | 0.385 | 0.632 | 0.25 | 25 |
| PCDH9 | 0.379 | 0.750 | 0.37 | 25 |
| NXPH2 | 0.370 | 0.603 | 0.23 | 25 |
| NLGN1.1 | 0.356 | 0.598 | 0.24 | 25 |
| NLGN1 | 0.351 | 0.615 | 0.26 | 25 |
| CHP1 | 0.351 | 0.457 | 0.11 | 25 |
| STX1a.1 | 0.349 | 0.767 | 0.42 | 25 |

**Mean distribution (108 targets):**

| Range | N targets |
|-------|-----------|
| < 0.25 | 50 (46%) |
| 0.25-0.30 | 34 (31%) |
| 0.30-0.40 | 20 (19%) |
| 0.40-0.50 | 4 (4%) |
| > 0.50 | 0 (0%) |

Median = 0.253 | Mean = 0.267 | SD = 0.058.

**No non-control candidate reaches mean 0.50.** TESC (0.464) remains the top hit, ~3.6 SD above median. AP2A2 is intracellular (adaptin complex subunit) — biologically implausible. **UNC5H4 (UNC5D, netrin repulsion receptor) is a new top-5 entrant** since the last report — biologically plausible as an axon guidance partner and worth tracking as D1-D6/D1-D7 data arrives. The neuroligin cluster (NLGN1/1.1/1.2) remains the most biologically coherent mid-tier signal.

---

## Results — Cross-Construct Comparison (33 targets with all 3 constructs)

Sorted by D1-D7 mean (descending). VEGFA shown for reference.

| Target | D1-D3 | D1-D6 | D1-D7 | Trend |
|--------|-------|-------|-------|-------|
| **VEGFA (ctrl)** | **0.790** | **0.605** | **0.561** | Falling |
| **PLXA4** | 0.240 | 0.270 | **0.338** | **Rising** |
| NLGN2_ECD | 0.315 | 0.323 | 0.319 | Flat |
| **LRP12** | 0.252 | 0.245 | **0.303** | **Rising** |
| **MER** | 0.259 | 0.195 | **0.296** | **Rising** |
| BASI | 0.412 | 0.371 | 0.284 | Falling |
| **KIRR3** | 0.182 | 0.183 | **0.279** | **Rising** |
| SVIP | 0.276 | 0.301 | 0.270 | Flat |
| BMPER | 0.206 | 0.276 | 0.268 | Flat |
| NXPH2 | 0.370 | 0.382 | 0.267 | Falling |
| NOE1 | 0.325 | 0.296 | 0.265 | Falling |
| Amyloid-like_protein_1 | 0.285 | 0.248 | 0.265 | Flat |
| BAI1 | 0.284 | 0.220 | 0.259 | Flat |
| PLXA1 | 0.243 | 0.239 | 0.257 | Flat |
| Contactin-5 | 0.215 | 0.237 | 0.245 | Rising |
| FSTL4 | 0.258 | 0.273 | 0.245 | Flat |
| NLGNX | 0.334 | 0.285 | 0.240 | Falling |
| MMP-16.1 | 0.220 | 0.261 | 0.237 | Flat |
| NLGN1.1 | 0.356 | 0.310 | 0.237 | Falling |
| NGL1 | 0.268 | 0.245 | 0.235 | Falling |
| NLGN1 | 0.351 | 0.317 | 0.235 | Falling |
| NLGN1.2 | 0.347 | 0.337 | 0.235 | Falling |
| LRP11.1 | 0.206 | 0.254 | 0.232 | Flat |
| Dtk | 0.192 | 0.191 | 0.223 | Rising |
| BAI3 | 0.251 | 0.228 | 0.222 | Flat |
| LRFN4 | 0.226 | 0.215 | 0.220 | Flat |
| CRIM1_ECD.1 | 0.273 | 0.244 | 0.220 | Falling |
| GPC5 | 0.235 | 0.213 | 0.219 | Flat |
| TWEAK | 0.244 | 0.282 | 0.214 | Flat |
| HPLN4 | 0.226 | 0.243 | 0.214 | Flat |
| sICAM-5 | 0.274 | 0.206 | 0.199 | Falling |
| sICAM-5.1 | 0.273 | 0.203 | 0.196 | Falling |
| Cadherin-12_ECD | 0.197 | 0.183 | 0.192 | Flat |

**Key cross-construct findings (33 targets vs. 24 at 2026-04-10):**

1. **PLXA4 — confirmed rising trend.** D1-D3 0.240 → D1-D6 0.270 → D1-D7 0.338. Plexin-A4 is a semaphorin receptor and a biologically plausible axon guidance candidate. D1-D7 mean sits ~1.7 SD above the D1-D7 background median (0.240). The monotonic cross-construct rise is the strongest non-control pattern in the dataset and is stable across the 9 additional targets added since the last report.

2. **LRP12 — rising trend confirmed.** D1-D3 0.252 → D1-D6 0.245 → D1-D7 0.303. LDL receptor-related protein 12; expressed in brain vasculature. Modest but consistent with the PLXA4 pattern (engagement requires domains beyond D1-D3).

3. **KIRR3 — rising from background.** 0.182 → 0.183 → 0.279. Ig-superfamily synapse protein, a clean "emerges only in the full ectodomain" pattern.

4. **MER (MERTK) — new rising entrant.** 0.259 → 0.195 → 0.296. TAM receptor tyrosine kinase; non-monotonic but D1-D7 exceeds D1-D3. Watch rather than elevate.

5. **Neuroligin degradation confirmed and broader.** NLGN1/1.1/1.2 and NLGNX all fall monotonically with construct length. Signature of a D1-D3-specific template artifact rather than true interaction. These should be deprioritized.

6. **BASI — now clearly falling** (0.412 → 0.371 → 0.284). Previously flagged as a strong D1-D3 hit; the full-construct signal collapses, matching the neuroligin artifact pattern.

---

## Interpretation

At 55% D1-D3 coverage and 16-17% D1-D6/D1-D7 coverage, the screen continues to produce a **clear negative result for high-confidence novel binary interactions**. No non-control target exceeds mean 0.50 in D1-D3. The background distribution is stable (median 0.253 vs 0.258 at n=87) — additional data has not shifted the background, confirming the previous distributional characterization.

The **positive control (VEGFA, mean 0.790 in D1-D3)** validates that the pipeline can detect a true interaction, placing an empirical ceiling on what a "real hit" looks like in this setup.

The cross-construct signal remains the most actionable layer:
- **PLXA4** is the strongest biologically plausible rising hit and is now supported by 33-target cross-construct N.
- **LRP12, KIRR3** confirm weaker but consistent rising trends.
- **Neuroligins and BASI** are confirmed template artifacts (D1-D3-specific).

---

## Caveats

- **D1-D3:** 108/196 scored (55%). 88 remain.
- **D1-D6/D1-D7:** 33 and 32 of 196 scored (~17%). Cross-construct interpretation still sensitive to which targets finish next.
- **84 targets** failed sequence fetch in prep — excluded from all constructs.
- **Binary screen only.** VEGF-bridged ternary interactions (NRP1, NRP2) are invisible to pairwise prediction and require a separate trimer run.
- **Template bias.** sFLT1 Ig-like domains have abundant PDB templates; poorly-templated partners can score artificially low.
- **Metric note.** Ranking uses `iptm+ptm` combined confidence from ranking_debug.json (same metric source as prior report). Absolute values are not directly comparable to pure ipTM from published screens.

---

## Addendum (2026-04-13 afternoon): Re-analysis with interface pLDDT

### Pipeline bug discovered

The post-AF2 extractor `09_extract_and_cleanup.py` was never invoked for any fullscreen job because of a path-resolution bug in the wrapper: it looked for `${WORK_DIR}/analysis/03_structural_prediction/09_extract_and_cleanup.py`, where `WORK_DIR=.../d1d{3,6,7}_fullscreen/` — three levels below the script's real location. The wrapper silently fell through to `WARNING: extract_and_cleanup.py not found, deleting pkl files directly`. As a result **every completed job wrote only `ranking_debug.json`, with the result_model_*.pkl files deleted before scores.json / pae_matrices.npz could be written.** We irretrievably lost: pure ipTM (separate from pTM), PAE matrices, interchain_pae, ipSAE, and LIS on all 175 completed targets.

**Fix applied:** `03_generate_lsf_jobs.py` now derives `PROJECT_ROOT="${WORK_DIR%/results/*}"` and resolves the extractor via project root. All 666 pending `.lsf` files were patched in place (base64 decode → replace → re-encode). Future jobs will produce `scores.json` + `pae_matrices.npz` correctly.

### What we recovered from surviving PDBs

`analysis/03_structural_prediction/12_recover_interface_from_pdb.py` computes from `ranked_*.pdb` B-factor (pLDDT) column and Cβ–Cβ 8 Å interface detection:
- `interface_plddt_joint_mean` — mean pLDDT over interface residues, averaged across 25 models
- `n_interface_residues_a/b` — interface residue counts per chain
- Whole-chain mean pLDDT per chain (sanity)

PAE-derived metrics cannot be recovered — they require the deleted pkl files.

### The ranking changes almost completely

**Pearson correlation between iptm+ptm and interface_plddt_joint on D1-D3 (n=110): r = 0.201.**

The two metrics rank targets nearly independently. **Every finding in the prior sections of this report was based on a metric that is essentially uncorrelated with actual interface quality.** The distributional analysis ("no hit exceeds mean 0.50") remains valid as a statement about the old metric, but it should not be interpreted as evidence of no interaction.

### Positive control (VEGFA) behavior

VEGFA shows the expected control signal on **both** metrics:
- D1-D3: iPLDDT 91.8 (4.3 SD above median 64.6) | iptm+ptm 0.790
- D1-D6: iPLDDT 89.7 | iptm+ptm 0.605
- D1-D7: iPLDDT 89.3 | iptm+ptm 0.561

Interface pLDDT is much more stable across constructs (only −2.5 from D1-D3 to D1-D7) than iptm+ptm (−0.23), consistent with interface pLDDT being a local metric unaffected by the growing unengaged flank.

### New D1-D7 top hits (ranked by interface pLDDT)

| Rank | Target | iPLDDT | iptm+ptm | Biological note |
|------|--------|--------|----------|-----------------|
| ctrl | **VEGFA** | **89.3** | 0.561 | Positive control |
| 1 | **PLXA4** | **73.3** | 0.338 | Plexin-A4, semaphorin receptor ✓ |
| 2 | **NLGN2_ECD** | 73.1 | 0.319 | Neuroligin-2 ectodomain |
| 3 | **PLXA1** | **70.3** | 0.257 | Plexin-A1 (PLXA4 paralog) ✓✓ |
| 4 | NLGNX | 70.2 | 0.240 | Neuroligin |
| 5 | NLGN1.2 | 69.9 | 0.235 | |
| 6 | NOE1 | 69.8 | 0.265 | Olfactomedin-1 |
| 7 | NLGN1.1 | 69.7 | 0.237 | |
| 8 | NLGN1 | 69.6 | 0.235 | |
| 9 | LRFN4 | 69.0 | 0.220 | LRR/Ig, synapse |
| 10 | NGL1 | 67.9 | 0.235 | Netrin-G ligand |

All 10 top non-control hits sit 0.7–1.5 SD above the D1-D7 non-control median (59.7) but 16 pts below VEGFA (89.3). This is a **clear three-tier separation**: positive control (~90) → candidate tier (~68–73) → background (~60).

### Cross-construct iPLDDT table (33 targets with all 3 constructs)

| Target | D1-D3 | D1-D6 | D1-D7 | Trend |
|--------|-------|-------|-------|-------|
| **VEGFA (ctrl)** | **91.8** | **89.7** | **89.3** | Stable |
| **PLXA4** | 66.2 | 60.3 | **73.3** | **Rising** |
| **NLGN2_ECD** | 67.8 | 72.0 | **73.1** | **Rising** |
| **PLXA1** | **75.5** | 72.6 | 70.3 | High & stable |
| NLGNX | 68.1 | 72.4 | 70.2 | Flat-high |
| NLGN1.2 | 68.7 | 72.8 | 69.9 | Flat-high |
| **NOE1** | 61.5 | 62.3 | **69.8** | **Rising** |
| NLGN1.1 | 68.7 | 73.1 | 69.7 | Flat-high |
| NLGN1 | 69.1 | 73.0 | 69.6 | Flat-high |
| LRFN4 | 72.0 | 69.3 | 69.0 | Flat-high |
| NGL1 | 70.7 | 71.2 | 67.9 | Flat-high |
| Amyloid-like_protein_1 | 60.4 | 59.2 | 62.9 | Flat |
| BAI1 | 63.3 | 63.7 | 62.7 | Flat |
| BAI3 | 64.0 | 64.6 | 62.3 | Flat |
| KIRR3 | 65.1 | 50.6 | 59.7 | Noisy |
| **MER** | 73.4 | 53.9 | **48.9** | **Collapsing** |
| **BASI** | 74.7 | 59.3 | **46.9** | **Collapsing** |
| GPC5 | 68.6 | 67.9 | 53.0 | Falling |
| Cadherin-12_ECD | 69.5 | 62.2 | 41.0 | Falling |

(background median D1-D7 iPLDDT ≈ 60; SD ≈ 9)

### Reinterpretation of prior findings

1. **Plexin-A signal is the real headline.** Both PLXA4 (rising 66→60→73) and PLXA1 (high-stable 76/73/70) sit in the candidate tier across all constructs. PLXA1 was entirely invisible in the old metric (iptm+ptm 0.243, median territory). These are the closest paralogs in the plexin-A family and both score well — this is the strongest biologically coherent signal in the dataset. **PLXA1 must be elevated to co-primary candidate alongside PLXA4.**

2. **Neuroligins are not template artifacts.** The declining iptm+ptm with construct length was a size-dilution artifact of the combined confidence metric, not evidence against interaction. NLGN1 family members and NLGN2_ECD all hold stable interface pLDDT 68–73 across D1-D3/D1-D6/D1-D7. NLGN2_ECD is actively rising (67.8→72.0→73.1). The previous recommendation to "deprioritize neuroligins" was wrong and is retracted.

3. **MER and BASI are the real artifacts.** Both were strong D1-D3 candidates under the old metric and under interface pLDDT, but collapse to 46–49 in D1-D7 — below background. These are D1-D3-specific and should be deprioritized.

4. **TESC is downgraded.** It was the old metric's #1 hit (iptm+ptm 0.464) but interface pLDDT is only 69.9 (rank 15 in D1-D3), and no D1-D6/D1-D7 data available yet. Probably a compact low-confidence contact, not a real interface.

5. **LRFN4, NGL1, NOE1, and NLGN2_ECD enter the candidate tier** for the first time. LRFN4 (LRR/Ig synapse adhesion), NGL1 (netrin-G ligand-1), NOE1 (olfactomedin-1), and NLGN2_ECD are all biologically plausible for axon/synapse contact biology and warrant follow-up.

### Caveats on interface pLDDT

- **Interface pLDDT alone is not sufficient** to call a hit. AF2 can produce locally confident structures that are wrong. The metric that actually discriminates true binders is ipSAE or interchain PAE — both require pkl files we no longer have for completed jobs.
- **The three-tier separation** (control ~90 / candidates ~68–73 / background ~60) is encouraging but should be treated as hypothesis-generating, not confirmatory.
- **Re-running key candidates** (PLXA1, PLXA4, NLGN2_ECD, NOE1, NGL1, LRFN4) with the **fixed** wrapper would cost ~6 × 3 constructs × ~25 GPU-hours = ~450 GPU-hours to recover ipSAE/LIS/PAE matrices for the shortlist. This is feasible and is the correct next step before any wet-lab follow-up.

---

## What Comes Next (revised)

1. **Re-run the candidate shortlist with the fixed wrapper** to recover ipSAE, LIS, and PAE matrices:
   - Tier 1 (highest priority): PLXA1, PLXA4, NLGN2_ECD, NOE1
   - Tier 2: NGL1, LRFN4, NLGN1, UNC5H4
   - Tier 3 (sanity): re-run VEGFA once to compare old vs new extractor output
   - Estimated cost: 9 targets × 3 constructs × ~20 GPU-hours ≈ 540 GPU-hours
2. **Continue the fullscreen runs** — fixed wrapper will write scores.json + pae_matrices.npz for all remaining jobs automatically. Outstanding: 88 D1-D3, 163 D1-D6, 164 D1-D7.
3. **PLXA1+PLXA4 interface residue mapping** from the existing 25 D1-D7 ranked PDBs (no re-run needed). Identify which sFLT1 domains (D1, D2, D3, D4…) engage the plexin sema domain.
4. **Deprioritize: MER, BASI, TESC** — construct-length collapse.
5. **Trimer predictions** (deferred, unchanged) — sFLT1 + VEGFA + NRP1/NRP2/SEMA3A.
