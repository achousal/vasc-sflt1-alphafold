# sFLT1 AlphaFold Multimer Interaction Screen — Progress Report

**Date:** 2026-04-01
**Project:** vasc-sflt1-alphafold
**PI:** Fanny Elahi
**Analyst:** Andres Chousal

---

## Objective

Test whether sFLT1 physically interacts with axon guidance proteins beyond its known VEGF ligand, using AlphaFold Multimer predictions on candidates surfaced by cross-cohort SomaScan proteomics.

---

## Screen Design

365 candidate proteins were selected from cross-cohort SomaScan overlap (MarkVCID, UCSF_AD, GNPC, WASHU) filtered for dual-somamer replication (validation against both sFLT1 somamers is better) and axon/semaphorin pathway membership. Each candidate is screened as a binary pair with sFLT1 across three construct lengths to probe domain-specific binding:

| Construct | sFLT1 residues | Length | Purpose |
|-----------|---------------|--------|---------|
| D1-D3 | 27-330 | 304 aa | VEGF-binding domains — known interaction surface |
| D1-D6 | 27-654 | 631 aa | Adds dimerization/co-receptor domains |
| D1-D7 | 27-750 | 721 aa | Full ectodomain |

Each target gets 25 models (5 predictions x 5 model seeds). 
Useful known controls: VEGFA, NRP1, NRP2, SEMA3A.

**Compute:** AlphaFold 2.3.2, multimer preset, reduced_dbs, NVIDIA A100 (Minerva HPC).

---

## Current Completion

| Batch | Scored | Total | Remaining |
|-------|--------|-------|-----------|
| D1-D3 (full screen) | **148** | 365 | 217 |
| D1-D6 (cross-construct subset) | 11 | — | — |
| D1-D7 (cross-construct subset) | 9 | — | — |

---

## Results so far (148 targets scored)

Since best-model ipTM is inflated by stochastic outliers, **mean ipTM across all 25 models is the reliable metric**.

**Top 12 targets by mean ipTM:**

| Target | ipTM mean | ipTM best | Delta | PAE (A) | Note |
|--------|-----------|-----------|--------|---------|------|
| **VEGFA** | **0.79** | 0.80 | 0.01 | — | Positive control |
| DDAH1 | 0.56 | 0.79 | 0.23 | 12.2 | Dimethylarginine hydrolase; not axon guidance |
| KBRS1 | 0.53 | 0.75 | 0.22 | 14.0 | Kinase; not axon guidance |
| CFC1 | 0.49 | 0.60 | 0.11 | 25.1 | Cripto; low delta but high PAE |
| EphB6 | 0.47 | 0.77 | 0.29 | 24.7 | Eph receptor (kinase-dead) |
| BT2A1 | 0.45 | 0.62 | 0.17 | 22.7 | Butyrophilin |
| DLL1 | 0.44 | 0.73 | 0.29 | 27.2 | Notch ligand |
| AP2A2 | 0.44 | 0.84 | 0.40 | 17.3 | Adaptin; intracellular — biologically implausible |
| GOLM1 | 0.43 | 0.60 | 0.18 | 24.6 | Golgi membrane protein |
| FGFR-3_ECD | 0.41 | 0.75 | 0.34 | 26.0 | FGF receptor |
| CNBP1 | 0.40 | 0.75 | 0.35 | 11.1 | Copine; low PAE but high delta |
| CAH11 | 0.39 | 0.77 | 0.37 | 12.7 | Carbonic anhydrase |

**Mean ipTM distribution across 148 targets:**

| Mean ipTM range | N targets |
|-----------------|-----------|
| < 0.25 | 77 |
| 0.25-0.30 | 41 |
| 0.30-0.40 | 18 |
| 0.40-0.50 | 9 |
| 0.50-0.60 | 2 (DDAH1, KBRS1) |
| 0.60-0.80 | 1 (VEGFA) |

Median mean ipTM = 0.245. 
80% of targets have mean ipTM between 0.20 and 0.30. 
DDAH1 and KBRS1 sit modestly above background with delta > 0.20.

30/148 targets (20%) cross ipTM > 0.6 by best model — With VEGFA as highest mean > 0.60. 
The remaining 29 are driven by 1-2 stochastic outlier models out of 25. This 20% false positive rate is consistent with published AF2 multimer benchmarks (3-20% FPR for random pairs).

### Cross-construct comparison

Using mean ipTM, NRP1 shows a cross-construct trend.

| Target | D1-D3 mean | D1-D6 mean | D1-D7 mean | Trend |
|--------|-----------|-----------|-----------|-------|
| **VEGFA** | **0.79** | **0.60** | **0.56** | Consistent; expected dilution with construct length |
| **NRP1** | 0.22 | 0.25 | **0.38** | Rising — only non-VEGFA target with real trend |
| NOE1 | 0.33 | 0.27 | 0.25 | Flat at background |
| NGL1 | 0.23 | 0.23 | 0.26 | Flat at background |
| BASI | 0.21 | 0.23 | 0.24 | Flat at background |
| Contactin-5 | 0.22 | 0.24 | 0.23 | Flat at background |
| SLIK4 | 0.23 | 0.25 | 0.23 | Flat at background |
| SEMA3A | 0.24 | 0.24 | 0.24 | Flat at background |

NRP1 D1-D7 mean of 0.38 is ~1.3 SD above the D1-D3 background (mean 0.27, SD 0.08).

NRP1 is a known VEGF co-receptor, is interaction with FLT1 via D4-D7 biologically plausible?

---

## Interpretation

The D1-D3 screen at 148/365 targets is producing a clear negative result for novel binary interactions: no candidate has mean ipTM above 0.60. The next-best targets (DDAH1 at 0.56, KBRS1 at 0.53) are not axon guidance proteins and sit within the tail of the background distribution.

The cross-construct data is more interesting by mean than by best. NRP1 is the only target with a genuine mean ipTM trend across constructs (0.22 → 0.25 → 0.38), suggesting a real but weak D4-D7 interaction.

---

## Caveats

- **148 of 365 targets scored.** 217 remain. The pattern could change with more data, though the false positive rate is already stable.
- **Binary screen only.** VEGF-bridged ternary interactions (NRP1, NRP2) are invisible.
- **Template bias.** sFLT1 Ig-like domains have abundant PDB templates. Poorly-templated partners may get artificially low scores.
- **Cross-construct data is a small subset.** D1-D6 (11 targets) and D1-D7 (9 targets) need consistency analysis and more targets.
- **Alternative scoring not yet applied.** ipSAE (PAE-filtered) and LIS (Local Interaction Score) may rescue candidates that raw ipTM misses.

---

## What Comes Next

 **Complete D1-D3 screen** — 217 targets remaining. Monitor orchestrator, handle OOM failures.

 **Trimer predictions** (deferred) — sFLT1 + VEGFA + NRP1/NRP2/SEMA3A. Tests VEGF-bridged binding directly. ~360 GPU-hours.
