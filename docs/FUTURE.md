# vasc-sflt1-alphafold: Future Work Plan

Post-completion development items. Organized by priority and dependency.

---

## Phase 1: Immediate post-AF2 analysis (blocked on 6 remaining GPU jobs)

### 1.1 Domain-resolved interface mapping
**Hypothesis link:** hyp-031 (domain-resolved sFLT1 binding mode classification)

Map predicted interface residues onto sFLT1 domain boundaries:
- D1: residues 1-110, D2: 111-220, D3: 221-338
- For each of 24 targets, compute per-domain inter-chain PAE and interface residue count
- Classify: D2-engaging (VEGF-competitive) vs D1/D3-engaging (non-competitive)
- Cross-tabulate against partner subcellular localization (secreted vs membrane-anchored)

**Deliverables:**
- `06_domain_interface_analysis.py` -- domain-resolved PAE decomposition
- `step03_domain_interface_scores.csv` -- per-target per-domain metrics
- Domain PAE heatmap faceted by localization class

**Acceptance:** VEGFA interface maps to D2; membrane-anchored partners show D1/D3 preference

### 1.2 Template bias quantification
**Hypothesis link:** hyp-015 (falsification-gated dual-filter)

Quantify whether PDB template abundance confounds ipTM rankings:
- For each target, count PDB template hits from AF2 MSA output (`msas/` dir)
- Compute Spearman correlation: ipTM rank vs template count
- If r > 0.3 (p < 0.05), template bias is a significant confounder

**Deliverables:**
- `07_template_bias_analysis.py` -- template count extraction + correlation
- `step03_template_bias_report.txt`
- Scatter plot: ipTM vs template count, colored by interaction call

**Acceptance:** Spearman test reported; if r > 0.3, flag ipTM rankings as template-confounded

### 1.3 Alternative scoring metrics
**From vault:** ipSAE, LIS (Local Interaction Score), SPOC-style classifiers

Implement PAE-filtered alternatives to raw ipTM:
- **ipSAE** (Dunbrack Lab): PAE-filtered ipTM that improves TP/FP separation for disordered partners
- **LIS** (PAE <= 12A interface-focused): better for partners with substantial intrinsic disorder (syndecans, NRP1)
- Compare all three (ipTM, ipSAE, LIS) rankings against VEGFA positive control

**Deliverables:**
- `08_alternative_scores.py` -- ipSAE and LIS computation
- `step03_score_comparison.csv` -- side-by-side rankings
- Rank concordance plot

---

## Phase 2: Trimer and expanded structural predictions

### 2.1 VEGF-bridged trimer predictions
**Hypothesis link:** hyp-028 (NRP1-dependent co-receptor sequestration)

NRP1 scored ipTM 0.22 in binary screen (expected false negative -- known VEGF-bridged interaction). Run trimer predictions:
- sFLT1-D1D3 + VEGF-A + NRP1 (b1/b2 domains, residues 275-586)
- sFLT1-D1D3 + VEGF-A + NRP2 (equivalent domains)
- sFLT1-D1D3 + VEGF-A + SEMA3A

AF2 multimer supports 3-chain inputs. Walltime: ~96-120h per trimer on A100 (larger complex).

**Deliverables:**
- 3 new LSF scripts, trimer FASTA files
- ipTM comparison: binary vs trimer for NRP1/NRP2/SEMA3A
- Interface mapping: does VEGF bridge sFLT1-NRP1 contact?

**Acceptance:** NRP1 trimer ipTM > 0.5 (recovers known interaction that binary missed)

### 2.2 sFLT1 D1-D7 full ectodomain screen
**From:** hyp-031 mechanism point about D4-D7 interactions being missed

Current screen uses D1-D3 only (338 aa). Some partners may engage D4-D7.
- Re-run top 5 hits + top 5 misses with full FLT1 ectodomain (residues 1-750)
- Compare: do any low-scoring targets improve with full ectodomain?
- Caveat: 750 + partner = >1000 residues, GPU memory and time increase substantially

**Deliverables:**
- 10 full-ectodomain LSF scripts
- `step03_d1d3_vs_full_comparison.csv`

---

## Phase 3: Biological validation integration

### 3.1 M2 surfaceome dual-filter
**Hypothesis link:** hyp-015 (dual-filter), hyp-006 (template bias + macrophage validation)
**From vault:** TMT surfaceome datasets PXD032801/PXD032823/PXD032967

Cross-reference AF2 structural hits with M2 macrophage surface proteome:
- Download TMT surfaceome data from PRIDE
- Identify which AF2 candidates are present on M2 macrophage surface
- Dual-filter: AF2 structural hit AND M2-surface-expressed = high-confidence candidate
- Map against SomaScan panel to identify overlap

**Deliverables:**
- `09_m2_surfaceome_filter.py` -- PRIDE data download + cross-reference
- `step03_dual_filter_candidates.csv` -- candidates passing both filters
- Venn diagram: AF2 hits vs M2-surface vs SomaScan panel

### 3.2 sFLT1-macrophage mediation analysis
**Hypothesis link:** hyp-004 (macrophage-mediated sFLT1 axonal injury), hyp-016 (staged validation)
**From vault:** sCD163/CCL2 mediation, NfL-myeloid feedback loop

Test the causal chain: sFLT1 -> macrophage activation (sCD163, CCL2) -> axonal injury (NfL):
- Requires VascBrain SIMOA data (cad-simoa-assembly project, currently incomplete)
- Structural equation modeling: sFLT1 -> sCD163 -> NfL path decomposition
- Control for NfL-myeloid reverse causation (2025 Cell Reports finding)
- Sex-stratified analysis (sCD163 is female-specific in late-stage PD)

**Blocked on:** cad-simoa-assembly completion, VascBrain SIMOA data harmonization

### 3.3 Brain PVM/microglia sFLT1 receptor validation
**From vault:** NRP1/HSPG surface profiling via SEA-AD and Allen Brain Atlas snRNA-seq

Query existing public scRNA-seq datasets for sFLT1 binding machinery expression:
- FLT1/NRP1/NRP2/SDC1/GPC1 expression in human brain PVMs and microglia
- SEA-AD (Allen Institute), Allen Brain Cell Atlas
- If expressed: supports in vivo engagement premise (hyp-016 Stage 1)

**Deliverables:**
- `10_scrna_receptor_expression.py` -- query SEA-AD API or download
- Expression heatmap: sFLT1 receptor/co-receptor by brain cell type

---

## Phase 4: Experimental validation planning

### 4.1 Crystallization propensity screening
**From docs/dev.md planned enhancements**

Run XtalPred or PPCpred on top AF2 hits to triage validation strategy:
- Small partners (<500 aa) with high ipTM: X-ray crystallography candidates
- Large complexes: cryo-EM candidates
- Disordered partners: cross-linking MS or HDX-MS candidates

**Deliverables:**
- Ranked experimental validation priority list
- Cost/feasibility matrix per technique per candidate

### 4.2 sFLT1 concentration-effect curve
**From vault:** sFLT1 threshold-dependent endothelial damage (<20 ng/ml minimal, 20-30 threshold, >50 severe)

Brain-specific sFLT1 dose-response characterization:
- Compare renal threshold data (preeclampsia literature) vs cerebrovascular context
- Define whether VascBrain plasma sFLT1 levels fall in the threshold-dependent range
- If yes: supports non-linear mechanism in hyp-005

### 4.3 Heparanase-mediated local sFLT1 release
**From vault:** inflammatory heparanase could liberate HSPG-bound sFLT1 in perivascular niche

Test whether heparanase activity in cerebrovascular tissue correlates with free sFLT1 levels:
- Mine CADASIL/CSVD proteomics for heparanase (HPSE) and sFLT1 co-elevation
- If correlated: supports local release mechanism bypassing BBB crossing

---

## Phase 5: SomaScan upstream refinements

### 5.1 ComBat-seq batch harmonization sensitivity
**From CLAUDE.md upstream constraints**

Systematic assessment of how ComBat-seq parameters affect sFLT1 association rankings:
- Re-run cross-cohort overlap (Step 1) with and without batch correction
- Compare: which candidates are stable vs batch-sensitive?
- Parameter sweep: ComBat-seq covariates, parametric vs non-parametric

### 5.2 Negative-direction overlap analysis
**From PLAN.md known risks**

Currently only positive-direction consensus forwarded to Steps 2-3. Negative overlap was sparse (3-9 proteins in validation cohorts). Revisit:
- Report negative overlap as finding (even if empty)
- If any proteins replicate negatively, add as AF2 candidates (sFLT1 anti-correlated = potential competitive displacement)

---

## Dependency Graph

```
Phase 1 (post-AF2, no new GPU) ─────────────────────────┐
  1.1 Domain interface    ──→ Phase 2.2 (full ectodomain) │
  1.2 Template bias       ──→ Phase 3.1 (dual-filter)     │
  1.3 Alternative scores  ──→ (informs all downstream)    │
                                                           │
Phase 2 (new GPU jobs, ~2 weeks) ─────────────────────────┤
  2.1 Trimer predictions  ──→ Phase 4.1 (xtal screening)  │
  2.2 Full ectodomain     ──→ Phase 4.1                   │
                                                           │
Phase 3 (biological validation, mixed blocking) ──────────┤
  3.1 M2 surfaceome       (unblocked, PRIDE data public)  │
  3.2 Mediation analysis  (blocked: cad-simoa-assembly)   │
  3.3 Brain scRNA-seq     (unblocked, SEA-AD public)      │
                                                           │
Phase 4 (experimental planning, needs Phase 1-3 results) ─┤
Phase 5 (upstream refinement, independent) ────────────────┘
```

---

## GPU Budget Estimate

| Item | N jobs | Est. GPU-hours | Timeline |
|------|--------|---------------|----------|
| Current batch (6 remaining) | 6 | ~400 | Mar 19-23 |
| Trimers (Phase 2.1) | 3 | ~360 | 1 week |
| Full ectodomain (Phase 2.2) | 10 | ~500 | 1-2 weeks |
| **Total remaining** | **19** | **~1260** | **~3 weeks** |

All constrained to A100 (AF2 2.3.2 CUDA incompatible with H100).
