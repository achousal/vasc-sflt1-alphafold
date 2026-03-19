# SFLT1_VEGF -- Development Notes

## Project Summary

sFLT1 (soluble VEGF receptor 1) decoy-receptor mechanism applied to axonogenesis: does sFLT1 sequester ligands needed for axonal guidance (semaphorins, neuropilins), preventing downstream signaling?

Pipeline: SomaScan cross-cohort overlap -> pathway enrichment -> AlphaFold multimer structural prediction.

## Completed

### Step 1: Cross-cohort overlap (R)
Consensus protein lists from 4 SomaScan cohorts (MarkVCID, UCSF_AD, GNPC, WASHU). Tier 1+2 proteins forwarded. UpSet plots, beta heatmap, tier summary generated.

### Step 2: Pathway enrichment (R)
ORA enrichment (GO BP/MF/CC, KEGG, Reactome) on consensus sets. Axon/semaphorin highlight plot confirms pathway signal.

### Step 3a-b: Candidate selection + AlphaFold submission
24 candidates selected (4 positive controls + 20 Tier 1 dual-somamer axon pathway proteins). FASTA generated, LSF jobs submitted.

### Key findings so far
- **VEGFA positive control: ipTM 0.816** -- pipeline calibrated correctly
- **NRP1/NRP2 false negatives (ipTM ~0.21)** -- explained by VEGF-bridged ternary interaction invisible to binary screen, domain truncation (D1-D3 only), and size asymmetry

## In Progress

### Step 3c: Result parsing + scoring (Python)
15/24 AlphaFold runs complete. 10 resubmitted with 72h wall time (2026-03-16).

| Status | Count | Targets |
|--------|-------|---------|
| Complete (25/25 pdb) | 15 | VEGFA, NRP1, NRP2, SEMA3A, STMN3, NOE1, FSTL4, Amyloid-like_protein_1, BASI, PCDH9, NRX1B.1, SLIK4, APLP2.1, Calcineurin_B_a, ROBO2 |
| Partial (wall time kill) | 3 | BAI1 (19/25), PLXA4 (20/25), PLXA1 (16/25) |
| MSA only (no predictions) | 6 | Contactin-5, MER, Dtk, NGL1, DSCAM, SLIT2 |
| Empty | 1 | ROBO2 |

**Next:**
- Run `af2_scores.py` on 15 completed to get preliminary ipTM/pTM rankings
- Build `04_parse_results.py` -- full scoring: ipTM, pTM, mean inter-chain PAE, mean interface pLDDT, N interface residues, interaction call
- Build `05_plot_results.py` -- PAE heatmaps per candidate, ipTM bar chart with VEGFA benchmark, summary table
- Calibrate PAE thresholds against VEGFA positive control (no universal threshold; per CLAUDE.md guardrails)
- Check template bias: correlate PAE with PDB template coverage for Ig-like domain inflation

## Future Work

### Subcellular localization stratification
Classify 24 partners by UniProt subcellular localization (secreted, single-pass TM, GPI-anchored, intracellular). Compare ipTM distributions across classes. sFLT1 is soluble (D1-D3, no TM), so its biology differs by compartment: free ligand trapping in plasma vs surface co-receptor engagement on endothelium/macrophages. Tests whether AF2 Multimer models surface-tethered interactions differently from soluble-soluble.

### Domain-resolved interface mapping
Extract per-residue interface contacts from AF2 predictions and map onto sFLT1 D1/D2/D3 boundaries. VEGFA canonically binds D2 -- validate pipeline recovery, then classify novel partners by which Ig domain they engage. Enables binding mode classification (see below).

### Binding mode classification
Cluster partners by predicted interface footprint on sFLT1:
- **VEGF-competitive** -- overlapping D2 interface with VEGFA
- **Non-competitive / distinct site** -- D1 or D3 engagement
- **Bridging** -- spanning multiple Ig domains
Determines whether novel partners could act as co-receptors vs competitors for VEGF binding.

### Ig-domain architecture bias test
Separate partners into Ig-domain-containing (NRP1, NRP2, PCDH9, NRX1B, Contactin-5) vs non-Ig (SEMA3A, SLIT2, ROBO2, STMN3, FSTL4). If Ig-domain partners score systematically higher, structural homology to sFLT1's own folds may inflate PAE confidence (extends template bias check).

### Macrophage surface expression cross-reference
Cross-reference high-ipTM partners against M2-macrophage surface proteome (van Aanhold 2025). Partners that are both structurally plausible AND expressed on M2 macrophages are candidate mediators of sFLT1's VEGF-independent immunomodulatory function. Supports hyp-006 dual-filter validation.

### PAE threshold parameter sweep
Systematic exploration of PAE score thresholds alongside ComBat-seq batch correction parameters. Goal: quantify how batch harmonization changes candidate rankings and assess filtering sensitivity across threshold choices.

### Trimer predictions
sFLT1 + NRP1 + VEGFA trimer to test the VEGF-bridged co-receptor hypothesis. NRP1/NRP2 false negatives in binary screen are expected if the interaction is VEGF-mediated. A successful trimer prediction (high ipTM for NRP1 chain in ternary complex) would confirm this mechanism and rescue known biology the binary screen misses.

### Crystallization propensity screening
Run XtalPred or PPCpred on shortlisted sFLT1 + candidate complexes to triage experimental validation strategy:
- X-ray crystallography: domain fragments with small partners
- Cryo-EM: full-length sFLT1 complexes
- Cross-linking MS / HDX-MS: conformational dynamics, lower resolution but broadly applicable

### Integration with SomaScan epidemiology
Overlay structural scores onto cross-cohort consensus and pathway enrichment. Final candidate table should rank by: (a) cross-cohort replication tier, (b) axonogenesis/semaphorin pathway membership, (c) structural plausibility (ipTM + PAE), (d) template bias flag.

### Upstream refinements
- ComBat-seq harmonization of raw SomaScan data (currently using pre-computed LM results)
- Comorbidity adjustment: sFLT1 associations may reflect shared vascular risk, not direct mechanism
- Background set: current ORA uses SomaScan v4.1 panel size (~7289) as approximate universe

## Data Schema Notes

Row-index column: MarkVCID/GNPC have it, UCSF_AD/WASHU don't. UCSF_AD files have `_age` suffix. GNPC uses snake_case headers. Neg annotation files for UCSF_AD/WASHU are single-column. GNPC neg has 3-5 proteins -- sparse negative overlap expected.

## Scoring Criteria

| Level | ipTM | Inter-chain PAE | Call |
|-------|------|-----------------|------|
| High confidence | > 0.8 | < 10 A | Predicted interaction |
| Moderate | > 0.6 | < 15 A | Predicted interaction |
| Low / no interaction | < 0.6 | > 15 A | Not supported |

VEGFA calibration: ipTM must exceed 0.7 or method has a calibration problem.
