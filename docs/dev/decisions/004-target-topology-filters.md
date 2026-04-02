# ADR-004: Target topology filters for AF2 sequence preparation

Date: 2026-04-02
Status: accepted

## Context

The fullscreen AF2 batch (365 targets) requires trimming each partner protein to the region accessible to extracellular sFLT1. UniProt annotates protein topology at multiple levels: keywords (subcellular location), features (Signal, Transmembrane, Topological domain), and named domains. Three questions arose during audit:

1. Should intracellular-only proteins (cytoplasm, nucleus, Golgi, ER) be excluded?
2. How should type II TM proteins be handled (N-terminal cytoplasmic, C-terminal extracellular)?
3. Should Golgi/ER lumenal domains be treated as sFLT1-accessible?

Topology audit (13_uniprot_explorer.R / 00_uniprot_api.R) found:
- 31 proteins with only intracellular keywords
- 18 proteins with no subcellular keywords
- 29 TM proteins with inferred ECD < 50 aa (mostly type II TM Golgi enzymes)
- 19 of those 29 have large annotated **Lumenal** domains (61-725 aa)
- 4 of those 29 have annotated **Extracellular** domains (rescued by annotated ECD strategy)

## Decisions

### Keep intracellular and no-keyword proteins
SomaScan detected these in CSF. Presence may reflect non-classical secretion, cell death/leakage, or exosomal release. AF2 scores structural complementarity regardless of biological context — a high score from a "cytoplasmic" protein is informative. Intracellular proteins also serve as implicit negative controls (expect low ipTM).

### Use annotated Extracellular domains as primary ECD source
UniProt `Topological domain: Extracellular` features provide ground-truth residue boundaries. This correctly handles type II TM proteins (TWEAK, AT1B1, AT1B2) where the SP→TM inference gives the wrong side of the membrane. Fall back to inferred SP→TM only when no annotation exists.

### Drop targets with only Lumenal domains (no Extracellular annotation)
Lumenal domains face the interior of the Golgi/ER, not the extracellular space. sFLT1 is extracellular and cannot access these compartments in vivo. Proteins with only `Topological domain: Lumenal` and no `Topological domain: Extracellular` annotation are excluded from AF2 modeling.

**Affected targets (19):** H6ST3, CHSTC, B3GA1, B3GA3, MGAT3, GALT1, LARGE, HS2ST, GOLM1, GOLM1.1, IMPA3, B4GT6, CASC4, B4GT2, XXLT1, B3GT2, PLD3, ENTP6, TM38B. All are Golgi-resident glycosyltransferases or ER enzymes.

### Skip TM proteins with ECD < 50 aa and no annotated Extracellular domain
After applying the annotated ECD rescue (which saves TWEAK, AT1B1, AT1B2), remaining TM proteins with < 50 aa of extracellular surface have insufficient interface area for meaningful AF2 prediction.

**Affected targets (6):** SYT5, QPCTL, CQ062, ZCD1, CS077, F163B.

### Use largest single extracellular segment for multi-pass TM proteins
For multi-pass TM proteins (7-12 TM helices), the union span of all `Topological domain: Extracellular` segments includes TM helices and cytoplasmic loops between them. Instead, use the **largest single** extracellular segment. For adhesion GPCRs (BAI1, AGRB2, LPHN3) this is the large N-terminal ECD (850-920 aa). For SV2A (12 TM) this is the biggest extracellular loop (130 aa).

**Rationale:** AF2 given TM helices in the input will try to fold them, wasting capacity and confusing interface prediction. The largest exposed segment is the most likely interaction surface.

## Consequence

- ~25 targets dropped (19 lumenal-only + 6 tiny ECD with no annotation)
- ~340 targets modeled (from 365)
- Type II TM proteins correctly handled via annotated ECD
- Intracellular proteins retained as negative controls and for dead-cell-leaking hypothesis
- Decision is encoded in `02_fetch_sequences.py` (`select_partner_region()`) and auditable via `00_uniprot_api.R`
