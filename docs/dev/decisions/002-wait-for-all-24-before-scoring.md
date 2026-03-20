# ADR-002: Wait for all 24 targets before final scoring

Date: 2026-03-19
Status: accepted

## Context

18/24 AF2 jobs are complete. Could run scoring pipeline now on partial results.

## Decision

Wait for all 24 to complete before running final `af2_scores.py` and downstream analysis. Preliminary analysis on 18 targets is fine for template bias (T1.2) and domain mapping (T1.1) since those will be trivially re-run.

## Consequence

Final results delayed ~4 days (DSCAM long pole ~Mar 22-23). Ensures all rankings, plots, and reports reflect the complete candidate set. Avoids publishing partial tables that need correction.
