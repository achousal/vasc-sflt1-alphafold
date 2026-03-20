# Test Strategy

## Smoke tests (must pass before any analysis)

- [x] All input CSVs load without NA in key columns (Target, beta, pvalue, padj)
- [x] Column names match CLAUDE.md schema for all 4 cohorts
- [x] VEGFA positive control: ipTM > 0.7 (calibration gate)

## Unit tests (`test_structural.py`)

- [ ] Domain boundary classification with synthetic PAE matrix (known D2 engagement)
- [ ] ipSAE and LIS calculations against hand-computed examples
- [ ] Template count parser against sample MSA output
- [ ] Dual-filter logic with mock data (known overlaps)

## Integration tests

- [ ] Each new script (06-10) runs on VEGFA output only (positive control smoke)
- [ ] VEGFA domain mapping hits D2 (biological ground truth)
- [ ] VEGFA remains top-3 under all alternative scoring metrics
- [ ] End-to-end: 04_parse -> 06_domain -> 08_alt_scores on VEGFA; check CSV schema + value ranges

## Reproducibility

- [x] Fixed random seeds where applicable
- [x] All scripts: non-interactive, exit 0 on success / non-zero on failure
- [x] No hardcoded paths; all paths relative to project root or via CLI args
- [ ] All new outputs include run metadata (timestamp, input file hash, script version)

## Existing tests

```bash
# R tests
Rscript -e "testthat::test_file('analysis/01_cross_cohort_overlap/test_overlap.R')"
Rscript -e "testthat::test_file('analysis/02_pathway_enrichment/test_enrichment.R')"

# Python tests
pytest analysis/03_structural_prediction/ -v
```
