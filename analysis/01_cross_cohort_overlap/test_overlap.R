# test_overlap.R -- Unit tests for Step 1: Cross-Cohort Overlap
#
# Run: Rscript -e "testthat::test_file('analysis/01_cross_cohort_overlap/test_overlap.R')"

library(testthat)

# Source helpers
script_dir <- dirname(sys.frame(1)$ofile %||% ".")
source(file.path(script_dir, "01_load_helpers.R"))
source(file.path(script_dir, "02_compute_overlap.R"))

# Resolve data dir relative to project root (2 levels up from this script)
project_root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
data_dir <- file.path(project_root, "data")

# --- Test load_lm_results() ----

test_that("load_lm_results handles files WITH row-index column", {
  skip_if_not(file.exists(file.path(data_dir, "MarkVCID/values_LM_VEGFsR1_sig_pos.csv")),
              "MarkVCID data not found")
  res <- load_lm_results(
    file.path(data_dir, "MarkVCID/values_LM_VEGFsR1_sig_pos.csv"), "MarkVCID"
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, LM_COLS)
  expect_true(all(res$cohort == "MarkVCID"))
  expect_true(all(res$direction == "Positive"))
  expect_true(is.numeric(res$beta))
  expect_true(is.numeric(res$padj))
})

test_that("load_lm_results handles files WITHOUT row-index column", {
  skip_if_not(file.exists(file.path(data_dir, "UCSF_AD/values_LM_VEGFsR1_sig_pos_age.csv")),
              "UCSF_AD data not found")
  res <- load_lm_results(
    file.path(data_dir, "UCSF_AD/values_LM_VEGFsR1_sig_pos_age.csv"), "UCSF_AD"
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, LM_COLS)
  expect_true(all(res$cohort == "UCSF_AD"))
})

test_that("load_lm_results handles WASHU files (no row-index)", {
  skip_if_not(file.exists(file.path(data_dir, "WASHU/values_LM_VEGFsR1_sig_pos.csv")),
              "WASHU data not found")
  res <- load_lm_results(
    file.path(data_dir, "WASHU/values_LM_VEGFsR1_sig_pos.csv"), "WASHU"
  )
  expect_named(res, LM_COLS)
})

test_that("load_lm_results handles GNPC files (with row-index)", {
  skip_if_not(file.exists(file.path(data_dir, "GNPC/values_LM_VEGFsR1_sig_pos.csv")),
              "GNPC data not found")
  res <- load_lm_results(
    file.path(data_dir, "GNPC/values_LM_VEGFsR1_sig_pos.csv"), "GNPC"
  )
  expect_named(res, LM_COLS)
})


# --- Test load_annotation() ----

test_that("load_annotation normalizes MarkVCID CamelCase", {
  skip_if_not(file.exists(
    file.path(data_dir, "MarkVCID/discovery_vasc_markvcid_pos_proteins_common.csv")),
    "MarkVCID annotation not found")
  res <- load_annotation(
    file.path(data_dir, "MarkVCID/discovery_vasc_markvcid_pos_proteins_common.csv"),
    "MarkVCID"
  )
  expect_true("target" %in% names(res))
  expect_true("uniprot" %in% names(res))
  expect_true("entrez_gene_symbol" %in% names(res))
  expect_true(is.character(res$target))
})

test_that("load_annotation normalizes GNPC snake_case", {
  skip_if_not(file.exists(
    file.path(data_dir, "GNPC/validation_gnpc_pos_proteins_common.csv")),
    "GNPC annotation not found")
  res <- load_annotation(
    file.path(data_dir, "GNPC/validation_gnpc_pos_proteins_common.csv"),
    "GNPC"
  )
  expect_true("target" %in% names(res))
  expect_true("uniprot" %in% names(res))
})

test_that("load_annotation handles UCSF_AD extras (SeqId, SomaId dropped)", {
  skip_if_not(file.exists(
    file.path(data_dir, "UCSF_AD/discovery_ucsfad_pos_proteins_common_adj_age.csv")),
    "UCSF_AD annotation not found")
  res <- load_annotation(
    file.path(data_dir, "UCSF_AD/discovery_ucsfad_pos_proteins_common_adj_age.csv"),
    "UCSF_AD"
  )
  expect_false("SeqId" %in% names(res))
  expect_false("SomaId" %in% names(res))
  expect_true("target" %in% names(res))
})

test_that("load_annotation handles single-column neg files (header 'x')", {
  skip_if_not(file.exists(
    file.path(data_dir, "UCSF_AD/discovery_ucsfad_neg_proteins_common_adj_age.csv")),
    "UCSF_AD neg annotation not found")
  res <- load_annotation(
    file.path(data_dir, "UCSF_AD/discovery_ucsfad_neg_proteins_common_adj_age.csv"),
    "UCSF_AD"
  )
  expect_true("target" %in% names(res))
  expect_equal(ncol(res), 1)
  expect_true(nrow(res) > 0)
})


# --- Test build_cohort_registry() ----

test_that("build_cohort_registry returns 16 entries", {
  skip_if_not(dir.exists(data_dir), "Data directory not found")
  reg <- build_cohort_registry(data_dir)
  expect_equal(nrow(reg), 16)
  expect_named(reg, c("cohort", "role", "somamer", "direction", "file_path"))
  expect_equal(sum(reg$role == "discovery"), 8)
  expect_equal(sum(reg$role == "validation"), 8)
})


# --- Test compute_tiers() with toy data ----

test_that("compute_tiers assigns correct tiers", {
  # Toy scenario: protein A in all 4, B in both disc + 1 val, C in disc only
  lm_list <- list(
    MarkVCID = tibble::tibble(target = c("A", "B", "C", "D")),
    UCSF_AD  = tibble::tibble(target = c("A", "B", "C")),
    GNPC     = tibble::tibble(target = c("A", "B")),
    WASHU    = tibble::tibble(target = c("A"))
  )
  tiers <- compute_tiers(lm_list)
  expect_equal(tiers$tier[tiers$target == "A"], 1L)
  expect_equal(tiers$tier[tiers$target == "B"], 2L)
  expect_equal(tiers$tier[tiers$target == "C"], 3L)
  # D is only in MarkVCID -- should not appear

  expect_false("D" %in% tiers$target)
})

test_that("compute_tiers handles empty inputs", {
  tiers <- compute_tiers(list())
  expect_equal(nrow(tiers), 0)
})


# --- Test no duplicates ----

test_that("no duplicate targets within (somamer, direction, cohort)", {
  skip_if_not(dir.exists(data_dir), "Data directory not found")
  reg <- build_cohort_registry(data_dir)
  for (i in seq_len(nrow(reg))) {
    r <- reg[i, ]
    if (!file.exists(r$file_path)) next
    dat <- load_lm_results(r$file_path, r$cohort)
    dup_count <- sum(duplicated(dat$target))
    expect_equal(dup_count, 0,
                 info = sprintf("Duplicates in %s %s %s", r$cohort, r$somamer, r$direction))
  }
})


# --- Test full pipeline output ----

test_that("full pipeline produces valid output CSVs", {
  skip_if_not(dir.exists(data_dir), "Data directory not found")
  reg <- build_cohort_registry(data_dir)
  annot <- build_annotation_table(data_dir)
  overlaps <- compute_all_overlaps(reg, annot)

  # Required columns present
  expected_cols <- c("target", "uniprot", "entrez_gene_symbol", "entrez_gene_id",
                     "direction", "somamer", "tier", "n_cohorts", "cohorts_sig",
                     "mean_beta", "min_padj", "dual_somamer")
  expect_true(all(expected_cols %in% names(overlaps)),
              info = paste("Missing:", setdiff(expected_cols, names(overlaps))))

  # No NA in target or direction
  expect_true(all(!is.na(overlaps$target)))
  expect_true(all(!is.na(overlaps$direction)))

  # Tier values valid
  expect_true(all(overlaps$tier %in% c(1L, 2L, 3L)))

  # Consensus (Tier 1+2) annotation check
  consensus <- dplyr::filter(overlaps, tier <= 2)
  if (nrow(consensus) > 0) {
    annot_rate <- sum(!is.na(consensus$uniprot)) / nrow(consensus)
    if (annot_rate < 0.9) {
      warning(sprintf("Annotation join rate %.1f%% < 90%%", annot_rate * 100))
    }
  }
})
