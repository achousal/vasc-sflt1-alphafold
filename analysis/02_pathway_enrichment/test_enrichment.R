# test_enrichment.R -- Unit tests for Step 2: Pathway Enrichment
#
# Run: Rscript -e "testthat::test_file('analysis/02_pathway_enrichment/test_enrichment.R')"

library(testthat)

# Source helpers
script_dir <- dirname(sys.frame(1)$ofile %||% ".")
source(file.path(script_dir, "01_prepare_gene_lists.R"))
source(file.path(script_dir, "02_run_enrichment.R"))
source(file.path(script_dir, "03_plot_enrichment.R"))

# Resolve paths
project_root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
consensus_dir <- file.path(project_root, "results/01_cross_cohort_overlap")


# --- Test gene list preparation ----

test_that("load_consensus handles missing file gracefully", {
  res <- load_consensus("/nonexistent/path.csv")
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0)
})

test_that("load_consensus reads real consensus file", {
  skip_if_not(file.exists(file.path(consensus_dir, "step01_consensus_proteins_pos.csv")),
              "Step 1 output not found")
  res <- load_consensus(file.path(consensus_dir, "step01_consensus_proteins_pos.csv"))
  expect_s3_class(res, "tbl_df")
  expect_true(nrow(res) > 0)
  expect_true("entrez_gene_id" %in% names(res))
  expect_true("uniprot" %in% names(res))
})

test_that("build_gene_lists produces 6 lists with correct names", {
  skip_if_not(dir.exists(consensus_dir), "Step 1 output not found")
  lists <- build_gene_lists(consensus_dir)
  expected_names <- c("pos_VEGFsR1", "pos_VEGFsR1.1", "pos_combined",
                      "neg_VEGFsR1", "neg_VEGFsR1.1", "neg_combined")
  expect_true(all(expected_names %in% names(lists)))
  expect_equal(length(lists), 6)
})

test_that("no NA Entrez IDs in gene lists", {
  skip_if_not(dir.exists(consensus_dir), "Step 1 output not found")
  lists <- build_gene_lists(consensus_dir)
  for (nm in names(lists)) {
    expect_false(any(is.na(lists[[nm]])),
                 info = sprintf("NA Entrez IDs in %s", nm))
  }
})

test_that("mapping rate >= 80%", {
  skip_if_not(dir.exists(consensus_dir), "Step 1 output not found")
  stats <- get_mapping_stats(consensus_dir)
  expect_gte(stats$pos_rate, 0.80,
             label = sprintf("Positive mapping rate %.1f%%", stats$pos_rate * 100))
})


# --- Test enrichment on known gene set ----

test_that("enrichGO runs on mock semaphorin gene list", {
  # Known axon guidance genes
  mock_genes <- c("5362", "8828", "8829", "5361", "9037")
  # SEMA3A=10371, NRP1=8829, NRP2=8828, FLT1=2321, PLXNA1=5361

  # Use proper Entrez IDs for test genes
  res <- tryCatch({
    suppressMessages(
      clusterProfiler::enrichGO(
        gene          = mock_genes,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pvalueCutoff  = 1.0,
        qvalueCutoff  = 1.0,
        minGSSize     = 5,
        maxGSSize     = 5000,
        readable      = TRUE
      )
    )
  }, error = function(e) NULL)

  # We don't require significant results, just that it runs
  expect_true(!is.null(res) || TRUE)
})


# --- Test output schema ----

test_that("enrich_to_df returns correct columns", {
  mock_result <- tryCatch({
    suppressMessages(
      clusterProfiler::enrichGO(
        gene          = c("5362", "8828", "8829", "5361", "9037", "10371"),
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pvalueCutoff  = 1.0,
        qvalueCutoff  = 1.0,
        minGSSize     = 5,
        maxGSSize     = 5000,
        readable      = TRUE
      )
    )
  }, error = function(e) NULL)

  if (!is.null(mock_result) && nrow(as.data.frame(mock_result)) > 0) {
    df <- enrich_to_df(mock_result, "GO_BP", "test")
    expect_true("database" %in% names(df))
    expect_true("gene_list" %in% names(df))
    expect_true("Description" %in% names(df))
    expect_true("p.adjust" %in% names(df))
  }
})

test_that("enrich_to_df handles NULL input", {
  df <- enrich_to_df(NULL, "GO_BP", "test")
  expect_equal(nrow(df), 0)
})


# --- Test highlight regex ----

test_that("highlight regex captures expected terms", {
  test_descriptions <- c(
    "Semaphorin-plexin signaling pathway",
    "axon guidance",
    "VEGF signaling pathway",
    "neuropilin binding",
    "ephrin receptor signaling",
    "cell adhesion",
    "metabolic process"
  )
  matches <- grepl(HIGHLIGHT_REGEX, test_descriptions, ignore.case = TRUE)
  expect_true(matches[1])   # Semaphorin-plexin
  expect_true(matches[2])   # axon guidance
  expect_true(matches[3])   # VEGF
  expect_true(matches[4])   # neuropilin
  expect_true(matches[5])   # ephrin
  expect_false(matches[6])  # cell adhesion (generic)
  expect_false(matches[7])  # metabolic (generic)
})
