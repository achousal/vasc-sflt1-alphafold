#!/usr/bin/env Rscript
# run_vegf_depletion.R -- T0.4: VEGF-Depletion Inversion Test
#
# Tests whether axon/semaphorin pathway enrichment depends on VEGF-pathway
# proteins or reflects independent sFLT1 biology.
#
# Method: remove VEGF-pathway genes from the consensus positive protein list,
# re-run ORA enrichment, and quantify attenuation of axon/guidance terms.
#
# If axon/semaphorin enrichment collapses: signal is VEGF-dependent.
# If signal persists: strongest evidence for direct sFLT1 activity.
#
# Usage:
#   Rscript run_vegf_depletion.R
#   Rscript run_vegf_depletion.R --consensus-dir results/01_cross_cohort_overlap/ \
#     --results-dir results/02_pathway_enrichment/vegf_depletion/

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tibble)
})

# --- CLI ----
option_list <- list(
  optparse::make_option("--consensus-dir", type = "character",
                        default = "results/01_cross_cohort_overlap/",
                        help = "Path to Step 1 consensus output [default: %default]"),
  optparse::make_option("--results-dir", type = "character",
                        default = "results/02_pathway_enrichment/vegf_depletion/",
                        help = "Output directory for depletion results [default: %default]"),
  optparse::make_option("--original-results-dir", type = "character",
                        default = "results/02_pathway_enrichment/",
                        help = "Path to original (undepleted) enrichment results [default: %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

consensus_dir    <- opt[["consensus-dir"]]
results_dir      <- opt[["results-dir"]]
orig_results_dir <- opt[["original-results-dir"]]

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Source enrichment helpers ----
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (file.exists("analysis/02_pathway_enrichment/01_prepare_gene_lists.R")) {
    return(normalizePath("analysis/02_pathway_enrichment"))
  }
  "."
}
script_dir <- get_script_dir()
source(file.path(script_dir, "01_prepare_gene_lists.R"))
source(file.path(script_dir, "02_run_enrichment.R"))

# ============================================================================
# Define VEGF-pathway genes to remove
# ============================================================================
# Broad VEGF/angiogenesis pathway: includes VEGF ligands, VEGF receptors,
# neuropilins (VEGF co-receptors), and canonical angiogenesis factors.
# Intentionally broad to test the strongest version of the inversion.

VEGF_PATHWAY_SYMBOLS <- c(
  # VEGF ligands
  "VEGFA", "VEGFB", "VEGFC", "VEGFD", "PGF",
  # VEGF receptors
  "FLT1", "KDR", "FLT4",
  # Neuropilins (VEGF co-receptors)
  "NRP1", "NRP2",
  # Related angiogenesis receptors
  "TIE1", "TIE2", "TEK",
  # PDGF (overlapping vascular pathway)
  "PDGFRA", "PDGFRB", "PDGFA", "PDGFB",
  # FGF (angiogenic overlap)
  "FGFR1", "FGFR2", "FGFR3", "FGFR4",
  # Angiopoietins
  "ANGPT1", "ANGPT2",
  # PlGF receptor alias
  "FLT1"
)

# ============================================================================
# Step 1: Load and deplete consensus proteins
# ============================================================================
message("=== T0.4: VEGF-Depletion Inversion Test ===\n")

consensus_path <- file.path(consensus_dir, "step01_consensus_proteins_pos.csv")
if (!file.exists(consensus_path)) {
  stop(sprintf("Consensus file not found: %s", consensus_path))
}

consensus <- readr::read_csv(consensus_path, show_col_types = FALSE)

# Identify VEGF-pathway entries
vegf_mask <- consensus$entrez_gene_symbol %in% VEGF_PATHWAY_SYMBOLS
n_total <- length(unique(consensus$target))
n_vegf <- length(unique(consensus$target[vegf_mask]))
n_remaining <- length(unique(consensus$target[!vegf_mask]))

message(sprintf("Consensus proteins: %d unique targets", n_total))
message(sprintf("VEGF-pathway proteins found: %d", n_vegf))
message(sprintf("Remaining after depletion: %d", n_remaining))

# Report which VEGF genes were found
found_vegf <- sort(unique(consensus$entrez_gene_symbol[vegf_mask]))
message(sprintf("  Removed: %s", paste(found_vegf, collapse = ", ")))

not_found <- setdiff(VEGF_PATHWAY_SYMBOLS, consensus$entrez_gene_symbol)
if (length(not_found) > 0) {
  message(sprintf("  Not in consensus (already absent): %s",
                  paste(sort(not_found), collapse = ", ")))
}

# Create depleted consensus
depleted <- consensus[!vegf_mask, ]

# Save depleted list for reference
depleted_path <- file.path(results_dir, "step01_consensus_proteins_pos_vegf_depleted.csv")
readr::write_csv(depleted, depleted_path)
message(sprintf("\nSaved depleted list: %s (%d entries)\n", depleted_path, nrow(depleted)))

# ============================================================================
# Step 2: Build gene lists from depleted consensus and run enrichment
# ============================================================================
message("[1/3] Building depleted gene lists...")

# We need to create a temporary consensus dir with the depleted file
# to reuse the existing build_gene_lists function
tmp_dir <- file.path(results_dir, "tmp_consensus")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(depleted_path, file.path(tmp_dir, "step01_consensus_proteins_pos.csv"),
          overwrite = TRUE)

# Also copy negative consensus (unchanged)
neg_path <- file.path(consensus_dir, "step01_consensus_proteins_neg.csv")
if (file.exists(neg_path)) {
  file.copy(neg_path, file.path(tmp_dir, "step01_consensus_proteins_neg.csv"),
            overwrite = TRUE)
}

gene_lists <- build_gene_lists(tmp_dir)

# Clean up temp dir
unlink(tmp_dir, recursive = TRUE)

# Run enrichment on pos_combined (depleted)
message("\n[2/3] Running enrichment on VEGF-depleted gene list...")
depleted_results <- run_all_enrichments(
  gene_ids = gene_lists$pos_combined,
  list_name = "pos_combined_vegf_depleted"
)
save_enrichment_results(depleted_results, "pos_combined_vegf_depleted", results_dir)

# ============================================================================
# Step 3: Compare original vs depleted enrichment
# ============================================================================
message("\n[3/3] Comparing original vs depleted enrichment...\n")

# Load original enrichment results
axon_regex <- "axon|semaphorin|plexin|neuropilin|guidance|slit|robo"

load_enrichment_csv <- function(dir, pattern = "pos_combined") {
  files <- list.files(dir, pattern = sprintf("step02_enrichment_.*_%s\\.csv", pattern),
                      full.names = TRUE)
  if (length(files) == 0) return(tibble::tibble())
  dfs <- lapply(files, function(f) {
    tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) tibble::tibble())
  })
  dplyr::bind_rows(dfs)
}

orig_all <- load_enrichment_csv(orig_results_dir, "pos_combined")
depl_all <- load_enrichment_csv(results_dir, "pos_combined_vegf_depleted")

if (nrow(orig_all) == 0) {
  stop("No original enrichment results found. Run Step 2 first.")
}

# Filter to axon/guidance terms
orig_axon <- orig_all |>
  dplyr::filter(grepl(axon_regex, Description, ignore.case = TRUE))

depl_axon <- depl_all |>
  dplyr::filter(grepl(axon_regex, Description, ignore.case = TRUE))

# Also get all terms for denominator
message(sprintf("Original enrichment: %d total terms, %d axon/guidance terms",
                nrow(orig_all), nrow(orig_axon)))
message(sprintf("Depleted enrichment: %d total terms, %d axon/guidance terms",
                nrow(depl_all), nrow(depl_axon)))

# Join on Description to compare term-by-term
comparison <- orig_axon |>
  dplyr::select(database, Description, ID,
                orig_padj = p.adjust, orig_count = Count,
                orig_genes = geneID) |>
  dplyr::left_join(
    depl_axon |>
      dplyr::select(Description,
                     depl_padj = p.adjust, depl_count = Count,
                     depl_genes = geneID),
    by = "Description"
  ) |>
  dplyr::mutate(
    orig_log_padj = -log10(as.numeric(orig_padj)),
    depl_log_padj = ifelse(is.na(depl_padj), 0, -log10(as.numeric(depl_padj))),
    attenuation_pct = ifelse(
      orig_log_padj > 0,
      round((1 - depl_log_padj / orig_log_padj) * 100, 1),
      NA
    ),
    status = dplyr::case_when(
      is.na(depl_padj) ~ "LOST",
      as.numeric(depl_padj) > 0.05 ~ "NON-SIG",
      TRUE ~ "RETAINED"
    )
  ) |>
  dplyr::arrange(dplyr::desc(orig_log_padj))

# Save comparison
comp_path <- file.path(results_dir, "vegf_depletion_comparison.csv")
readr::write_csv(comparison, comp_path)

# ============================================================================
# Report
# ============================================================================
report_lines <- c(
  "# T0.4: VEGF-Depletion Inversion Test Results",
  sprintf("# Date: %s", Sys.Date()),
  "",
  "## Summary",
  sprintf("- Consensus proteins: %d unique targets", n_total),
  sprintf("- VEGF-pathway proteins removed: %d (%s)",
          n_vegf, paste(found_vegf, collapse = ", ")),
  sprintf("- Remaining: %d (%.1f%% of original)",
          n_remaining, n_remaining / n_total * 100),
  "",
  "## Enrichment Comparison",
  sprintf("- Original axon/guidance terms: %d", nrow(orig_axon)),
  sprintf("- Depleted axon/guidance terms: %d", nrow(depl_axon)),
  "",
  sprintf("- Terms RETAINED (padj < 0.05): %d", sum(comparison$status == "RETAINED")),
  sprintf("- Terms NON-SIG (padj > 0.05): %d", sum(comparison$status == "NON-SIG")),
  sprintf("- Terms LOST (absent from depleted): %d", sum(comparison$status == "LOST")),
  ""
)

# Top terms detail
report_lines <- c(report_lines,
  "## Top Axon/Guidance Terms (sorted by original significance)",
  sprintf("%-55s %8s %8s %6s %8s %s",
          "Term", "Orig_p", "Depl_p", "Count", "Atten%", "Status"),
  strrep("-", 100)
)

for (i in seq_len(min(nrow(comparison), 25))) {
  row <- comparison[i, ]
  orig_p <- formatC(as.numeric(row$orig_padj), format = "e", digits = 1)
  depl_p <- if (is.na(row$depl_padj)) "N/A" else formatC(as.numeric(row$depl_padj), format = "e", digits = 1)
  count_str <- sprintf("%s->%s",
                       as.character(row$orig_count),
                       ifelse(is.na(row$depl_count), "0", as.character(row$depl_count)))
  atten_str <- ifelse(is.na(row$attenuation_pct), "N/A",
                       sprintf("%.0f%%", row$attenuation_pct))
  desc_trunc <- substr(row$Description, 1, 55)
  report_lines <- c(report_lines,
    sprintf("%-55s %8s %8s %6s %8s %s",
            desc_trunc, orig_p, depl_p, count_str, atten_str, row$status)
  )
}

# Verdict
n_retained <- sum(comparison$status == "RETAINED")
n_total_axon <- nrow(comparison)
retention_rate <- if (n_total_axon > 0) n_retained / n_total_axon * 100 else 0

report_lines <- c(report_lines, "",
  "## Verdict",
  sprintf("Axon/guidance term retention rate: %.0f%% (%d/%d)",
          retention_rate, n_retained, n_total_axon)
)

if (retention_rate >= 80) {
  report_lines <- c(report_lines,
    "",
    "RESULT: Axon/guidance enrichment is ROBUST to VEGF-pathway depletion.",
    "The inversion hypothesis ('sflt1-associated reductions in axonogenesis",
    "pathway scores reflect VEGF pathway suppression rather than direct",
    "sflt1 activity') is NOT SUPPORTED.",
    "Signal persists: evidence for direct sFLT1 activity in axon pathways."
  )
} else if (retention_rate >= 50) {
  report_lines <- c(report_lines,
    "",
    "RESULT: Axon/guidance enrichment is PARTIALLY attenuated by VEGF depletion.",
    "Mixed evidence: some signal is VEGF-dependent, some appears independent.",
    "Interpretation should be cautious."
  )
} else {
  report_lines <- c(report_lines,
    "",
    "RESULT: Axon/guidance enrichment COLLAPSES after VEGF-pathway removal.",
    "The inversion hypothesis IS SUPPORTED: the signal is VEGF-dependent.",
    "This undermines the direct sFLT1 mechanism (hyp-004/005/016/028/031)."
  )
}

# Write report
report_path <- file.path(results_dir, "vegf_depletion_report.txt")
writeLines(report_lines, report_path)
cat(paste(report_lines, collapse = "\n"))
cat("\n")

message(sprintf("\nReport: %s", report_path))
message(sprintf("Comparison: %s", comp_path))
message("Done.")
