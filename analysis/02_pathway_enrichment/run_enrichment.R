#!/usr/bin/env Rscript
# run_enrichment.R -- Step 2 entrypoint: Pathway Enrichment Analysis
#
# Usage:
#   Rscript run_enrichment.R [--consensus-dir DIR] [--results-dir DIR] [--background FILE]

suppressPackageStartupMessages(library(optparse))

# --- CLI ----
option_list <- list(
  optparse::make_option("--consensus-dir", type = "character",
                        default = "results/01_cross_cohort_overlap/",
                        help = "Path to Step 1 consensus output [default: %default]"),
  optparse::make_option("--results-dir", type = "character",
                        default = "results/02_pathway_enrichment/",
                        help = "Path to results output directory [default: %default]"),
  optparse::make_option("--background", type = "character", default = NULL,
                        help = "Optional background gene list CSV [default: none]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

consensus_dir <- opt[["consensus-dir"]]
results_dir   <- opt[["results-dir"]]
bg_path       <- opt$background

if (!dir.exists(consensus_dir)) {
  stop(sprintf("Consensus directory not found: %s\nRun Step 1 first.", consensus_dir))
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# --- Source helpers ----
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
source(file.path(script_dir, "03_plot_enrichment.R"))


# --- Main ----
message("=== Step 2: Pathway Enrichment ===")
message(sprintf("Consensus: %s", consensus_dir))
message(sprintf("Results:   %s", results_dir))

# 1. Build universe
message("\n[1/5] Building background universe...")
universe <- build_universe(bg_path)

# 2. Prepare gene lists
message("\n[2/5] Preparing gene lists...")
gene_lists <- build_gene_lists(consensus_dir)
mapping_stats <- get_mapping_stats(consensus_dir)
message(sprintf("  Positive mapping rate: %d/%d (%.1f%%)",
                mapping_stats$pos_mapped, mapping_stats$pos_total,
                mapping_stats$pos_rate * 100))
if (mapping_stats$neg_total > 0) {
  message(sprintf("  Negative mapping rate: %d/%d (%.1f%%)",
                  mapping_stats$neg_mapped, mapping_stats$neg_total,
                  mapping_stats$neg_rate * 100))
}

# 3. Run enrichment for each gene list
message("\n[3/5] Running enrichment...")
all_results <- list()
all_dfs <- list()

# Focus on combined lists (pos_combined, neg_combined) plus per-somamer
for (list_name in names(gene_lists)) {
  ids <- gene_lists[[list_name]]
  if (length(ids) < 3) {
    message(sprintf("\n  Skipping '%s': only %d genes", list_name, length(ids)))
    next
  }
  res <- run_all_enrichments(ids, list_name, universe)
  all_results[[list_name]] <- res
  df <- save_enrichment_results(res, list_name, results_dir)
  if (nrow(df) > 0) all_dfs <- c(all_dfs, list(df))
}

# 4. Save combined enrichment objects
message("\n[4/5] Saving enrichment objects...")
saveRDS(all_results, file.path(results_dir, "step02_enrichment_objects.rds"))
message("  Saved: step02_enrichment_objects.rds")

# 5. Generate plots
message("\n[5/5] Generating plots...")

# Dot plots for pos_combined (primary analysis)
if ("pos_combined" %in% names(all_results)) {
  plot_all_dotplots(all_results$pos_combined, "pos_combined", results_dir)
  plot_combined_dotplot(all_results$pos_combined, "pos_combined", results_dir)

  # Concept network for GO:BP
  plot_cnetplot(all_results$pos_combined$GO_BP, "pos_combined", results_dir)

  # Highlight axon/semaphorin terms
  highlighted <- plot_highlight(all_results$pos_combined, "pos_combined", results_dir)
}

# Dot plots for neg_combined (if results exist)
if ("neg_combined" %in% names(all_results)) {
  plot_all_dotplots(all_results$neg_combined, "neg_combined", results_dir)
}

# 6. Stats report
message("\nGenerating stats report...")
report_lines <- c(
  "Step 02: Pathway Enrichment -- Summary Statistics",
  sprintf("Generated: %s", Sys.time()),
  ""
)

report_lines <- c(report_lines, "== Gene list sizes ==")
for (nm in names(gene_lists)) {
  report_lines <- c(report_lines, sprintf("  %s: %d genes", nm, length(gene_lists[[nm]])))
}

report_lines <- c(report_lines, "",
  sprintf("== Mapping rates =="),
  sprintf("  Positive: %d/%d (%.1f%%)",
          mapping_stats$pos_mapped, mapping_stats$pos_total,
          mapping_stats$pos_rate * 100),
  sprintf("  Negative: %d/%d (%.1f%%)",
          mapping_stats$neg_mapped, mapping_stats$neg_total,
          if (is.na(mapping_stats$neg_rate)) 0 else mapping_stats$neg_rate * 100),
  ""
)

report_lines <- c(report_lines, "== Significant terms per database ==")
for (list_name in names(all_results)) {
  report_lines <- c(report_lines, sprintf("  --- %s ---", list_name))
  for (db_name in names(all_results[[list_name]])) {
    res <- all_results[[list_name]][[db_name]]
    n_sig <- if (!is.null(res)) nrow(as.data.frame(res)) else 0
    report_lines <- c(report_lines, sprintf("    %s: %d terms", db_name, n_sig))

    # Top 5 terms
    if (n_sig > 0) {
      df <- as.data.frame(res) |> dplyr::arrange(p.adjust)
      top5 <- utils::head(df, 5)
      for (j in seq_len(nrow(top5))) {
        report_lines <- c(report_lines, sprintf(
          "      %d. %s (padj=%.2e, %s)",
          j, top5$Description[j], top5$p.adjust[j], top5$GeneRatio[j]
        ))
      }
    }
  }
}

# Axon/semaphorin highlight summary
report_lines <- c(report_lines, "", "== Axon/semaphorin pathway highlight ==")
if (exists("highlighted") && nrow(highlighted) > 0) {
  report_lines <- c(report_lines,
                    sprintf("  %d terms matched highlight regex", nrow(highlighted)))
  for (i in seq_len(min(10, nrow(highlighted)))) {
    report_lines <- c(report_lines, sprintf(
      "  %d. %s [%s] (padj=%.2e)",
      i, highlighted$Description[i], highlighted$database[i], highlighted$p.adjust[i]
    ))
  }
} else {
  report_lines <- c(report_lines,
    "  No axon/semaphorin/neuropilin terms found (explicit null result)")
}

writeLines(report_lines, file.path(results_dir, "step02_stats_report.txt"))
message("  Saved: step02_stats_report.txt")

message("\n=== Step 2 complete ===")
