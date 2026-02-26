#!/usr/bin/env Rscript
# run_overlap.R -- Step 1 entrypoint: Cross-Cohort Overlap Analysis
#
# Usage:
#   Rscript run_overlap.R [--data-dir data/] [--results-dir results/01_cross_cohort_overlap/]

suppressPackageStartupMessages(library(optparse))

# --- CLI ----
option_list <- list(
  optparse::make_option("--data-dir", type = "character",
                        default = "data/",
                        help = "Path to data directory [default: %default]"),
  optparse::make_option("--results-dir", type = "character",
                        default = "results/01_cross_cohort_overlap/",
                        help = "Path to results output directory [default: %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

data_dir    <- opt[["data-dir"]]
results_dir <- opt[["results-dir"]]

# Resolve relative paths from script location if running from project root
if (!dir.exists(data_dir)) {
  stop(sprintf("Data directory not found: %s", data_dir))
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# --- Source helpers ----
# Resolve script directory robustly for both Rscript and source()
get_script_dir <- function() {
  # Try commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  # Fallback: current working directory + known relative path
  if (file.exists("analysis/01_cross_cohort_overlap/01_load_helpers.R")) {
    return(normalizePath("analysis/01_cross_cohort_overlap"))
  }
  "."
}
script_dir <- get_script_dir()
source(file.path(script_dir, "01_load_helpers.R"))
source(file.path(script_dir, "02_compute_overlap.R"))
source(file.path(script_dir, "03_plot_overlap.R"))


# --- Main ----
message("=== Step 1: Cross-Cohort Overlap ===")
message(sprintf("Data:    %s", data_dir))
message(sprintf("Results: %s", results_dir))

# 1. Build registry and annotation
message("\n[1/5] Building file registry...")
registry <- build_cohort_registry(data_dir)
message(sprintf("  Found %d file entries", nrow(registry)))

# Verify all files exist
missing_files <- registry$file_path[!file.exists(registry$file_path)]
if (length(missing_files) > 0) {
  warning(sprintf("Missing %d data files:\n  %s",
                  length(missing_files), paste(missing_files, collapse = "\n  ")))
}

message("\n[2/5] Building annotation table...")
annot_table <- build_annotation_table(data_dir)
message(sprintf("  %d unique protein annotations loaded", nrow(annot_table)))

# 3. Compute overlaps
message("\n[3/5] Computing overlaps...")
overlaps <- compute_all_overlaps(registry, annot_table)
message(sprintf("  %d total overlap entries (all tiers)", nrow(overlaps)))

# 4. Write outputs
message("\n[4/5] Writing output files...")

# All overlaps
readr::write_csv(overlaps, file.path(results_dir, "step01_all_overlaps.csv"))
message("  Saved: step01_all_overlaps.csv")

# Tier 1+2 consensus: positive
consensus_pos <- overlaps |>
  dplyr::filter(tier <= 2, direction == "pos")
readr::write_csv(consensus_pos, file.path(results_dir, "step01_consensus_proteins_pos.csv"))
message(sprintf("  Saved: step01_consensus_proteins_pos.csv (%d rows)", nrow(consensus_pos)))

# Tier 1+2 consensus: negative
consensus_neg <- overlaps |>
  dplyr::filter(tier <= 2, direction == "neg")
readr::write_csv(consensus_neg, file.path(results_dir, "step01_consensus_proteins_neg.csv"))
message(sprintf("  Saved: step01_consensus_proteins_neg.csv (%d rows)", nrow(consensus_neg)))

# Stats report
report_lines <- generate_stats_report(overlaps, registry, annot_table)
writeLines(report_lines, file.path(results_dir, "step01_stats_report.txt"))
message("  Saved: step01_stats_report.txt")

# 5. Plots
message("\n[5/5] Generating plots...")
plot_all_upsets(registry, results_dir)
plot_beta_heatmap(overlaps, registry, results_dir)
plot_tier_summary(overlaps, results_dir)

message("\n=== Step 1 complete ===")
