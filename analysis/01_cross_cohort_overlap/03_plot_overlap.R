# 03_plot_overlap.R -- Visualizations for cross-cohort overlap results
#
# Generates:
#   - 4 UpSet plots (somamer x direction)
#   - Beta heatmap for Tier 1+2 consensus proteins
#   - Tier summary bar chart

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(UpSetR)
})


#' Build a binary membership matrix for UpSet plot
#'
#' @param registry Tibble from build_cohort_registry().
#' @param somamer Character. Somamer to filter.
#' @param direction Character. Direction to filter ("pos" or "neg").
#' @return A data.frame with columns: target, one logical column per cohort.
build_upset_matrix <- function(registry, somamer, direction) {
  stratum <- dplyr::filter(registry, somamer == !!somamer, direction == !!direction)
  cohort_names <- c("MarkVCID", "UCSF_AD", "GNPC", "WASHU")

  targets_list <- list()
  for (coh in cohort_names) {
    fpath <- stratum$file_path[stratum$cohort == coh]
    if (length(fpath) == 0 || !file.exists(fpath)) next
    dat <- load_lm_results(fpath, coh)
    targets_list[[coh]] <- unique(dat$target)
  }

  all_targets <- sort(unique(unlist(targets_list)))
  mat <- data.frame(target = all_targets, stringsAsFactors = FALSE)
  for (coh in cohort_names) {
    mat[[coh]] <- as.integer(mat$target %in% targets_list[[coh]])
  }
  mat
}


#' Generate a single UpSet plot
#'
#' @param mat Data.frame from build_upset_matrix().
#' @param title Character. Plot title.
#' @param results_dir Character. Output directory.
#' @param filename Character. Output filename.
plot_upset <- function(mat, title, results_dir, filename) {
  cohort_cols <- setdiff(names(mat), "target")
  if (length(cohort_cols) < 2 || nrow(mat) < 2) {
    message(sprintf("  Skipping UpSet %s: insufficient data (%d targets, %d cohorts)",
                    title, nrow(mat), length(cohort_cols)))
    return(invisible(NULL))
  }

  pdf(file.path(results_dir, filename), width = 10, height = 6)
  print(
    UpSetR::upset(
      mat[, cohort_cols, drop = FALSE],
      sets = cohort_cols,
      order.by = "freq",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Proteins per Cohort",
      text.scale = 1.3,
      main.bar.color = "#2C3E50",
      sets.bar.color = "#34495E",
      point.size = 3,
      line.size = 1
    )
  )
  # Add title via grid (UpSetR uses grid graphics, not base)
  grid::grid.text(
    title, x = 0.65, y = 0.95,
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()
  message(sprintf("  Saved: %s", filename))
}


#' Generate all 4 UpSet plots
#'
#' @param registry Tibble from build_cohort_registry().
#' @param results_dir Character. Output directory.
plot_all_upsets <- function(registry, results_dir) {
  combos <- expand.grid(
    somamer = c("VEGFsR1", "VEGFsR1.1"),
    direction = c("pos", "neg"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(combos))) {
    som <- combos$somamer[i]
    dir <- combos$direction[i]
    mat <- build_upset_matrix(registry, som, dir)
    title <- sprintf("Cross-Cohort Overlap: %s %s", som, dir)
    # Replace "." in somamer for filename safety
    fname <- sprintf("step01_upset_%s_%s.pdf",
                     gsub("\\.", "_", som), dir)
    plot_upset(mat, title, results_dir, fname)
  }
}


#' Generate beta heatmap for Tier 1+2 consensus proteins
#'
#' @param overlaps Tibble from compute_all_overlaps().
#' @param registry Tibble from build_cohort_registry().
#' @param results_dir Character. Output directory.
plot_beta_heatmap <- function(overlaps, registry, results_dir) {
  consensus <- dplyr::filter(overlaps, tier <= 2)
  if (nrow(consensus) == 0) {
    message("  Skipping beta heatmap: no Tier 1+2 proteins")
    return(invisible(NULL))
  }

  # Collect beta values for consensus targets
  beta_long <- list()
  for (i in seq_len(nrow(registry))) {
    r <- registry[i, ]
    if (!file.exists(r$file_path)) next
    dat <- load_lm_results(r$file_path, r$cohort)
    dat$somamer <- r$somamer
    dat$direction_label <- r$direction
    beta_long <- c(beta_long, list(
      dat |> dplyr::select(target, beta, cohort, somamer, direction_label)
    ))
  }
  beta_df <- dplyr::bind_rows(beta_long)

  # Filter to consensus targets
  consensus_targets <- unique(consensus$target)
  beta_df <- beta_df |>
    dplyr::filter(target %in% consensus_targets)

  if (nrow(beta_df) == 0) {
    message("  Skipping beta heatmap: no beta data for consensus targets")
    return(invisible(NULL))
  }

  # Order targets by mean absolute beta
  target_order <- beta_df |>
    dplyr::group_by(target) |>
    dplyr::summarise(mean_abs_beta = mean(abs(beta), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(desc(mean_abs_beta)) |>
    dplyr::pull(target)

  beta_df$target <- factor(beta_df$target, levels = rev(target_order))
  beta_df$cohort <- factor(beta_df$cohort,
                           levels = c("MarkVCID", "UCSF_AD", "GNPC", "WASHU"))

  p <- ggplot2::ggplot(beta_df, ggplot2::aes(x = cohort, y = target, fill = beta)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "Beta"
    ) +
    ggplot2::facet_wrap(~ somamer, scales = "free_y") +
    ggplot2::labs(
      title = "Effect sizes for replicated proteins (Tier 1+2)",
      x = "Cohort", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    )

  # Adjust height based on number of targets
  n_targets <- length(unique(beta_df$target))
  height <- max(6, n_targets * 0.2 + 2)

  ggplot2::ggsave(
    file.path(results_dir, "step01_beta_heatmap.pdf"),
    p, width = 12, height = height, limitsize = FALSE
  )
  message("  Saved: step01_beta_heatmap.pdf")
}


#' Generate tier summary bar chart
#'
#' @param overlaps Tibble from compute_all_overlaps().
#' @param results_dir Character. Output directory.
plot_tier_summary <- function(overlaps, results_dir) {
  tier_counts <- overlaps |>
    dplyr::group_by(somamer, direction, tier) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(
      tier_label = paste0("Tier ", tier),
      facet_label = paste(somamer, direction, sep = " | ")
    )

  if (nrow(tier_counts) == 0) {
    message("  Skipping tier summary: no data")
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(tier_counts,
                       ggplot2::aes(x = tier_label, y = n, fill = tier_label)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.3, size = 3.5) +
    ggplot2::facet_wrap(~ facet_label, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(
      "Tier 1" = "#1B9E77",
      "Tier 2" = "#D95F02",
      "Tier 3" = "#7570B3"
    )) +
    ggplot2::labs(
      title = "Overlap tier counts by somamer and direction",
      x = NULL, y = "Number of proteins"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    file.path(results_dir, "step01_tier_summary.pdf"),
    p, width = 10, height = 6
  )
  message("  Saved: step01_tier_summary.pdf")
}
