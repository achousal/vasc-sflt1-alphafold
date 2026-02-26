# 03_plot_enrichment.R -- Visualization for pathway enrichment results
#
# Generates:
#   - Dot plots: top 20 terms per database per direction
#   - Combined dot plot faceted by database for positive consensus
#   - Highlight plot: axon/semaphorin/neuropilin/plexin/VEGF/angiogenesis terms
#   - Concept network plot (cnetplot) for top GO:BP terms

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(enrichplot)
})

# Regex for axon guidance / semaphorin / neuropilin / VEGF pathway terms
HIGHLIGHT_REGEX <- "axon|semaphorin|neuropilin|plexin|VEGF|angiogenesis|slit|robo|ephrin|netrin|guidance"


#' Generate a dot plot for a single enrichment result
#'
#' @param enrich_result enrichResult object or NULL.
#' @param title Character. Plot title.
#' @param results_dir Character. Output directory.
#' @param filename Character. Output filename.
#' @param n_show Integer. Max terms to show.
plot_dotplot <- function(enrich_result, title, results_dir, filename, n_show = 20) {
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    message(sprintf("  Skipping dotplot %s: no significant terms", filename))
    return(invisible(NULL))
  }

  p <- enrichplot::dotplot(enrich_result, showCategory = n_show) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )

  n_terms <- min(n_show, nrow(as.data.frame(enrich_result)))
  height <- max(5, n_terms * 0.3 + 2)

  ggplot2::ggsave(
    file.path(results_dir, filename),
    p, width = 10, height = height, limitsize = FALSE
  )
  message(sprintf("  Saved: %s", filename))
}


#' Generate individual dot plots for all databases
#'
#' @param results Named list from run_all_enrichments().
#' @param list_name Character. Gene list name.
#' @param results_dir Character. Output directory.
plot_all_dotplots <- function(results, list_name, results_dir) {
  for (db_name in names(results)) {
    if (is.null(results[[db_name]])) next
    title <- sprintf("%s enrichment: %s", db_name, list_name)
    fname <- sprintf("step02_dotplot_%s_%s.pdf", db_name, list_name)
    plot_dotplot(results[[db_name]], title, results_dir, fname)
  }
}


#' Generate combined dot plot faceted by database
#'
#' Combines top terms from all databases into a single faceted plot.
#'
#' @param results Named list from run_all_enrichments().
#' @param list_name Character. Gene list name.
#' @param results_dir Character. Output directory.
#' @param n_per_db Integer. Top terms per database.
plot_combined_dotplot <- function(results, list_name, results_dir, n_per_db = 10) {
  dfs <- list()
  for (db_name in names(results)) {
    df <- enrich_to_df(results[[db_name]], db_name, list_name)
    if (nrow(df) > 0) {
      df <- df |>
        dplyr::arrange(p.adjust) |>
        dplyr::slice_head(n = n_per_db)
      dfs <- c(dfs, list(df))
    }
  }

  if (length(dfs) == 0) {
    message("  Skipping combined dotplot: no significant terms")
    return(invisible(NULL))
  }

  combined <- dplyr::bind_rows(dfs) |>
    dplyr::mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      })
    )

  p <- ggplot2::ggplot(
    combined,
    ggplot2::aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num),
                 size = Count, color = p.adjust)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "#B2182B", high = "#2166AC", name = "p.adjust") +
    ggplot2::facet_wrap(~ database, scales = "free_y") +
    ggplot2::labs(
      title = sprintf("Combined enrichment: %s (top %d per DB)", list_name, n_per_db),
      x = "Gene Ratio", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )

  n_terms <- nrow(combined)
  height <- max(8, n_terms * 0.25 + 3)
  n_facets <- length(unique(combined$database))
  width <- max(10, n_facets * 5)

  ggplot2::ggsave(
    file.path(results_dir, sprintf("step02_dotplot_combined_%s.pdf", list_name)),
    p, width = width, height = height, limitsize = FALSE
  )
  message(sprintf("  Saved: step02_dotplot_combined_%s.pdf", list_name))
}


#' Generate highlight plot for axon/semaphorin/neuropilin terms
#'
#' Filters enrichment results for terms matching the highlight regex
#' and produces a focused bar chart.
#'
#' @param results Named list from run_all_enrichments().
#' @param list_name Character. Gene list name.
#' @param results_dir Character. Output directory.
#' @return Tibble of highlighted terms (for stats reporting).
plot_highlight <- function(results, list_name, results_dir) {
  dfs <- list()
  for (db_name in names(results)) {
    df <- enrich_to_df(results[[db_name]], db_name, list_name)
    if (nrow(df) > 0) dfs <- c(dfs, list(df))
  }

  if (length(dfs) == 0) {
    message("  No enrichment results to highlight")
    return(tibble::tibble())
  }

  combined <- dplyr::bind_rows(dfs)

  # Filter for highlight terms
  highlighted <- combined |>
    dplyr::filter(grepl(HIGHLIGHT_REGEX, Description, ignore.case = TRUE))

  if (nrow(highlighted) == 0) {
    message("  No axon/semaphorin/neuropilin terms found in enrichment results")
    message("  (This is an explicit null result, not an error)")
    return(tibble::tibble())
  }

  message(sprintf("  Found %d axon/semaphorin-related terms", nrow(highlighted)))

  highlighted <- highlighted |>
    dplyr::mutate(
      neg_log_padj = -log10(p.adjust),
      label = sprintf("%s [%s]", Description, database)
    ) |>
    dplyr::arrange(desc(neg_log_padj))

  p <- ggplot2::ggplot(
    highlighted,
    ggplot2::aes(x = neg_log_padj, y = reorder(label, neg_log_padj), fill = database)
  ) +
    ggplot2::geom_col(show.legend = TRUE) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Database") +
    ggplot2::labs(
      title = sprintf("Axon/semaphorin/neuropilin pathway terms: %s", list_name),
      x = "-log10(p.adjust)", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )

  height <- max(4, nrow(highlighted) * 0.4 + 2)
  ggplot2::ggsave(
    file.path(results_dir, "step02_highlight_axon_semaphorin.pdf"),
    p, width = 12, height = height, limitsize = FALSE
  )
  message("  Saved: step02_highlight_axon_semaphorin.pdf")

  highlighted
}


#' Generate concept network plot for top GO:BP terms
#'
#' @param go_bp enrichResult or NULL. GO:BP enrichment result.
#' @param list_name Character. Gene list name.
#' @param results_dir Character. Output directory.
#' @param n_show Integer. Number of terms to show.
plot_cnetplot <- function(go_bp, list_name, results_dir, n_show = 10) {
  if (is.null(go_bp) || nrow(as.data.frame(go_bp)) == 0) {
    message("  Skipping cnetplot: no GO:BP results")
    return(invisible(NULL))
  }

  p <- tryCatch({
    enrichplot::cnetplot(go_bp, showCategory = n_show, node_label = "category") +
      ggplot2::labs(title = sprintf("GO:BP concept network: %s", list_name)) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12))
  }, error = function(e) {
    message(sprintf("  cnetplot failed: %s", conditionMessage(e)))
    NULL
  })

  if (is.null(p)) return(invisible(NULL))

  ggplot2::ggsave(
    file.path(results_dir, sprintf("step02_cnetplot_GO_BP_%s.pdf", list_name)),
    p, width = 14, height = 12
  )
  message(sprintf("  Saved: step02_cnetplot_GO_BP_%s.pdf", list_name))
}
