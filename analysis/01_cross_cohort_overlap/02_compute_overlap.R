# 02_compute_overlap.R -- Core overlap logic
#
# Computes tiered consensus protein lists across 4 cohorts for each
# (somamer, direction) combination.
#
# Tier 1: significant in all 4 cohorts
# Tier 2: both discovery (MarkVCID AND UCSF_AD) + at least 1 validation (GNPC or WASHU)
# Tier 3: both discovery only (no validation) -- reported but NOT forwarded

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


#' Compute overlap tiers for a single (somamer, direction) stratum
#'
#' @param lm_list Named list of tibbles (one per cohort) with at least `target` column.
#' @param discovery Character vector. Discovery cohort names.
#' @param validation Character vector. Validation cohort names.
#' @return Tibble with columns: target, tier, n_cohorts, cohorts_sig
compute_tiers <- function(lm_list, discovery = c("MarkVCID", "UCSF_AD"),
                          validation = c("GNPC", "WASHU")) {
  all_cohorts <- c(discovery, validation)
  present_cohorts <- intersect(names(lm_list), all_cohorts)

  if (length(present_cohorts) == 0) {
    return(tibble::tibble(
      target = character(), tier = integer(),
      n_cohorts = integer(), cohorts_sig = character()
    ))
  }

  # Build presence matrix: target x cohort
  target_cohort <- lapply(present_cohorts, function(coh) {
    tibble::tibble(target = unique(lm_list[[coh]]$target), cohort = coh)
  }) |> dplyr::bind_rows()

  # Pivot wide: one row per target, TRUE/FALSE per cohort
  presence <- target_cohort |>
    dplyr::mutate(present = TRUE) |>
    tidyr::pivot_wider(
      names_from = cohort, values_from = present, values_fill = FALSE
    )

  # Compute tier
  disc_present <- intersect(discovery, present_cohorts)
  val_present <- intersect(validation, present_cohorts)

  presence <- presence |>
    dplyr::rowwise() |>
    dplyr::mutate(
      in_all_disc = all(dplyr::c_across(dplyr::all_of(disc_present))),
      n_val = sum(dplyr::c_across(dplyr::all_of(val_present))),
      n_cohorts = sum(dplyr::c_across(dplyr::all_of(present_cohorts))),
      cohorts_sig = paste(
        present_cohorts[as.logical(dplyr::c_across(dplyr::all_of(present_cohorts)))],
        collapse = ";"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      tier = dplyr::case_when(
        n_cohorts == length(present_cohorts) & length(present_cohorts) == 4 ~ 1L,
        in_all_disc & n_val >= 1 ~ 2L,
        in_all_disc & n_val == 0 ~ 3L,
        TRUE ~ NA_integer_
      )
    ) |>
    dplyr::filter(!is.na(tier)) |>
    dplyr::select(target, tier, n_cohorts, cohorts_sig)

  presence
}


#' Compute full overlap table across all (somamer, direction) strata
#'
#' @param registry Tibble from build_cohort_registry().
#' @param annot_table Tibble from build_annotation_table().
#' @return Tibble: target, uniprot, entrez_gene_symbol, entrez_gene_id,
#'         direction, somamer, tier, n_cohorts, cohorts_sig, mean_beta,
#'         min_padj, dual_somamer
compute_all_overlaps <- function(registry, annot_table) {
  somamers   <- unique(registry$somamer)
  directions <- unique(registry$direction)

  results <- list()
  for (som in somamers) {
    for (dir in directions) {
      stratum <- dplyr::filter(registry, somamer == som, direction == dir)

      # Load all cohort LM results for this stratum
      lm_list <- list()
      beta_data <- list()
      padj_data <- list()
      for (i in seq_len(nrow(stratum))) {
        coh <- stratum$cohort[i]
        fpath <- stratum$file_path[i]
        if (!file.exists(fpath)) {
          message(sprintf("  SKIP: %s not found", fpath))
          next
        }
        dat <- load_lm_results(fpath, coh)
        lm_list[[coh]] <- dat
        beta_data[[coh]] <- dat |> dplyr::select(target, beta)
        padj_data[[coh]] <- dat |> dplyr::select(target, padj)
      }

      if (length(lm_list) == 0) next

      tiers <- compute_tiers(lm_list)
      if (nrow(tiers) == 0) next

      # Compute mean beta and min padj across significant cohorts
      all_beta <- dplyr::bind_rows(beta_data, .id = "cohort")
      all_padj <- dplyr::bind_rows(padj_data, .id = "cohort")

      summary_stats <- tiers |>
        dplyr::left_join(
          all_beta |>
            dplyr::group_by(target) |>
            dplyr::summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop"),
          by = "target"
        ) |>
        dplyr::left_join(
          all_padj |>
            dplyr::group_by(target) |>
            dplyr::summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop"),
          by = "target"
        )

      summary_stats$somamer   <- som
      summary_stats$direction <- dir
      results <- c(results, list(summary_stats))
    }
  }

  if (length(results) == 0) {
    return(tibble::tibble(
      target = character(), uniprot = character(),
      entrez_gene_symbol = character(), entrez_gene_id = integer(),
      direction = character(), somamer = character(), tier = integer(),
      n_cohorts = integer(), cohorts_sig = character(),
      mean_beta = numeric(), min_padj = numeric(), dual_somamer = logical()
    ))
  }

  all_results <- dplyr::bind_rows(results)

  # Flag dual-somamer proteins (replicated in BOTH VEGFsR1 and VEGFsR1.1)
  dual <- all_results |>
    dplyr::filter(tier <= 2) |>
    dplyr::group_by(target, direction) |>
    dplyr::summarise(n_somamers = dplyr::n_distinct(somamer), .groups = "drop") |>
    dplyr::filter(n_somamers == 2) |>
    dplyr::select(target, direction) |>
    dplyr::mutate(dual_somamer = TRUE)

  all_results <- all_results |>
    dplyr::left_join(dual, by = c("target", "direction")) |>
    dplyr::mutate(dual_somamer = dplyr::coalesce(dual_somamer, FALSE))

  # Join annotation
  if (nrow(annot_table) > 0 && "uniprot" %in% names(annot_table)) {
    annot_joined <- join_annotation(all_results$target, annot_table)
    annot_cols <- intersect(
      c("uniprot", "entrez_gene_symbol", "entrez_gene_id", "target_full_name"),
      names(annot_joined)
    )

    all_results <- all_results |>
      dplyr::left_join(
        annot_joined |> dplyr::select(target, dplyr::all_of(annot_cols)),
        by = "target"
      )
  }

  # Ensure all expected columns exist
  for (col in c("uniprot", "entrez_gene_symbol", "entrez_gene_id")) {
    if (!col %in% names(all_results)) {
      all_results[[col]] <- NA
    }
  }

  # Reorder
  all_results |>
    dplyr::select(
      target, uniprot, entrez_gene_symbol, entrez_gene_id,
      direction, somamer, tier, n_cohorts, cohorts_sig,
      mean_beta, min_padj, dual_somamer
    ) |>
    dplyr::arrange(tier, direction, somamer, desc(abs(mean_beta)))
}


#' Generate stats report text
#'
#' @param overlaps Tibble from compute_all_overlaps().
#' @param registry Tibble from build_cohort_registry().
#' @param annot_table Tibble from build_annotation_table().
#' @return Character vector of report lines.
generate_stats_report <- function(overlaps, registry, annot_table) {
  lines <- c(
    "Step 01: Cross-Cohort Overlap -- Summary Statistics",
    paste0("Generated: ", Sys.time()),
    ""
  )

  # Per-file counts
  lines <- c(lines, "== Per-file protein counts ==")
  for (i in seq_len(nrow(registry))) {
    r <- registry[i, ]
    if (file.exists(r$file_path)) {
      n <- nrow(readr::read_csv(r$file_path, show_col_types = FALSE)) - 0
      # Subtract header already handled by read_csv
      lines <- c(lines, sprintf("  %s | %s | %s: %d proteins",
                                r$cohort, r$somamer, r$direction, n))
    }
  }

  lines <- c(lines, "")

  # Tier counts
  lines <- c(lines, "== Overlap tier counts ==")
  tier_summary <- overlaps |>
    dplyr::group_by(somamer, direction, tier) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  for (i in seq_len(nrow(tier_summary))) {
    r <- tier_summary[i, ]
    lines <- c(lines, sprintf("  %s | %s | Tier %d: %d proteins",
                              r$somamer, r$direction, r$tier, r$n))
  }

  lines <- c(lines, "")

  # Dual-somamer counts
  n_dual <- sum(overlaps$dual_somamer & overlaps$tier <= 2)
  lines <- c(lines, sprintf("== Dual-somamer proteins (Tier 1+2): %d ==", n_dual))
  if (n_dual > 0) {
    dual_targets <- unique(overlaps$target[overlaps$dual_somamer & overlaps$tier <= 2])
    lines <- c(lines, paste("  ", dual_targets))
  }

  lines <- c(lines, "")

  # Annotation join rate
  consensus <- dplyr::filter(overlaps, tier <= 2)
  if (nrow(consensus) > 0) {
    n_annot <- sum(!is.na(consensus$uniprot))
    rate <- round(100 * n_annot / nrow(consensus), 1)
    lines <- c(lines, sprintf("== Annotation join rate (Tier 1+2): %d/%d (%.1f%%) ==",
                              n_annot, nrow(consensus), rate))
    if (rate < 90) {
      lines <- c(lines, "  WARNING: annotation join rate < 90%")
      missing <- consensus$target[is.na(consensus$uniprot)]
      lines <- c(lines, paste("  Missing annotation:", paste(unique(missing), collapse = ", ")))
    }
  }

  lines
}
