# 01_load_helpers.R -- Shared loading and normalization functions for sFLT1 cohort data
#
# Functions:
#   load_lm_results(path, cohort)   -- Read LM CSV, normalize columns
#   load_annotation(path, cohort)   -- Read annotation CSV, normalize to snake_case
#   build_cohort_registry(data_dir) -- Build tibble mapping all 16 LM files
#   find_annotation_files(data_dir) -- Locate annotation files per cohort

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})


# --- Column name mappings ----

# Canonical output columns for LM results
LM_COLS <- c("target", "beta", "pvalue", "padj", "minus_log_padj", "direction", "cohort")

# Map raw annotation headers -> canonical snake_case
ANNOT_RENAME <- c(
  "Target"            = "target",
  "target"            = "target",
  "TargetFullName"    = "target_full_name",
  "target_full_name"  = "target_full_name",
  "UniProt"           = "uniprot",
  "uni_prot"          = "uniprot",
  "EntrezGeneID"      = "entrez_gene_id",
  "entrez_gene_id"    = "entrez_gene_id",
  "EntrezGeneSymbol"  = "entrez_gene_symbol",
  "entrez_gene_symbol" = "entrez_gene_symbol"
)


# --- Loading functions ----

#' Load a single LM results CSV
#'
#' Detects and drops a leading row-index column (empty header "").
#' Normalizes to canonical output schema.
#'
#' @param path Character. File path to CSV.
#' @param cohort Character. Cohort name to tag rows.
#' @return Tibble with columns: target, beta, pvalue, padj, minus_log_padj, direction, cohort
load_lm_results <- function(path, cohort) {
  stopifnot(file.exists(path))
  raw <- readr::read_csv(path, show_col_types = FALSE)

  # Drop leading row-index column (header is "" or "X1" or ...1)
  first_col <- names(raw)[1]
  if (first_col == "" || grepl("^\\.{3}\\d+$", first_col) || first_col == "X1") {
    raw <- raw[, -1, drop = FALSE]
  }

  # Validate required columns
  required <- c("Target", "beta", "pvalue", "padj", "minus_log_padj", "Significant_Direction")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop(sprintf("File %s missing columns: %s", basename(path), paste(missing, collapse = ", ")))
  }

  tibble::tibble(
    target         = as.character(raw$Target),
    beta           = as.numeric(raw$beta),
    pvalue         = as.numeric(raw$pvalue),
    padj           = as.numeric(raw$padj),
    minus_log_padj = as.numeric(raw$minus_log_padj),
    direction      = as.character(raw$Significant_Direction),
    cohort         = cohort
  )
}


#' Load a protein annotation CSV
#'
#' Handles: row-index column, snake_case vs CamelCase, UCSF_AD extras (SeqId, SomaId),
#' single-column neg files (header "x").
#'
#' @param path Character. File path to annotation CSV.
#' @param cohort Character. Cohort name.
#' @return Tibble with canonical columns: target, target_full_name, uniprot,
#'         entrez_gene_id, entrez_gene_symbol. Single-column neg files return
#'         tibble with only `target`.
load_annotation <- function(path, cohort) {
  stopifnot(file.exists(path))
  raw <- readr::read_csv(path, show_col_types = FALSE)

  # Handle single-column neg files (header "x")
  if (ncol(raw) == 1 && names(raw)[1] == "x") {
    return(tibble::tibble(target = as.character(raw$x)))
  }

  # Drop leading row-index column
  first_col <- names(raw)[1]
  if (first_col == "" || grepl("^\\.{3}\\d+$", first_col) || first_col == "X1") {
    raw <- raw[, -1, drop = FALSE]
  }

  # Drop UCSF_AD extras
  raw <- raw[, !names(raw) %in% c("SeqId", "SomaId"), drop = FALSE]

  # Rename to canonical
  old_names <- names(raw)
  new_names <- old_names
  for (i in seq_along(old_names)) {
    if (old_names[i] %in% names(ANNOT_RENAME)) {
      new_names[i] <- ANNOT_RENAME[old_names[i]]
    }
  }
  names(raw) <- new_names

  # Select canonical columns (tolerate missing)
  canon <- c("target", "target_full_name", "uniprot", "entrez_gene_id", "entrez_gene_symbol")
  present <- intersect(canon, names(raw))
  out <- raw[, present, drop = FALSE]

  # Coerce types
  if ("target" %in% names(out)) out$target <- as.character(out$target)
  if ("uniprot" %in% names(out)) out$uniprot <- as.character(out$uniprot)
  if ("entrez_gene_id" %in% names(out)) {
    out$entrez_gene_id <- suppressWarnings(as.integer(out$entrez_gene_id))
  }
  if ("entrez_gene_symbol" %in% names(out)) out$entrez_gene_symbol <- as.character(out$entrez_gene_symbol)

  tibble::as_tibble(out)
}


#' Build a registry of all 16 LM result files
#'
#' Programmatically discovers files under data_dir/{cohort}/ and maps them to
#' cohort, role, somamer, direction.
#'
#' @param data_dir Character. Path to the data/ directory.
#' @return Tibble with columns: cohort, role, somamer, direction, file_path
build_cohort_registry <- function(data_dir) {
  cohorts <- list(
    list(name = "MarkVCID", role = "discovery", suffix = ""),
    list(name = "UCSF_AD", role = "discovery", suffix = "_age"),
    list(name = "GNPC",    role = "validation", suffix = ""),
    list(name = "WASHU",   role = "validation", suffix = "")
  )

  somamers   <- c("VEGFsR1", "VEGFsR1.1")
  directions <- c("pos", "neg")

  registry <- tibble::tibble(
    cohort    = character(),
    role      = character(),
    somamer   = character(),
    direction = character(),
    file_path = character()
  )

  for (coh in cohorts) {
    for (som in somamers) {
      for (dir in directions) {
        fname <- sprintf("values_LM_%s_sig_%s%s.csv", som, dir, coh$suffix)
        fpath <- file.path(data_dir, coh$name, fname)
        registry <- dplyr::bind_rows(registry, tibble::tibble(
          cohort    = coh$name,
          role      = coh$role,
          somamer   = som,
          direction = dir,
          file_path = fpath
        ))
      }
    }
  }
  registry
}


#' Find annotation files per cohort
#'
#' Returns a tibble mapping cohort, direction -> annotation file path.
#' GNPC has pos-only annotation. UCSF_AD/WASHU neg are single-column files.
#'
#' @param data_dir Character. Path to the data/ directory.
#' @return Tibble with columns: cohort, direction, annot_path
find_annotation_files <- function(data_dir) {
  annot_map <- list(
    list(cohort = "MarkVCID", direction = "pos",
         file = "discovery_vasc_markvcid_pos_proteins_common.csv"),
    list(cohort = "MarkVCID", direction = "neg",
         file = "discovery_vasc_markvcid_neg_proteins_common.csv"),
    list(cohort = "UCSF_AD", direction = "pos",
         file = "discovery_ucsfad_pos_proteins_common_adj_age.csv"),
    list(cohort = "UCSF_AD", direction = "neg",
         file = "discovery_ucsfad_neg_proteins_common_adj_age.csv"),
    list(cohort = "GNPC", direction = "pos",
         file = "validation_gnpc_pos_proteins_common.csv"),
    # GNPC neg annotation does not exist
    list(cohort = "WASHU", direction = "pos",
         file = "validation_washu_pos_proteins_common_between_the_2_somamers_adj_for_age.csv"),
    list(cohort = "WASHU", direction = "neg",
         file = "validation_washu_neg_proteins_common_between_the_2_somamers_adj_for_age.csv")
  )

  out <- tibble::tibble(
    cohort    = character(),
    direction = character(),
    annot_path = character()
  )
  for (entry in annot_map) {
    fpath <- file.path(data_dir, entry$cohort, entry$file)
    if (file.exists(fpath)) {
      out <- dplyr::bind_rows(out, tibble::tibble(
        cohort     = entry$cohort,
        direction  = entry$direction,
        annot_path = fpath
      ))
    }
  }
  out
}


#' Build a merged annotation table from all available annotation files
#'
#' Priority: MarkVCID -> UCSF_AD -> GNPC -> WASHU.
#' Performs exact match on target name, then fuzzy match stripping .N suffixes.
#'
#' @param data_dir Character. Path to the data/ directory.
#' @return Tibble with unique target -> annotation mapping
build_annotation_table <- function(data_dir) {
  annot_files <- find_annotation_files(data_dir)

  # Priority order
  priority_order <- c("MarkVCID", "UCSF_AD", "GNPC", "WASHU")

  all_annot <- list()
  for (coh in priority_order) {
    coh_files <- dplyr::filter(annot_files, cohort == coh)
    for (i in seq_len(nrow(coh_files))) {
      a <- load_annotation(coh_files$annot_path[i], coh)
      if (ncol(a) > 1) {
        all_annot <- c(all_annot, list(a))
      }
    }
  }

  if (length(all_annot) == 0) return(tibble::tibble())

  combined <- dplyr::bind_rows(all_annot)

  # Keep first occurrence per target (priority order preserved)
  combined <- combined[!duplicated(combined$target), ]

  tibble::as_tibble(combined)
}


#' Join targets with annotation, with fuzzy fallback
#'
#' First tries exact match on target. For misses, strips trailing .N suffix
#' (e.g., "14-3-3 eta.1" -> "14-3-3 eta") and retries.
#'
#' @param targets Character vector of target names.
#' @param annot_table Tibble from build_annotation_table().
#' @return Tibble with target + annotation columns. NA for unmatched.
join_annotation <- function(targets, annot_table) {
  tbl <- tibble::tibble(target = unique(targets))

  # Exact join
  joined <- dplyr::left_join(tbl, annot_table, by = "target")

  # Fuzzy fallback for misses: strip .N suffix

  missing_mask <- is.na(joined$uniprot)
  if (any(missing_mask)) {
    missing_targets <- joined$target[missing_mask]
    stripped <- stringr::str_replace(missing_targets, "\\.\\d+$", "")

    # Build lookup with stripped names
    annot_stripped <- annot_table
    annot_stripped$target_stripped <- stringr::str_replace(annot_stripped$target, "\\.\\d+$", "")

    fuzzy_join <- tibble::tibble(
      target = missing_targets,
      target_stripped = stripped
    ) |>
      dplyr::left_join(
        annot_stripped |> dplyr::select(-target),
        by = "target_stripped"
      ) |>
      dplyr::select(-target_stripped)

    # Fill in the missing rows
    annot_cols <- setdiff(names(annot_table), "target")
    for (col in annot_cols) {
      if (col %in% names(fuzzy_join)) {
        joined[[col]][missing_mask] <- fuzzy_join[[col]]
      }
    }
  }

  joined
}
