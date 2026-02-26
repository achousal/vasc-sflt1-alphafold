# 01_prepare_gene_lists.R -- Prepare gene lists for ORA enrichment
#
# Reads consensus CSVs from Step 1, extracts Entrez IDs, rescues missing IDs
# via UniProt-to-Entrez mapping. Produces 6 gene lists:
#   pos/neg x VEGFsR1/VEGFsR1.1/combined

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})


#' Load consensus CSV and extract unique Entrez IDs
#'
#' @param path Character. Path to consensus CSV from Step 1.
#' @return Tibble with columns: target, uniprot, entrez_gene_symbol, entrez_gene_id
load_consensus <- function(path) {
  if (!file.exists(path)) {
    message(sprintf("  Consensus file not found: %s", path))
    return(tibble::tibble(
      target = character(), uniprot = character(),
      entrez_gene_symbol = character(), entrez_gene_id = integer()
    ))
  }
  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::select(target, uniprot, entrez_gene_symbol, entrez_gene_id,
                  direction, somamer, tier, dual_somamer)
}


#' Rescue missing Entrez IDs via UniProt-to-Entrez mapping
#'
#' Uses org.Hs.eg.db to map UniProt accessions to Entrez IDs for proteins
#' where entrez_gene_id is NA.
#'
#' @param df Tibble with at least uniprot and entrez_gene_id columns.
#' @return Same tibble with rescued entrez_gene_id values filled in.
rescue_entrez_ids <- function(df) {
  missing_mask <- is.na(df$entrez_gene_id) & !is.na(df$uniprot)
  if (!any(missing_mask)) return(df)

  missing_uniprot <- unique(df$uniprot[missing_mask])
  message(sprintf("  Attempting to rescue %d missing Entrez IDs via UniProt mapping...",
                  length(missing_uniprot)))

  mapped <- tryCatch({
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = missing_uniprot,
      keytype = "UNIPROT",
      columns = c("ENTREZID", "SYMBOL")
    )
  }, error = function(e) {
    message(sprintf("  UniProt mapping failed: %s", conditionMessage(e)))
    return(data.frame(UNIPROT = character(), ENTREZID = character(),
                      SYMBOL = character()))
  })

  if (nrow(mapped) == 0) return(df)

  # Deduplicate: keep first mapping per UniProt
  mapped <- mapped[!duplicated(mapped$UNIPROT), ]

  lookup <- stats::setNames(as.integer(mapped$ENTREZID), mapped$UNIPROT)
  symbol_lookup <- stats::setNames(mapped$SYMBOL, mapped$UNIPROT)

  for (i in which(missing_mask)) {
    up <- df$uniprot[i]
    if (up %in% names(lookup) && !is.na(lookup[up])) {
      df$entrez_gene_id[i] <- lookup[up]
    }
    if (is.na(df$entrez_gene_symbol[i]) && up %in% names(symbol_lookup)) {
      df$entrez_gene_symbol[i] <- symbol_lookup[up]
    }
  }

  n_rescued <- sum(!is.na(df$entrez_gene_id[missing_mask]))
  message(sprintf("  Rescued %d / %d Entrez IDs", n_rescued, length(missing_uniprot)))
  df
}


#' Build the 6 gene lists for enrichment
#'
#' @param consensus_dir Character. Directory containing Step 1 consensus CSVs.
#' @return Named list of integer vectors (Entrez IDs). Names:
#'         pos_VEGFsR1, pos_VEGFsR1.1, pos_combined,
#'         neg_VEGFsR1, neg_VEGFsR1.1, neg_combined
build_gene_lists <- function(consensus_dir) {
  pos <- load_consensus(file.path(consensus_dir, "step01_consensus_proteins_pos.csv"))
  neg <- load_consensus(file.path(consensus_dir, "step01_consensus_proteins_neg.csv"))

  # Rescue missing Entrez IDs
  pos <- rescue_entrez_ids(pos)
  neg <- rescue_entrez_ids(neg)

  # Build lists per somamer + combined
  extract_ids <- function(df, direction_val, somamer_val = NULL) {
    sub <- df
    if (!is.null(somamer_val)) {
      sub <- dplyr::filter(sub, somamer == somamer_val)
    }
    ids <- unique(sub$entrez_gene_id[!is.na(sub$entrez_gene_id)])
    as.character(ids)
  }

  lists <- list(
    pos_VEGFsR1     = extract_ids(pos, "pos", "VEGFsR1"),
    "pos_VEGFsR1.1" = extract_ids(pos, "pos", "VEGFsR1.1"),
    pos_combined     = extract_ids(pos, "pos"),
    neg_VEGFsR1     = extract_ids(neg, "neg", "VEGFsR1"),
    "neg_VEGFsR1.1" = extract_ids(neg, "neg", "VEGFsR1.1"),
    neg_combined     = extract_ids(neg, "neg")
  )

  # Report
  for (nm in names(lists)) {
    message(sprintf("  Gene list '%s': %d Entrez IDs", nm, length(lists[[nm]])))
  }

  lists
}


#' Get mapping rate statistics
#'
#' @param consensus_dir Character. Directory containing Step 1 consensus CSVs.
#' @return Named list with mapping_rate_pos, mapping_rate_neg, n_total, n_mapped
get_mapping_stats <- function(consensus_dir) {
  pos <- load_consensus(file.path(consensus_dir, "step01_consensus_proteins_pos.csv"))
  neg <- load_consensus(file.path(consensus_dir, "step01_consensus_proteins_neg.csv"))

  pos <- rescue_entrez_ids(pos)
  neg <- rescue_entrez_ids(neg)

  pos_unique <- pos[!duplicated(pos$target), ]
  neg_unique <- neg[!duplicated(neg$target), ]

  list(
    pos_total   = nrow(pos_unique),
    pos_mapped  = sum(!is.na(pos_unique$entrez_gene_id)),
    pos_rate    = if (nrow(pos_unique) > 0) sum(!is.na(pos_unique$entrez_gene_id)) / nrow(pos_unique) else NA,
    neg_total   = nrow(neg_unique),
    neg_mapped  = sum(!is.na(neg_unique$entrez_gene_id)),
    neg_rate    = if (nrow(neg_unique) > 0) sum(!is.na(neg_unique$entrez_gene_id)) / nrow(neg_unique) else NA
  )
}
