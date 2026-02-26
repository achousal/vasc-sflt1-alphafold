# 02_run_enrichment.R -- Run ORA enrichment across GO, KEGG, Reactome
#
# Over-Representation Analysis on discrete significant protein sets.
# Databases: GO:BP, GO:MF, GO:CC, KEGG, Reactome.
# Background: SomaScan v4.1 panel size (~7289) or user-provided.

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(readr)
})

# Default SomaScan v4.1 approximate panel size
SOMASCAN_V41_SIZE <- 7289L


#' Run GO enrichment (BP, MF, CC) for a gene list
#'
#' @param gene_ids Character vector of Entrez IDs.
#' @param universe Character vector of background Entrez IDs (or NULL for default).
#' @param padj_cutoff Numeric. Adjusted p-value cutoff.
#' @param qvalue_cutoff Numeric. q-value cutoff.
#' @param min_gs Integer. Minimum gene set size.
#' @param max_gs Integer. Maximum gene set size.
#' @return Named list with elements BP, MF, CC (each an enrichResult or NULL).
run_go_enrichment <- function(gene_ids, universe = NULL,
                              padj_cutoff = 0.05, qvalue_cutoff = 0.1,
                              min_gs = 10, max_gs = 500) {
  if (length(gene_ids) < 3) {
    message("  Skipping GO: fewer than 3 genes")
    return(list(BP = NULL, MF = NULL, CC = NULL))
  }

  ontologies <- c("BP", "MF", "CC")
  results <- list()

  for (ont in ontologies) {
    message(sprintf("  Running GO:%s ...", ont))
    res <- tryCatch({
      clusterProfiler::enrichGO(
        gene          = gene_ids,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        universe      = universe,
        pAdjustMethod = "BH",
        pvalueCutoff  = padj_cutoff,
        qvalueCutoff  = qvalue_cutoff,
        minGSSize     = min_gs,
        maxGSSize     = max_gs,
        readable      = TRUE
      )
    }, error = function(e) {
      message(sprintf("  GO:%s failed: %s", ont, conditionMessage(e)))
      NULL
    })
    results[[ont]] <- res
  }
  results
}


#' Run KEGG enrichment
#'
#' @param gene_ids Character vector of Entrez IDs.
#' @param universe Character vector of background Entrez IDs (or NULL).
#' @param padj_cutoff Numeric.
#' @param qvalue_cutoff Numeric.
#' @param min_gs Integer.
#' @param max_gs Integer.
#' @return enrichResult or NULL.
run_kegg_enrichment <- function(gene_ids, universe = NULL,
                                padj_cutoff = 0.05, qvalue_cutoff = 0.1,
                                min_gs = 10, max_gs = 500) {
  if (length(gene_ids) < 3) {
    message("  Skipping KEGG: fewer than 3 genes")
    return(NULL)
  }

  message("  Running KEGG ...")
  tryCatch({
    clusterProfiler::enrichKEGG(
      gene          = gene_ids,
      organism      = "hsa",
      universe      = universe,
      pAdjustMethod = "BH",
      pvalueCutoff  = padj_cutoff,
      qvalueCutoff  = qvalue_cutoff,
      minGSSize     = min_gs,
      maxGSSize     = max_gs
    )
  }, error = function(e) {
    message(sprintf("  KEGG failed: %s", conditionMessage(e)))
    NULL
  })
}


#' Run Reactome enrichment
#'
#' @param gene_ids Character vector of Entrez IDs.
#' @param universe Character vector of background Entrez IDs (or NULL).
#' @param padj_cutoff Numeric.
#' @param qvalue_cutoff Numeric.
#' @param min_gs Integer.
#' @param max_gs Integer.
#' @return enrichResult or NULL.
run_reactome_enrichment <- function(gene_ids, universe = NULL,
                                    padj_cutoff = 0.05, qvalue_cutoff = 0.1,
                                    min_gs = 10, max_gs = 500) {
  if (length(gene_ids) < 3) {
    message("  Skipping Reactome: fewer than 3 genes")
    return(NULL)
  }

  message("  Running Reactome ...")
  tryCatch({
    ReactomePA::enrichPathway(
      gene          = gene_ids,
      organism      = "human",
      universe      = universe,
      pAdjustMethod = "BH",
      pvalueCutoff  = padj_cutoff,
      qvalueCutoff  = qvalue_cutoff,
      minGSSize     = min_gs,
      maxGSSize     = max_gs,
      readable      = TRUE
    )
  }, error = function(e) {
    message(sprintf("  Reactome failed: %s", conditionMessage(e)))
    NULL
  })
}


#' Run all enrichment databases for a single gene list
#'
#' @param gene_ids Character vector of Entrez IDs.
#' @param list_name Character. Name of the gene list (for logging).
#' @param universe Character vector of background Entrez IDs (or NULL).
#' @return Named list: GO_BP, GO_MF, GO_CC, KEGG, Reactome (each enrichResult or NULL)
run_all_enrichments <- function(gene_ids, list_name, universe = NULL) {
  message(sprintf("\n== Enrichment for '%s' (%d genes) ==", list_name, length(gene_ids)))

  go <- run_go_enrichment(gene_ids, universe)
  kegg <- run_kegg_enrichment(gene_ids, universe)
  reactome <- run_reactome_enrichment(gene_ids, universe)

  list(
    GO_BP    = go$BP,
    GO_MF    = go$MF,
    GO_CC    = go$CC,
    KEGG     = kegg,
    Reactome = reactome
  )
}


#' Extract enrichment results as a tidy data.frame
#'
#' @param enrich_result An enrichResult object or NULL.
#' @param database Character. Database name tag.
#' @param list_name Character. Gene list name tag.
#' @return Tibble with enrichment results, or empty tibble if NULL/no results.
enrich_to_df <- function(enrich_result, database, list_name) {
  if (is.null(enrich_result)) {
    return(tibble::tibble())
  }
  df <- as.data.frame(enrich_result)
  if (nrow(df) == 0) return(tibble::tibble())
  df$database  <- database
  df$gene_list <- list_name
  tibble::as_tibble(df)
}


#' Save all enrichment results for a gene list
#'
#' @param results Named list from run_all_enrichments().
#' @param list_name Character. Gene list name (used in filenames).
#' @param results_dir Character. Output directory.
#' @return Tibble of all combined results.
save_enrichment_results <- function(results, list_name, results_dir) {
  all_dfs <- list()
  for (db_name in names(results)) {
    df <- enrich_to_df(results[[db_name]], db_name, list_name)
    if (nrow(df) > 0) {
      fname <- sprintf("step02_enrichment_%s_%s.csv", db_name, list_name)
      readr::write_csv(df, file.path(results_dir, fname))
      message(sprintf("  Saved: %s (%d terms)", fname, nrow(df)))
      all_dfs <- c(all_dfs, list(df))
    }
  }
  if (length(all_dfs) > 0) dplyr::bind_rows(all_dfs) else tibble::tibble()
}


#' Build a background universe from a file or use default panel size
#'
#' @param background_path Character or NULL. Path to CSV with Entrez IDs.
#' @return Character vector of Entrez IDs, or NULL if using default.
build_universe <- function(background_path = NULL) {
  if (is.null(background_path) || !file.exists(background_path)) {
    message(sprintf("  Using default SomaScan v4.1 universe (no explicit background)"))
    return(NULL)
  }

  message(sprintf("  Loading background from: %s", background_path))
  bg <- readr::read_csv(background_path, show_col_types = FALSE)

  # Try to find entrez_gene_id column
  id_col <- intersect(c("entrez_gene_id", "EntrezGeneID", "ENTREZID"), names(bg))
  if (length(id_col) == 0) {
    message("  WARNING: No Entrez ID column found in background file, using default")
    return(NULL)
  }

  ids <- unique(as.character(bg[[id_col[1]]]))
  ids <- ids[!is.na(ids) & ids != ""]
  message(sprintf("  Background universe: %d Entrez IDs", length(ids)))
  ids
}
