#!/usr/bin/env Rscript
# 14_uniprot_api.R -- UniProt explorer for AF2 target grounding
#
# Two modes:
#   1. Fetch any accession(s), see what UniProt knows.
#      What feature types exist? What topology? What domains?
#   2. GROUNDING: fetch all our targets, see how they map onto
#      UniProt categories and what trimming decisions follow.
#
# Usage (interactive -- source then explore):
#   source("analysis/03_structural_prediction/14_uniprot_api.R")
#
#   # --- Explore one protein ---
#   explore("P17948")                   # sFLT1 itself
#   explore("O14786")                   # NRP1
#   explore("P21860")                   # ERBB3
#
#   # --- What does UniProt give us? ---
#   x <- fetch_one("O14786")
#   names(x)                            # raw JSON fields
#   map_features(x)                     # all feature types + counts
#   map_keywords(x)                     # all keywords
#   map_comments(x)                     # comment types (function, subcell, etc.)
#
#   # --- Ground our targets ---
#   fetch_candidates("results/03_structural_prediction/step03_candidates.csv")
#   View(prot)       # protein-level decisions
#   View(feat)       # all features flat
#   View(domains)    # named structural domains
#   landscape()      # what categories exist across all targets
#   audit()          # trimming decisions + flags
#
# Usage (batch):
#   Rscript 14_uniprot_api.R <candidates.csv> [output_prefix]

library(data.table)
library(httr)
library(jsonlite)


# =========================================================================
# PART 1: RAW API -- see what UniProt gives you
# =========================================================================

fetch_one <- function(accession) {
  #' Fetch raw UniProt JSON. Returns the parsed list -- poke around.
  url <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".json")
  resp <- GET(url, timeout(30))
  if (status_code(resp) == 429) { Sys.sleep(5); resp <- GET(url, timeout(30)) }
  if (status_code(resp) != 200) { warning("HTTP ", status_code(resp)); return(NULL) }
  content(resp, as = "parsed", simplifyVector = FALSE)
}


map_features <- function(data) {
  #' What feature types does this protein have? Counts + total residues.
  feats <- data$features %||% list()
  rows <- lapply(feats, function(f) {
    s <- f$location$start$value; e <- f$location$end$value
    if (is.null(s) || is.null(e)) return(NULL)
    data.table(type = f$type %||% "?", len = as.integer(e - s + 1),
               desc = substr(f$description %||% "", 1, 60))
  })
  dt <- rbindlist(Filter(Negate(is.null), rows))
  if (nrow(dt) == 0) return(data.table())
  dt[, .(n = .N, total_aa = sum(len), examples = paste(head(unique(desc), 3), collapse = " | ")),
     by = type][order(-n)]
}


map_keywords <- function(data) {
  #' All UniProt keywords for this protein.
  kws <- data$keywords %||% list()
  data.table(
    category = vapply(kws, function(k) k$category %||% "", ""),
    keyword  = vapply(kws, function(k) k$name %||% "", "")
  )[order(category, keyword)]
}


map_comments <- function(data) {
  #' What comment types (function, subcellular location, tissue specificity...)?
  cmts <- data$comments %||% list()
  rows <- lapply(cmts, function(c) {
    ctype <- c$commentType %||% "?"
    # Extract text from nested structures
    texts <- c$texts %||% list()
    txt <- if (length(texts) > 0) {
      paste(vapply(texts, function(t) substr(t$value %||% "", 1, 120), ""), collapse = " ")
    } else ""
    # Subcellular location has a different structure
    if (ctype == "SUBCELLULAR LOCATION") {
      locs <- c$subcellularLocations %||% list()
      loc_names <- vapply(locs, function(l) l$location$value %||% "", "")
      txt <- paste(loc_names, collapse = "; ")
    }
    data.table(type = ctype, text = txt)
  })
  rbindlist(rows)
}


explore <- function(accession) {
  #' One-shot: fetch + print everything useful about a protein.
  cat("Fetching", accession, "...\n")
  data <- fetch_one(accession)
  if (is.null(data)) { cat("Failed.\n"); return(invisible()) }

  # Identity
  name <- data$proteinDescription$recommendedName$fullName$value %||%
    data$proteinDescription$submittedName[[1]]$fullName$value %||% "?"
  gene <- if (length(data$genes) > 0) data$genes[[1]]$geneName$value %||% "?" else "?"
  org  <- data$organism$scientificName %||% "?"
  slen <- data$sequence$length %||% 0

  cat("\n", strrep("-", 60), "\n")
  cat(" ", name, "\n")
  cat("  Gene:", gene, " | Accession:", accession, " | Length:", slen, "aa\n")
  cat("  Organism:", org, "\n")
  cat(strrep("-", 60), "\n")

  # Keywords by category
  kw <- map_keywords(data)
  if (nrow(kw) > 0) {
    cat("\nKeywords:\n")
    for (cat_name in unique(kw$category)) {
      cat("  ", cat_name, ": ", paste(kw[category == cat_name, keyword], collapse = ", "), "\n")
    }
  }

  # Subcellular location (from comments)
  cmts <- map_comments(data)
  subcell <- cmts[type == "SUBCELLULAR LOCATION"]
  if (nrow(subcell) > 0) {
    cat("\nSubcellular location:\n  ", subcell$text[1], "\n")
  }

  # Feature map
  cat("\nFeature types:\n")
  fm <- map_features(data)
  if (nrow(fm) > 0) print(fm)

  # Topology detail
  feats <- data$features %||% list()
  topo_types <- c("Signal", "Transmembrane", "Topological domain",
                  "Chain", "Lipidation", "Intramembrane", "Propeptide")
  topo <- lapply(feats, function(f) {
    if (!(f$type %||% "") %in% topo_types) return(NULL)
    s <- f$location$start$value; e <- f$location$end$value
    if (is.null(s) || is.null(e)) return(NULL)
    data.table(type = f$type, start = s, end = e, len = e - s + 1,
               desc = f$description %||% "")
  })
  topo <- rbindlist(Filter(Negate(is.null), topo))
  if (nrow(topo) > 0) {
    cat("\nTopology:\n")
    print(topo)
  }

  # Named domains
  doms <- lapply(feats, function(f) {
    if (!identical(f$type, "Domain")) return(NULL)
    s <- f$location$start$value; e <- f$location$end$value
    if (is.null(s) || is.null(e)) return(NULL)
    data.table(domain = f$description %||% "?", start = s, end = e, len = e - s + 1)
  })
  doms <- rbindlist(Filter(Negate(is.null), doms))
  if (nrow(doms) > 0) {
    cat("\nNamed domains:\n")
    print(doms)
    cat("  Domain coverage:", round(sum(doms$len) / slen * 100, 1), "% of sequence\n")
  }

  # AF2 modeling recommendation
  cat("\nAF2 recommendation:\n")
  has_tm <- nrow(topo[type == "Transmembrane"]) > 0
  has_sp <- nrow(topo[type == "Signal"]) > 0
  has_ecd <- nrow(topo[type == "Topological domain" & grepl("Extracellular", desc)]) > 0

  if (has_tm && has_ecd) {
    ecd_rows <- topo[type == "Topological domain" & grepl("Extracellular", desc)]
    ecd_range <- paste0(min(ecd_rows$start), "-", max(ecd_rows$end))
    ecd_len <- sum(ecd_rows$len)
    cat("  TM protein with annotated ECD:", ecd_range, "(", ecd_len, "aa )\n")
    if (ecd_len < 50) cat("  WARNING: ECD < 50 aa -- skip for AF2\n")
  } else if (has_tm) {
    sp_end <- if (has_sp) max(topo[type == "Signal", end]) else 0
    tm_start <- min(topo[type == "Transmembrane", start])
    inf_len <- tm_start - sp_end - 2
    cat("  TM protein, NO annotated ECD. Inferred:", sp_end + 1, "-", tm_start - 1,
        "(", inf_len, "aa )\n")
    if (inf_len < 50) cat("  WARNING: inferred ECD < 50 aa -- skip for AF2\n")
  } else if (has_sp) {
    sp_end <- max(topo[type == "Signal", end])
    cat("  Soluble/secreted. Remove signal peptide (1-", sp_end, "), model",
        sp_end + 1, "-", slen, "(", slen - sp_end, "aa )\n")
  } else {
    cat("  Soluble, no signal peptide. Model full-length (", slen, "aa )\n")
  }

  invisible(data)
}


# =========================================================================
# PART 2: BATCH FETCH + GROUNDING
# =========================================================================

parse_features_flat <- function(data, accession, target = NA_character_) {
  #' Parse all features into a flat table.
  seq_len <- data$sequence$length %||% 0L
  feats <- data$features %||% list()
  rows <- lapply(feats, function(f) {
    s <- f$location$start$value; e <- f$location$end$value
    if (is.null(s) || is.null(e)) return(NULL)
    data.table(target = target, uniprot = accession, seq_length = as.integer(seq_len),
               feature_type = f$type %||% "", start = as.integer(s), end = as.integer(e),
               length = as.integer(e - s + 1), description = f$description %||% "")
  })
  rbindlist(Filter(Negate(is.null), rows))
}


parse_protein_row <- function(data, accession, target = NA_character_) {
  #' Parse one protein into a summary row with trimming decision.
  seq_len <- data$sequence$length %||% 0L
  feats <- data$features %||% list()

  kws <- vapply(data$keywords %||% list(), function(k) k$name %||% "", "")
  is_gpi <- "GPI-anchor" %in% kws

  # Subcellular keywords
  loc_terms <- c("Cell membrane", "Secreted", "Cell surface",
                 "Extracellular space", "Extracellular matrix",
                 "Cytoplasm", "Nucleus", "Mitochondrion",
                 "Endoplasmic reticulum", "Golgi apparatus", "Lysosome", "Peroxisome")
  accessible_kws <- c("Cell membrane", "Secreted", "Cell surface",
                       "Extracellular space", "Extracellular matrix")
  intracellular_kws <- c("Cytoplasm", "Nucleus", "Mitochondrion",
                          "Endoplasmic reticulum", "Golgi apparatus", "Lysosome", "Peroxisome")
  subcell <- kws[kws %in% loc_terms]

  # Parse features by type
  sp <- tm <- ecd <- chain <- NULL
  for (f in feats) {
    s <- f$location$start$value; e <- f$location$end$value
    if (is.null(s) || is.null(e)) next
    if (f$type == "Signal" && is.null(sp)) sp <- c(s, e)
    if (f$type == "Transmembrane") tm <- c(tm, list(c(s, e)))
    if (f$type == "Topological domain" && grepl("Extracellular", f$description %||% ""))
      ecd <- c(ecd, list(c(s, e)))
    if (f$type == "Chain" && is.null(chain)) chain <- c(s, e)
  }

  n_tm <- length(tm)
  protein_type <- if (n_tm > 1) "multi_tm"
    else if (n_tm == 1) "type_i_tm"
    else if (is_gpi) "gpi_anchored"
    else "soluble"

  sp_end <- if (!is.null(sp)) sp[2] else NA_integer_

  # Annotated ECD
  ecd_ann_start <- ecd_ann_end <- ecd_ann_len <- NA_integer_
  if (length(ecd) > 0) {
    ecd_ann_start <- min(vapply(ecd, `[`, 0L, 1))
    ecd_ann_end   <- max(vapply(ecd, `[`, 0L, 2))
    ecd_ann_len   <- sum(vapply(ecd, function(x) x[2] - x[1] + 1L, 0L))
  }

  # Inferred ECD
  inf_start <- inf_end <- inf_len <- NA_integer_
  if (n_tm > 0) {
    inf_start <- if (!is.na(sp_end)) sp_end + 1L else 1L
    inf_end   <- tm[[1]][1] - 1L
    if (inf_end >= inf_start) inf_len <- inf_end - inf_start + 1L
  }

  # Best ECD
  best_start <- if (!is.na(ecd_ann_start)) ecd_ann_start else inf_start
  best_end   <- if (!is.na(ecd_ann_end))   ecd_ann_end   else inf_end
  best_len   <- if (!is.na(ecd_ann_len))   ecd_ann_len   else inf_len
  best_src   <- if (!is.na(ecd_ann_start)) "annotated"
    else if (!is.na(inf_start)) "inferred" else NA_character_

  # Mature range
  mature_start <- if (!is.na(sp_end)) sp_end + 1L
    else if (!is.null(chain)) chain[1] else 1L
  mature_end <- if (!is.null(chain)) chain[2] else as.integer(seq_len)

  # Accessibility
  is_accessible <- any(subcell %in% accessible_kws) ||
    protein_type %in% c("type_i_tm", "multi_tm", "gpi_anchored")
  is_intra <- length(subcell) > 0 && all(subcell %in% intracellular_kws) && !is_accessible

  # Model decision
  model_region <- if (is_intra) "SKIP:intracellular"
    else if (protein_type %in% c("type_i_tm", "multi_tm") && !is.na(best_len) && best_len >= 50)
      paste0("ecd_", best_start, "-", best_end)
    else if (protein_type %in% c("type_i_tm", "multi_tm") && (is.na(best_len) || best_len < 50))
      "SKIP:ecd_too_short"
    else if (protein_type == "gpi_anchored")
      paste0("mature_", mature_start, "-", mature_end)
    else if (!is.na(sp_end))
      paste0("mature_", mature_start, "-", mature_end)
    else "full_length"

  model_length <- if (grepl("^SKIP", model_region)) 0L
    else if (grepl("^ecd_", model_region)) as.integer(best_len)
    else as.integer(mature_end - mature_start + 1L)

  # Named domain count
  n_domains <- sum(vapply(feats, function(f) identical(f$type, "Domain"), FALSE))

  data.table(
    target = target, uniprot = accession, seq_length = as.integer(seq_len),
    protein_type = protein_type, keywords = paste(subcell, collapse = "; "),
    is_accessible = is_accessible, is_intracellular = is_intra,
    sp_end = sp_end, n_tm = n_tm,
    ecd_ann_len = ecd_ann_len, inf_len = inf_len,
    best_ecd_start = best_start, best_ecd_end = best_end,
    best_ecd_len = best_len, best_ecd_source = best_src,
    n_domains = n_domains,
    model_region = model_region, model_length = model_length
  )
}


fetch_candidates <- function(csv_path, uniprot_col = "uniprot",
                             target_col = "target", delay = 0.4) {
  #' Fetch UniProt topology for all proteins in a CSV.
  #' Populates global: prot, feat, domains
  cands <- fread(csv_path)
  if (!uniprot_col %in% names(cands)) stop("No '", uniprot_col, "' column")
  cands <- cands[!duplicated(get(uniprot_col))]
  n <- nrow(cands)
  cat("Fetching", n, "proteins from UniProt...\n")

  prot_list <- feat_list <- list()
  for (i in seq_len(n)) {
    acc <- cands[[uniprot_col]][i]
    tgt <- if (target_col %in% names(cands)) cands[[target_col]][i] else acc
    data <- fetch_one(acc)
    if (is.null(data)) { cat("  FAILED:", tgt, "\n"); next }
    prot_list[[i]] <- parse_protein_row(data, acc, tgt)
    feat_list[[i]] <- parse_features_flat(data, acc, tgt)
    if (i %% 25 == 0) cat("  ", i, "/", n, "\n")
    Sys.sleep(delay)
  }

  proteins <- rbindlist(Filter(Negate(is.null), prot_list))
  features <- rbindlist(Filter(Negate(is.null), feat_list))

  assign("prot", proteins, envir = .GlobalEnv)
  assign("feat", features, envir = .GlobalEnv)
  assign("domains", features[feature_type == "Domain"], envir = .GlobalEnv)

  cat("\nDone.", nrow(proteins), "proteins,", nrow(features), "features.\n")
  landscape()
  invisible(proteins)
}


# =========================================================================
# PART 3: LANDSCAPE + AUDIT -- see where our targets land
# =========================================================================

landscape <- function() {
  #' What UniProt categories exist across our fetched targets?
  if (!exists("prot", envir = .GlobalEnv)) stop("Run fetch_candidates() first")
  p <- get("prot", envir = .GlobalEnv)
  f <- get("feat", envir = .GlobalEnv)

  cat("\n", strrep("=", 60), "\n LANDSCAPE: what we're working with\n",
      strrep("=", 60), "\n")

  # Protein types
  cat("\nProtein types:\n")
  print(p[, .(n = .N, pct = round(.N / nrow(p) * 100, 1),
    avg_len = round(mean(seq_length)),
    skip = sum(grepl("SKIP", model_region))), by = protein_type][order(-n)])

  # Feature types across all proteins
  cat("\nFeature type coverage (how many proteins have each):\n")
  cov <- f[, .(n_proteins = uniqueN(uniprot), total_entries = .N,
    avg_per_protein = round(.N / uniqueN(uniprot), 1)),
    by = feature_type][order(-n_proteins)]
  print(cov[1:min(20, nrow(cov))])

  # Subcellular keywords
  cat("\nSubcellular locations:\n")
  kw_list <- strsplit(p$keywords, "; ")
  all_kw <- unlist(kw_list)
  all_kw <- all_kw[all_kw != ""]
  if (length(all_kw) > 0) {
    kw_dt <- data.table(kw = all_kw)[, .N, by = kw][order(-N)]
    kw_dt[, pct := round(N / nrow(p) * 100, 1)]
    print(kw_dt)
  }

  # Annotation completeness
  cat("\nAnnotation completeness:\n")
  cat("  Has signal peptide:     ", sum(!is.na(p$sp_end)), "/", nrow(p), "\n")
  cat("  Has TM helix:           ", sum(p$n_tm > 0), "/", nrow(p), "\n")
  cat("  Has annotated ECD:      ", sum(!is.na(p$ecd_ann_len)), "/", nrow(p), "\n")
  cat("  Has named domains:      ", sum(p$n_domains > 0), "/", nrow(p), "\n")
  cat("  Accessible to sFLT1:    ", sum(p$is_accessible), "/", nrow(p), "\n")
  cat("  Intracellular only:     ", sum(p$is_intracellular), "/", nrow(p), "\n")

  # Domain families
  d <- get("domains", envir = .GlobalEnv)
  if (nrow(d) > 0) {
    cat("\nTop domain families:\n")
    # Normalize: "Ig-like C2-type 1" -> "Ig-like C2-type"
    d[, family := gsub(" [0-9]+$", "", description)]
    fam <- d[, .(n_domains = .N, n_proteins = uniqueN(uniprot)), by = family][order(-n_proteins)]
    print(fam[1:min(15, nrow(fam))])
  }
}


audit <- function() {
  #' Trimming decisions + flags for review.
  if (!exists("prot", envir = .GlobalEnv)) stop("Run fetch_candidates() first")
  p <- get("prot", envir = .GlobalEnv)

  cat("\n", strrep("=", 60), "\n AUDIT: model region decisions\n",
      strrep("=", 60), "\n")

  # Decision summary
  cat("\nDecisions:\n")
  p[, decision := fifelse(grepl("^SKIP", model_region), model_region,
    fifelse(grepl("^ecd_", model_region), "ecd",
    fifelse(grepl("^mature_", model_region), "mature", "full_length")))]
  print(p[, .N, by = decision][order(-N)])

  # ECD source for TM proteins
  tm <- p[protein_type %in% c("type_i_tm", "multi_tm") & !grepl("SKIP", model_region)]
  if (nrow(tm) > 0) {
    cat("\nECD source (TM proteins):\n")
    print(tm[, .N, by = best_ecd_source])

    # Mismatches
    mm <- tm[!is.na(ecd_ann_len) & !is.na(inf_len) & abs(ecd_ann_len - inf_len) > 20]
    if (nrow(mm) > 0) {
      cat("\nAnnotated vs inferred ECD mismatches (>20 aa):\n")
      print(mm[, .(target, uniprot, protein_type,
        annotated = ecd_ann_len, inferred = inf_len,
        diff = ecd_ann_len - inf_len, used = best_ecd_source)][order(-abs(diff))])
    }
  }

  # Skips
  skips <- p[grepl("SKIP", model_region)]
  if (nrow(skips) > 0) {
    cat("\nSkipped targets (", nrow(skips), "):\n")
    print(skips[, .(target, uniprot, protein_type, seq_length, keywords, model_region)])
  }

  # Size distribution for modeled proteins
  modeled <- p[!grepl("SKIP", model_region)]
  cat("\nModel length distribution (", nrow(modeled), " proteins ):\n")
  cat("  Min:    ", min(modeled$model_length), "aa\n")
  cat("  Median: ", median(modeled$model_length), "aa\n")
  cat("  Mean:   ", round(mean(modeled$model_length)), "aa\n")
  cat("  Max:    ", max(modeled$model_length), "aa\n")

  sflt1 <- 304L
  modeled[, total := model_length + sflt1]
  modeled[, walltime := fifelse(total < 700, "48h",
    fifelse(total < 1200, "72h", fifelse(total < 1600, "96h", "144h")))]
  cat("\nWalltime tiers (with sFLT1 D1-D3 = 304 aa):\n")
  print(modeled[, .N, by = walltime][order(walltime)])
}


# =========================================================================
# INTERACTIVE HELPERS
# =========================================================================

show_protein <- function(x) {
  #' Full detail for one protein by accession or target name.
  p <- get("prot", envir = .GlobalEnv)
  row <- p[uniprot == x | target == x]
  if (nrow(row) == 0) { cat("Not found:", x, "\n"); return(invisible()) }
  t(row)
}

show_domains <- function(x) {
  #' Named domains for one protein.
  f <- get("feat", envir = .GlobalEnv)
  f[feature_type == "Domain" & (uniprot == x | target == x),
    .(description, start, end, length)]
}

show_topology <- function(x) {
  #' Topology features for one protein.
  f <- get("feat", envir = .GlobalEnv)
  topo_types <- c("Signal", "Transmembrane", "Topological domain", "Chain",
                  "Lipidation", "Intramembrane", "Transit peptide", "Propeptide")
  f[feature_type %in% topo_types & (uniprot == x | target == x),
    .(feature_type, start, end, length, description)]
}


# =========================================================================
# CLI
# =========================================================================
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat("Usage: Rscript 14_uniprot_api.R <candidates.csv> [output_prefix]\n")
    quit(status = 1)
  }
  csv_path <- args[1]
  prefix <- if (length(args) >= 2) args[2] else sub("\\.csv$", "", csv_path)

  fetch_candidates(csv_path)
  audit()

  fwrite(get("prot", envir = .GlobalEnv), paste0(prefix, "_topology.csv"))
  fwrite(get("feat", envir = .GlobalEnv), paste0(prefix, "_features.csv"))
  cat("\nSaved:", paste0(prefix, "_topology.csv"), "\n")
  cat("Saved:", paste0(prefix, "_features.csv"), "\n")
}
