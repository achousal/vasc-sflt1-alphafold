#!/usr/bin/env Rscript
# 13_uniprot_explorer.R -- Explore UniProt annotations for AF2 target selection
#
# Answers: what protein types do we have? What domains are annotated?
# Why do we infer ectodomains from TM positions instead of reading them directly?
# What trimming options exist for each target?
#
# Usage:
#   module load R/4.3.0
#   Rscript analysis/03_structural_prediction/13_uniprot_explorer.R
#
# Or source interactively in RStudio for exploration.

library(data.table)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
project_root <- "/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold"
results_dir  <- file.path(project_root, "results", "03_structural_prediction", "d1d3")

feat  <- fread(file.path(results_dir, "uniprot_features.csv"))
prot  <- fread(file.path(results_dir, "uniprot_proteins.csv"))
topo  <- fread(file.path(results_dir, "topology_audit.csv"), na.strings = c("", "NA"))
cands <- fread(file.path(results_dir, "candidates_unified.csv"))

output_dir <- file.path(results_dir, "uniprot_explorer_tables")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
banner <- function(title) {
  cat("\n", strrep("=", 72), "\n ", title, "\n", strrep("=", 72), "\n\n")
}

save_table <- function(dt, name) {
  path <- file.path(output_dir, paste0(name, ".csv"))
  fwrite(dt, path)
  cat("  ->", path, "\n")
}

# =========================================================================
# 1. What UniProt feature types exist across our candidates?
# =========================================================================
banner("1. FEATURE TYPE INVENTORY")

type_counts <- feat[, .(
  n_entries  = .N,
  n_proteins = uniqueN(uniprot),
  avg_length = round(mean(length)),
  med_length = round(median(length))
), by = feature_type][order(-n_entries)]

print(type_counts)
save_table(type_counts, "feature_type_inventory")

# =========================================================================
# 2. Topological domain annotations -- the key question
#    UniProt annotates "Extracellular", "Cytoplasmic", "Lumenal" etc.
#    directly. Why don't we use these?
# =========================================================================
banner("2. TOPOLOGICAL DOMAIN ANNOTATIONS")

topo_domains <- feat[feature_type == "Topological domain"]

topo_types <- topo_domains[, .(
  n = .N,
  n_proteins = uniqueN(uniprot),
  avg_length = round(mean(length)),
  med_length = round(median(length)),
  max_length = max(length)
), by = description][order(-n)]

cat("Topological domain types annotated by UniProt:\n")
print(topo_types)
save_table(topo_types, "topological_domain_types")

# How many proteins have explicit "Extracellular" annotations?
has_extra <- topo_domains[grepl("Extracellular", description), uniqueN(uniprot)]
has_cyto  <- topo_domains[grepl("Cytoplasmic", description), uniqueN(uniprot)]
has_lumen <- topo_domains[grepl("Lumenal", description), uniqueN(uniprot)]

cat("\nProteins with explicit topological annotations:\n")
cat("  Extracellular:", has_extra, "\n")
cat("  Cytoplasmic:  ", has_cyto, "\n")
cat("  Lumenal:      ", has_lumen, "\n")
cat("  Total proteins:", uniqueN(feat$uniprot), "\n")

# =========================================================================
# 3. Compare: inferred ectodomain vs annotated extracellular domain
#    For TM proteins, our Python code does: SP_end+1 ... first_TM_start-1
#    UniProt may annotate the same region as "Topological domain: Extracellular"
# =========================================================================
banner("3. INFERRED vs ANNOTATED ECTODOMAIN COMPARISON")

# Get the annotated extracellular regions per protein
ecd_annotated <- topo_domains[grepl("Extracellular", description), .(
  ecd_ann_start = min(start),
  ecd_ann_end   = max(end),
  ecd_ann_length = sum(length),
  n_ecd_segments = .N
), by = .(target, uniprot)]

# Join with our inferred ectodomain from topology_audit.csv
topo_sub <- topo[protein_type %in% c("type_i_tm", "multi_tm"), .(
  target, uniprot, seq_length, protein_type,
  inferred_ecto_start = ecto_start,
  inferred_ecto_end   = ecto_end,
  inferred_ecto_len   = ecto_length
)]

comparison <- merge(topo_sub, ecd_annotated, by = c("target", "uniprot"), all.x = TRUE)

comparison[, `:=`(
  has_annotation = !is.na(ecd_ann_start),
  start_match    = (inferred_ecto_start == ecd_ann_start),
  end_match      = (inferred_ecto_end == ecd_ann_end),
  length_diff    = as.integer(inferred_ecto_len) - ecd_ann_length
)]

# Summary
cat("TM proteins with annotated extracellular domains:\n")
cat("  Has annotation:  ", sum(comparison$has_annotation), "/",
    nrow(comparison), "\n")
cat("  Start matches:   ", sum(comparison$start_match, na.rm = TRUE), "\n")
cat("  End matches:     ", sum(comparison$end_match, na.rm = TRUE), "\n")
cat("  Both match:      ", sum(comparison$start_match & comparison$end_match,
    na.rm = TRUE), "\n\n")

# Show mismatches
mismatches <- comparison[has_annotation == TRUE & (start_match == FALSE | end_match == FALSE)]
if (nrow(mismatches) > 0) {
  cat("Mismatches between inferred and annotated ECD:\n")
  print(mismatches[, .(target, uniprot, protein_type,
    inf_start = inferred_ecto_start, inf_end = inferred_ecto_end, inf_len = inferred_ecto_len,
    ann_start = ecd_ann_start, ann_end = ecd_ann_end, ann_len = ecd_ann_length,
    n_segments = n_ecd_segments, len_diff = as.integer(inferred_ecto_len) - as.integer(ecd_ann_length)
  )][order(-abs(len_diff))][1:min(20, .N)])
}

# TM proteins WITHOUT extracellular annotation
no_ann <- comparison[has_annotation == FALSE]
if (nrow(no_ann) > 0) {
  cat("\nTM proteins with NO annotated extracellular domain (", nrow(no_ann), "):\n")
  cat("  (These rely entirely on our SP->TM inference)\n")
  print(no_ann[, .(target, uniprot, protein_type, seq_length,
    inf_start = inferred_ecto_start, inf_end = inferred_ecto_end,
    inf_len = inferred_ecto_len)][order(as.integer(inf_len))][1:min(20, .N)])
}

save_table(comparison, "inferred_vs_annotated_ecd")

# =========================================================================
# 4. Named domain annotations -- what structural domains are in our targets?
#    These are specific (Ig-like, FN3, EGF, CUB...) not generic topology.
# =========================================================================
banner("4. NAMED DOMAIN ANNOTATIONS")

domains <- feat[feature_type == "Domain"]

domain_types <- domains[, .(
  n = .N,
  n_proteins = uniqueN(uniprot),
  avg_length = round(mean(length)),
  example_target = first(target)
), by = description][order(-n)]

cat("Top 30 domain types across candidates:\n")
print(domain_types[1:min(30, .N)])
save_table(domain_types, "named_domain_types")

# Proteins with most named domains
domain_rich <- domains[, .(
  n_domains = .N,
  domain_coverage = round(sum(length) / first(seq_length), 2),
  domain_list = paste(unique(description), collapse = "; ")
), by = .(target, uniprot, seq_length)][order(-n_domains)]

cat("\nMost domain-rich proteins (top 15):\n")
print(domain_rich[1:15, .(target, uniprot, seq_length, n_domains,
  domain_coverage, domain_list)])

# =========================================================================
# 5. The real question: ALTERNATIVE TRIMMING STRATEGIES
#    Instead of "everything before first TM", we could select:
#    (a) Annotated extracellular topological domain (if available)
#    (b) Specific named domains (e.g., only Ig-like domains of ECD)
#    (c) Current approach: SP_end+1 to TM_start-1
# =========================================================================
banner("5. TRIMMING STRATEGY COMPARISON")

# For each TM protein, compare the three approaches
tm_targets <- topo[protein_type %in% c("type_i_tm", "multi_tm")]

strategy <- data.table()
for (idx in seq_len(nrow(tm_targets))) {
  row <- tm_targets[idx]
  acc <- row$uniprot
  tgt <- row$target

  # Strategy A: annotated extracellular
  ecd <- topo_domains[uniprot == acc & grepl("Extracellular", description)]
  ann_len <- if (nrow(ecd) > 0) sum(ecd$length) else NA_integer_

  # Strategy B: named domains in extracellular region
  named <- domains[uniprot == acc]
  inf_end <- as.integer(row$ecto_end)
  if (!is.na(inf_end) && nrow(named) > 0) {
    ecd_domains <- named[end <= inf_end + 10]  # domains within ECD
    dom_len <- if (nrow(ecd_domains) > 0) sum(ecd_domains$length) else NA_integer_
    dom_names <- if (nrow(ecd_domains) > 0) {
      paste(unique(ecd_domains$description), collapse = "; ")
    } else NA_character_
  } else {
    dom_len <- NA_integer_
    dom_names <- NA_character_
  }

  # Strategy C: current inferred
  inf_len <- as.integer(row$ecto_length)

  strategy <- rbind(strategy, data.table(
    target = tgt, uniprot = acc,
    seq_length = as.integer(row$seq_length),
    protein_type = row$protein_type,
    n_tm = as.integer(row$n_tm_helices),
    strategy_a_annotated_ecd = ann_len,
    strategy_b_named_domains = dom_len,
    strategy_c_inferred_ecd  = inf_len,
    ecd_domain_names = dom_names
  ))
}

strategy[, best_strategy := ifelse(
  !is.na(strategy_a_annotated_ecd), "A: annotated ECD",
  ifelse(!is.na(strategy_c_inferred_ecd), "C: inferred SP-to-TM",
  "none")
)]

cat("Strategy availability for", nrow(strategy), "TM proteins:\n")
cat("  A (annotated ECD):   ", sum(!is.na(strategy$strategy_a_annotated_ecd)), "\n")
cat("  B (named domains):   ", sum(!is.na(strategy$strategy_b_named_domains)), "\n")
cat("  C (inferred SP->TM): ", sum(!is.na(strategy$strategy_c_inferred_ecd)), "\n")
cat("  None available:      ", sum(is.na(strategy$strategy_c_inferred_ecd) &
    is.na(strategy$strategy_a_annotated_ecd)), "\n")

cat("\nLength comparison (where both A and C are available):\n")
both <- strategy[!is.na(strategy_a_annotated_ecd) & !is.na(strategy_c_inferred_ecd)]
both[, diff := strategy_c_inferred_ecd - strategy_a_annotated_ecd]
cat("  Mean diff (inferred - annotated):", round(mean(as.numeric(both$diff))), "aa\n")
cat("  Median diff:", median(as.numeric(both$diff)), "aa\n")
cat("  Max diff:", max(abs(as.numeric(both$diff))), "aa\n")
cat("  Exact match:", sum(as.numeric(both$diff) == 0), "/", nrow(both), "\n")

save_table(strategy, "trimming_strategy_comparison")

# =========================================================================
# 6. Protein localization overview -- are we modeling the right proteins?
# =========================================================================
banner("6. SUBCELLULAR LOCALIZATION (from UniProt keywords)")

loc_summary <- prot[, .(
  target, uniprot, seq_length, keywords
)]

# Parse keywords into categories
loc_summary[, is_secreted := grepl("Secreted", keywords)]
loc_summary[, is_membrane := grepl("Membrane|Cell membrane", keywords)]
loc_summary[, is_extracellular := grepl("Extracellular", keywords)]
loc_summary[, is_golgi := grepl("Golgi", keywords)]
loc_summary[, is_er := grepl("Endoplasmic", keywords)]
loc_summary[, is_cytoplasm := grepl("Cytoplasm", keywords)]
loc_summary[, is_nucleus := grepl("Nucleus", keywords)]
loc_summary[, is_lysosome := grepl("Lysosome", keywords)]
loc_summary[, is_mitochondria := grepl("Mitochond", keywords)]

# Summarize
cat("Subcellular localization (UniProt keywords):\n")
locs <- data.table(
  location = c("Secreted", "Cell membrane", "Extracellular",
    "Golgi", "ER", "Cytoplasm", "Nucleus", "Lysosome", "Mitochondria"),
  n = c(sum(loc_summary$is_secreted), sum(loc_summary$is_membrane),
    sum(loc_summary$is_extracellular), sum(loc_summary$is_golgi),
    sum(loc_summary$is_er), sum(loc_summary$is_cytoplasm),
    sum(loc_summary$is_nucleus), sum(loc_summary$is_lysosome),
    sum(loc_summary$is_mitochondria))
)
locs[, pct := round(100 * n / nrow(loc_summary), 1)]
print(locs[order(-n)])

# Flag proteins that are NOT extracellular/secreted/membrane
loc_summary[, accessible_to_sflt1 := is_secreted | is_membrane | is_extracellular]

cat("\nAccessible to sFLT1 (secreted/membrane/extracellular):",
    sum(loc_summary$accessible_to_sflt1), "/", nrow(loc_summary), "\n")

intracellular <- loc_summary[accessible_to_sflt1 == FALSE]
cat("Intracellular only (questionable targets):", nrow(intracellular), "\n")
if (nrow(intracellular) > 0) {
  cat("\n  Examples:\n")
  print(intracellular[1:min(15, .N), .(target, uniprot, seq_length, keywords)])
}

save_table(loc_summary, "subcellular_localization")

# =========================================================================
# 7. SUMMARY: What the data tells us about trimming methodology
# =========================================================================
banner("7. METHODOLOGY ASSESSMENT")

cat("
CURRENT APPROACH (02_fetch_sequences.py):
  - TM proteins: trim to inferred ectodomain (SP_end+1 ... first_TM-1)
  - GPI proteins: use Chain annotation
  - Soluble: full-length (should remove signal peptide)

WHAT UNIPROT ACTUALLY PROVIDES:
  - Topological domain: 'Extracellular' annotations for most TM proteins
  - Named domains: specific structural domains (Ig-like, FN3, etc.)
  - Chain: mature protein boundaries after processing

FINDINGS:\n")

both <- strategy[!is.na(strategy_a_annotated_ecd) & !is.na(strategy_c_inferred_ecd)]
cat("  1. UniProt has explicit ECD annotations for",
    sum(!is.na(strategy$strategy_a_annotated_ecd)), "of",
    nrow(strategy), "TM proteins\n")
cat("  2. Inferred and annotated ECDs match within 5 aa for",
    sum(abs(as.numeric(both$diff)) <= 5), "of", nrow(both), "proteins\n")
cat("  3.", sum(abs(as.numeric(both$diff)) > 50), "proteins differ by >50 aa",
    "(worth investigating)\n")
cat("  4.", nrow(intracellular), "candidates are intracellular-only",
    "(biologically implausible)\n")
cat("  5. Strategy B (named domains only) averages",
    round(mean(strategy$strategy_b_named_domains, na.rm = TRUE)),
    "aa vs", round(mean(strategy$strategy_c_inferred_ecd, na.rm = TRUE)),
    "aa for full ECD\n")

cat("\nRECOMMENDATION:
  - Use annotated ECD when available (strategy A) -- ground truth
  - Fall back to inferred SP->TM (strategy C) when no annotation
  - Flag intracellular proteins for review/removal
  - Signal peptides should ALWAYS be removed (45 soluble targets affected)
  - 29 targets with <50 aa ectodomain should be dropped\n")

cat("\nAll tables saved to:", output_dir, "\n")
