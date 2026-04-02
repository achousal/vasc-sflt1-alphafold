#!/usr/bin/env Rscript
# 12_topology_audit.R -- Diagnostic tables for AF2 target topology and trimming
#
# Reads topology_audit.csv (from UniProt API fetch) and generates readable
# summary tables for reviewing target selection methodology.
#
# Usage:
#   module load R/4.3.0
#   Rscript analysis/03_structural_prediction/12_topology_audit.R

library(data.table)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
project_root <- "/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold"
results_dir  <- file.path(project_root, "results", "03_structural_prediction", "d1d3")
audit_path   <- file.path(results_dir, "topology_audit.csv")
output_dir   <- file.path(results_dir, "topology_audit_tables")

if (!file.exists(audit_path)) {
  stop("topology_audit.csv not found at: ", audit_path)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
dt <- fread(audit_path, na.strings = c("", "NA"))

sflt1_length <- 304L

dt[, total_current := current_target_residues + sflt1_length]
dt[, total_corrected := ifelse(
  needs_correction == TRUE & !is.na(suggested_length),
  as.integer(suggested_length) + sflt1_length,
  total_current
)]

walltime_tier <- function(total_res) {
  ifelse(total_res < 700, "48h",
  ifelse(total_res < 1200, "72h",
  ifelse(total_res < 1600, "96h", "144h")))
}

dt[, walltime_current   := walltime_tier(total_current)]
dt[, walltime_corrected := walltime_tier(total_corrected)]
dt[, walltime_changed   := walltime_current != walltime_corrected]

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
banner <- function(title, con = stdout()) {
  cat("\n", strrep("=", 70), "\n", file = con)
  cat(" ", title, "\n", file = con)
  cat(strrep("=", 70), "\n\n", file = con)
}

# ---------------------------------------------------------------------------
# Table 1: Protein Type Breakdown
# ---------------------------------------------------------------------------
type_summary <- dt[, .(
  n          = .N,
  need_fix   = sum(needs_correction == TRUE, na.rm = TRUE),
  has_sp     = sum(!is.na(signal_peptide_start)),
  has_tm     = sum(n_tm_helices > 0),
  avg_full   = round(mean(seq_length, na.rm = TRUE)),
  avg_current = round(mean(current_target_residues)),
  avg_after  = round(mean(ifelse(
    needs_correction == TRUE & !is.na(suggested_length),
    as.numeric(suggested_length), current_target_residues
  )))
), by = protein_type][order(-n)]

# ---------------------------------------------------------------------------
# Table 2: Correction Impact
# ---------------------------------------------------------------------------
correction_impact <- dt[needs_correction == TRUE, .(
  n          = .N,
  total_saved = sum(as.numeric(residues_saved), na.rm = TRUE),
  avg_saved  = round(mean(as.numeric(residues_saved), na.rm = TRUE)),
  med_saved  = round(median(as.numeric(residues_saved), na.rm = TRUE)),
  max_saved  = max(as.numeric(residues_saved), na.rm = TRUE),
  wall_drops = sum(walltime_changed, na.rm = TRUE)
), by = protein_type][order(-total_saved)]

# ---------------------------------------------------------------------------
# Table 3: Top 30 biggest corrections
# ---------------------------------------------------------------------------
biggest <- dt[needs_correction == TRUE][order(-as.numeric(residues_saved))][
  1:min(30, .N), .(
    target, uniprot, type = protein_type,
    full = seq_length, current = current_target_residues,
    corrected = suggested_length, saved = residues_saved,
    region_now = current_region, region_fix = suggested_region,
    wall_now = walltime_current, wall_fix = walltime_corrected
)]

# ---------------------------------------------------------------------------
# Table 4: Signal peptide audit (soluble with SP included)
# ---------------------------------------------------------------------------
sp_audit <- dt[protein_type == "soluble" & needs_correction == TRUE, .(
  target, uniprot, full = seq_length,
  sp = paste0(signal_peptide_start, "-", signal_peptide_end),
  sp_len = signal_peptide_end - signal_peptide_start + 1,
  current = current_target_residues,
  fix = suggested_length, fix_region = suggested_region
)][order(-sp_len)]

# ---------------------------------------------------------------------------
# Table 5: TM proteins at full-length (top 40)
# ---------------------------------------------------------------------------
tm_full <- dt[protein_type %in% c("type_i_tm", "multi_tm") &
              current_region == "full-length", .(
  target, uniprot, type = protein_type,
  full = seq_length, n_tm = n_tm_helices,
  tm1 = paste0(first_tm_start, "-", first_tm_end),
  ecto = paste0(ecto_start, "-", ecto_end), ecto_len = ecto_length,
  current = current_target_residues, saved = residues_saved,
  fix_region = suggested_region
)][order(-as.numeric(saved))][1:min(40, .N)]

# ---------------------------------------------------------------------------
# Table 6: GPI-anchored proteins
# ---------------------------------------------------------------------------
gpi <- dt[protein_type == "gpi_anchored", .(
  target, uniprot, full = seq_length,
  gpi = gpi_anchor_site,
  chain = paste0(chain_start, "-", chain_end),
  current = current_target_residues,
  needs_fix = needs_correction,
  fix_region = suggested_region, saved = residues_saved
)][order(-as.numeric(saved))]

# ---------------------------------------------------------------------------
# Table 7: Walltime changes
# ---------------------------------------------------------------------------
wall_changes <- dt[walltime_changed == TRUE, .(
  target, uniprot, type = protein_type,
  total_now = total_current, total_fix = total_corrected,
  wall_now = walltime_current, wall_fix = walltime_corrected,
  saved = residues_saved
)][order(-as.numeric(saved))][1:min(30, .N)]

# ---------------------------------------------------------------------------
# Table 8: Full decision table
# ---------------------------------------------------------------------------
full_table <- dt[, .(
  target, uniprot, type = protein_type, full = seq_length,
  sp = ifelse(!is.na(signal_peptide_start),
    paste0(signal_peptide_start, "-", signal_peptide_end), "-"),
  n_tm = n_tm_helices,
  region_now = current_region, res_now = current_target_residues,
  fix = needs_correction,
  region_fix = ifelse(!is.na(suggested_region), suggested_region, "-"),
  res_fix = ifelse(!is.na(suggested_length), suggested_length, current_target_residues),
  saved = ifelse(!is.na(residues_saved), residues_saved, 0L),
  total_fix = total_corrected,
  wall = walltime_corrected
)][order(type, -as.numeric(saved))]

# ---------------------------------------------------------------------------
# Write report
# ---------------------------------------------------------------------------
report_path <- file.path(output_dir, "topology_audit_report.txt")
con <- file(report_path, "w")

cat("TOPOLOGY AUDIT REPORT\n", file = con)
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n", file = con)
cat("Targets:", nrow(dt), "\n", file = con)
cat("Need correction:", sum(dt$needs_correction == TRUE, na.rm = TRUE), "/", nrow(dt), "\n", file = con)
cat("Total residues saved:", sum(as.numeric(dt$residues_saved), na.rm = TRUE), "\n\n", file = con)

for (tbl_info in list(
  list(type_summary,  "1. PROTEIN TYPE BREAKDOWN"),
  list(correction_impact, "2. CORRECTION IMPACT BY TYPE"),
  list(biggest,       "3. TOP 30 BIGGEST CORRECTIONS"),
  list(sp_audit,      "4. SIGNAL PEPTIDE AUDIT (soluble proteins)"),
  list(tm_full,       "5. TM PROTEINS AT FULL-LENGTH (top 40)"),
  list(gpi,           "6. GPI-ANCHORED PROTEINS"),
  list(wall_changes,  "7. WALLTIME CHANGES (top 30)")
)) {
  banner(tbl_info[[2]], con)
  lines <- capture.output(print(tbl_info[[1]], nrows = 500))
  cat(paste(lines, collapse = "\n"), "\n", file = con)
}

close(con)

# CSVs
fwrite(full_table, file.path(output_dir, "full_decision_table.csv"))
fwrite(dt[needs_correction == TRUE, .(
  target, uniprot, protein_type,
  current_region, current_target_residues,
  suggested_region, suggested_length, residues_saved
)], file.path(output_dir, "corrections_needed.csv"))

cat("Report:", report_path, "\n")
cat("CSV:", file.path(output_dir, "full_decision_table.csv"), "\n")
cat("CSV:", file.path(output_dir, "corrections_needed.csv"), "\n")

# ---------------------------------------------------------------------------
# Stdout summaries
# ---------------------------------------------------------------------------
banner("PROTEIN TYPE BREAKDOWN")
print(type_summary)

banner("CORRECTION IMPACT BY TYPE")
print(correction_impact)

banner("TOP 10 BIGGEST CORRECTIONS")
print(biggest[1:10])

banner("SIGNAL PEPTIDE AUDIT (soluble, top 10)")
print(sp_audit[1:min(10, nrow(sp_audit))])

cat("\n\n--- Summary ---\n")
cat("TM proteins at full-length:  ", nrow(dt[protein_type %in% c("type_i_tm","multi_tm") & current_region == "full-length"]), "\n")
cat("GPI proteins needing fix:    ", sum(gpi$needs_fix == TRUE, na.rm = TRUE), "\n")
cat("Soluble with SP included:    ", nrow(sp_audit), "\n")
cat("Total corrections needed:    ", sum(dt$needs_correction == TRUE, na.rm = TRUE), "\n")
cat("Total residues saved:        ", sum(as.numeric(dt$residues_saved), na.rm = TRUE), "\n")
cat("Targets with walltime change:", sum(dt$walltime_changed, na.rm = TRUE), "\n")
