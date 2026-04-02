#!/bin/bash
set -euo pipefail
#
# 11_batch_extract_and_merge.sh
#
# Run on Minerva (login node, no GPU needed).
# Walks all completed AF2 targets, extracts scores.json from pkl files
# via 09_extract_and_cleanup.py, then merges into a unified CSV via
# 10_merge_scores.py.
#
# NOTE: 09_extract_and_cleanup.py deletes result_model_*.pkl files after
# extraction to reclaim disk. features.pkl is retained. PDB structures
# and scores.json are retained. If you need result pkl files for
# re-analysis, back them up first.
#
# Usage:
#   bash analysis/03_structural_prediction/11_batch_extract_and_merge.sh [--force]
#
#   --force   Re-extract even if scores.json already exists (e.g., after
#             upgrading 09_extract_and_cleanup.py to add new metrics).
#
# Run from the project root:
#   /sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold

FORCE=0
if [[ "${1:-}" == "--force" ]]; then
    FORCE=1
    echo "[$(date '+%F %T')] FORCE mode: re-extracting all targets"
fi

HPC_ROOT="/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold"
SCRIPT_DIR="${HPC_ROOT}/analysis/03_structural_prediction"
RESULTS_DIR="${HPC_ROOT}/results/03_structural_prediction/d1d3"
CANDIDATES="${RESULTS_DIR}/candidates_unified.csv"
OUTPUT_CSV="${RESULTS_DIR}/step03_interaction_scores_all.csv"

SFLT1_LENGTH=304  # D1-D3 construct

# Ensure numpy is available (needed by 09_extract_and_cleanup.py)
module purge 2>/dev/null || true
module load python/3.10.2 2>/dev/null || true

cd "${HPC_ROOT}"

echo "[$(date '+%F %T')] Batch extract + merge"
echo "Results dir: ${RESULTS_DIR}"
echo "Candidates:  ${CANDIDATES}"
echo ""

# ---------------------------------------------------------------------------
# Step 1: Extract scores.json from completed targets
# ---------------------------------------------------------------------------
n_extracted=0
n_skipped=0
n_no_ranking=0

for target_dir in "${RESULTS_DIR}"/results/*/; do
    [ -d "$target_dir" ] || continue
    target_name=$(basename "$target_dir")

    # AF2 puts output in a subdirectory named after the FASTA stem
    af_subdir=$(find "$target_dir" -name "ranking_debug.json" -printf "%h" -quit 2>/dev/null)

    if [ -z "$af_subdir" ]; then
        ((n_no_ranking++)) || true
        continue
    fi

    if [ "$FORCE" -eq 0 ] && [ -f "${af_subdir}/scores.json" ]; then
        ((n_skipped++)) || true
        continue
    fi

    echo "  Extracting: ${target_name}"
    python3 "${SCRIPT_DIR}/09_extract_and_cleanup.py" "$af_subdir" "$SFLT1_LENGTH"
    ((n_extracted++)) || true
done

echo ""
echo "[$(date '+%F %T')] Extraction complete"
echo "  Extracted: ${n_extracted}"
echo "  Already had scores.json: ${n_skipped}"
echo "  No ranking_debug.json (incomplete/failed): ${n_no_ranking}"
echo ""

# ---------------------------------------------------------------------------
# Step 2: Merge all scores.json into unified CSV
# ---------------------------------------------------------------------------

# The merge script expects the AF output dir (parent of per-target subdirs)
# and the candidates CSV
if [ ! -f "$CANDIDATES" ]; then
    echo "WARN: candidates_unified.csv not found at ${CANDIDATES}"
    echo "Falling back to fullscreen candidates"
    CANDIDATES="${HPC_ROOT}/results/03_structural_prediction/d1d3_fullscreen/candidates_unified.csv"
fi

echo "[$(date '+%F %T')] Merging scores -> ${OUTPUT_CSV}"
python3 "${SCRIPT_DIR}/10_merge_scores.py" \
    "${RESULTS_DIR}/results" \
    "${CANDIDATES}" \
    "${OUTPUT_CSV}"

echo ""
echo "[$(date '+%F %T')] Done. Output: ${OUTPUT_CSV}"
echo ""

# Quick summary (header-based column lookup)
if [ -f "$OUTPUT_CSV" ]; then
    n_total=$(tail -n +2 "$OUTPUT_CSV" | wc -l)

    # Find iptm_best column index from header
    iptm_col=$(head -1 "$OUTPUT_CSV" | tr ',' '\n' | grep -n "^iptm_best$" | cut -d: -f1)
    call_col=$(head -1 "$OUTPUT_CSV" | tr ',' '\n' | grep -n "^interaction_call$" | cut -d: -f1)
    bias_col=$(head -1 "$OUTPUT_CSV" | tr ',' '\n' | grep -n "^template_bias_flag$" | cut -d: -f1 || echo "")

    n_scored=0
    n_above_06=0
    n_biased=0

    if [ -n "$call_col" ]; then
        n_scored=$(awk -F, -v c="$call_col" 'NR>1 && $c != "not_run" && $c != "" {n++} END{print n+0}' "$OUTPUT_CSV")
    fi
    if [ -n "$iptm_col" ]; then
        n_above_06=$(awk -F, -v c="$iptm_col" 'NR>1 && $c+0 > 0.6 {n++} END{print n+0}' "$OUTPUT_CSV")
    fi
    if [ -n "$bias_col" ]; then
        n_biased=$(awk -F, -v c="$bias_col" 'NR>1 && ($c == "True" || $c == "true" || $c == "1") {n++} END{print n+0}' "$OUTPUT_CSV")
    fi

    echo "Summary:"
    echo "  Total targets:         ${n_total}"
    echo "  With scores:           ${n_scored}"
    echo "  ipTM > 0.6:            ${n_above_06}"
    [ -n "$bias_col" ] && echo "  Template bias flagged:  ${n_biased}"
fi
