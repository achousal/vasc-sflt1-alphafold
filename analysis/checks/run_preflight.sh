#!/bin/bash
# run_preflight.sh -- Run all preflight checks before AF2 submission
#
# Usage:
#   bash analysis/checks/run_preflight.sh              # local checks
#   bash analysis/checks/run_preflight.sh --check-hpc  # on Minerva, also verify HPC paths
#   bash analysis/checks/run_preflight.sh --skip-uniprot  # offline, no UniProt fetch
#
# Writes report to analysis/checks/preflight_report.txt

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPORT="$SCRIPT_DIR/preflight_report.txt"
CHECK_HPC=""
SKIP_UNIPROT=""

for arg in "$@"; do
    case "$arg" in
        --check-hpc) CHECK_HPC="--check-hpc-paths" ;;
        --skip-uniprot) SKIP_UNIPROT="--skip-uniprot" ;;
    esac
done

passed=0
failed=0
results=()

# Tee all output to both terminal and report file
exec > >(tee "$REPORT") 2>&1

echo "============================================================"
echo "AF2 Preflight Report"
echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Host: $(hostname)"
echo "============================================================"

run_check() {
    local name="$1"; shift
    echo ""
    echo ">>> $name"
    if python3 "$@"; then
        passed=$((passed + 1))
        results+=("PASS  $name")
    else
        failed=$((failed + 1))
        results+=("FAIL  $name")
    fi
}

run_check "FASTA integrity" "$SCRIPT_DIR/check_fasta_integrity.py" $SKIP_UNIPROT
run_check "Candidate completeness" "$SCRIPT_DIR/check_candidates_complete.py"
run_check "AF2 inputs" "$SCRIPT_DIR/check_af2_inputs.py" $CHECK_HPC

echo ""
echo "============================================================"
echo "SUMMARY: $passed passed, $failed failed"
echo "============================================================"
for r in "${results[@]}"; do
    echo "  $r"
done

echo ""
echo "Report written to: $REPORT"

if [ "$failed" -gt 0 ]; then
    exit 1
fi
