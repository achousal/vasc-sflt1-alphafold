#!/bin/bash
# run_preflight.sh -- Run all preflight checks before AF2 submission
#
# Usage:
#   bash analysis/checks/run_preflight.sh              # local checks
#   bash analysis/checks/run_preflight.sh --check-hpc  # on Minerva, also verify HPC paths

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHECK_HPC=""
if [[ "${1:-}" == "--check-hpc" ]]; then
    CHECK_HPC="--check-hpc-paths"
fi

passed=0
failed=0
results=()

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

echo "============================================================"
echo "AF2 Preflight Checks"
echo "============================================================"

run_check "FASTA integrity" "$SCRIPT_DIR/check_fasta_integrity.py"
run_check "Candidate completeness" "$SCRIPT_DIR/check_candidates_complete.py"
run_check "AF2 inputs" "$SCRIPT_DIR/check_af2_inputs.py" $CHECK_HPC

echo ""
echo "============================================================"
echo "SUMMARY: $passed passed, $failed failed"
echo "============================================================"
for r in "${results[@]}"; do
    echo "  $r"
done

if [ "$failed" -gt 0 ]; then
    exit 1
fi
