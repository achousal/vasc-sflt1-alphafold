#!/usr/bin/env python3
"""07_generate_fullscreen_batch.py -- Generate Tier 1 full-screen AF2 batch.

Scales the AF2 Multimer pipeline from the curated protein list to all Tier 1
consensus targets. Supports all sFLT1 constructs (D1-D3, D1-D6, D1-D7).

Output:
  results/03_structural_prediction/{construct}_fullscreen/
    candidates.csv
    batch_summary.json
    fasta/          ~340 two-chain FASTA files
    jobs/           ~340 .lsf scripts + submit_all.sh + manifest.json + wrapper.sh

Usage:
    python 07_generate_fullscreen_batch.py --construct d1d3 --include-existing
    python 07_generate_fullscreen_batch.py --construct d1d6 --include-existing
    python 07_generate_fullscreen_batch.py --construct d1d7 --include-existing
    python 07_generate_fullscreen_batch.py --construct d1d3 --dry-run
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import importlib

_mod02 = importlib.import_module("02_fetch_sequences")
_mod03 = importlib.import_module("03_generate_lsf_jobs")

CONSTRUCTS = _mod02.CONSTRUCTS
fetch_all_sequences = _mod02.fetch_all_sequences


# ---------------------------------------------------------------------------
# Walltime heuristics (A100, AF2 2.3.2, 25 models)
# ---------------------------------------------------------------------------

def _walltime_for_residues(total_residues: int) -> str:
    """Assign walltime based on total complex residues."""
    if total_residues < 800:
        return "48:00"
    if total_residues < 1200:
        return "72:00"
    if total_residues < 1800:
        return "96:00"
    return "144:00"


# ---------------------------------------------------------------------------
# Target selection
# ---------------------------------------------------------------------------

def select_tier1_targets(
    consensus_path: Path,
    existing_candidates_path: Path | None = None,
    include_existing: bool = False,
) -> pd.DataFrame:
    """Select all Tier 1 targets, optionally excluding already-run targets.

    Parameters
    ----------
    consensus_path : Path
        Path to step01_consensus_proteins_pos.csv.
    existing_candidates_path : Path or None
        Path to d1d3/candidates.csv from the 24-protein run. If provided and
        include_existing is False, these targets are excluded.
    include_existing : bool
        If True, include targets that overlap with the existing batch.

    Returns
    -------
    pd.DataFrame
        Deduplicated Tier 1 targets with columns from consensus CSV.
    """
    raw = pd.read_csv(consensus_path)

    # Deduplicate: keep best tier per target
    raw = raw.sort_values(["tier", "dual_somamer"], ascending=[True, False])
    raw = raw.drop_duplicates(subset="target", keep="first").reset_index(drop=True)

    # Filter to Tier 1
    tier1 = raw[raw["tier"] == 1].copy()
    logger.info("Tier 1 targets after dedup: %d", len(tier1))

    # Exclude existing batch targets
    if existing_candidates_path and existing_candidates_path.exists() and not include_existing:
        existing = pd.read_csv(existing_candidates_path)
        existing_targets = set(existing["target"])
        before = len(tier1)
        tier1 = tier1[~tier1["target"].isin(existing_targets)].reset_index(drop=True)
        excluded = before - len(tier1)
        logger.info("Excluded %d targets already in d1d3 batch -> %d new targets",
                     excluded, len(tier1))

    # Validate UniProt coverage
    missing_up = tier1[tier1["uniprot"].isna() | (tier1["uniprot"].str.len() < 3)]
    if len(missing_up) > 0:
        logger.error("Targets with missing UniProt: %s", missing_up["target"].tolist())
        logger.error("Run patch_missing_uniprot.py first")
        sys.exit(1)

    return tier1


# ---------------------------------------------------------------------------
# Batch generation (follows 06_generate_ec_batches.py pattern)
# ---------------------------------------------------------------------------

def generate_fullscreen_batch(
    targets_df: pd.DataFrame,
    base_dir: Path,
    construct_key: str = "d1d3",
    hpc_root: str = "",
    project_account: str = "acc_vascbrain",
    gpu_type: str = "",
    dry_run: bool = False,
) -> dict:
    """Generate FASTAs and LSF jobs for a fullscreen Tier 1 batch.

    Parameters
    ----------
    construct_key : str
        sFLT1 construct variant: d1d3, d1d6, or d1d7.

    Returns summary dict with counts and walltime estimates.
    """
    construct = CONSTRUCTS[construct_key]
    batch_name = f"{construct_key}_fullscreen"

    batch_dir = base_dir / batch_name
    fasta_dir = batch_dir / "fasta"
    jobs_dir = batch_dir / "jobs"
    candidates_csv = batch_dir / "candidates.csv"

    batch_dir.mkdir(parents=True, exist_ok=True)

    n_targets = len(targets_df)

    logger.info("=" * 60)
    logger.info("Batch: %s (%s, %d targets)", batch_name, construct.name, n_targets)
    logger.info("  Construct: %s (%d aa, residues %d-%d)",
                construct.name, construct.length, construct.start, construct.end)
    logger.info("=" * 60)

    # 1. Write candidates CSV
    targets_df.to_csv(candidates_csv, index=False)
    logger.info("Wrote %d candidates to %s", n_targets, candidates_csv.name)

    # 2. Fetch sequences + write FASTAs
    results = fetch_all_sequences(candidates_csv, fasta_dir, construct=construct)

    ok_targets = {k: v for k, v in results.items() if v["status"] == "ok"}
    failed = {k: v for k, v in results.items() if v["status"] != "ok"}
    if failed:
        logger.warning("Failed fetches (%d): %s", len(failed), list(failed.keys()))

    logger.info("Fetched %d / %d sequences successfully", len(ok_targets), n_targets)

    # 3. Compute per-target walltimes
    walltime_map = {}
    for target, info in ok_targets.items():
        total = info["total_residues"]
        wt = _walltime_for_residues(total)
        walltime_map[target] = {"total_residues": total, "walltime": wt}

    # Walltime distribution summary
    wt_counts = {}
    total_gpu_h = 0
    for info in walltime_map.values():
        wt = info["walltime"]
        wt_counts[wt] = wt_counts.get(wt, 0) + 1
        total_gpu_h += int(wt.replace(":00", ""))

    logger.info("Walltime distribution:")
    for wt in sorted(wt_counts.keys()):
        h = int(wt.replace(":00", ""))
        n = wt_counts[wt]
        logger.info("  %s: %d jobs = %d GPU-h", wt, n, n * h)
    logger.info("  Total: %d jobs = %d GPU-h", len(ok_targets), total_gpu_h)

    if dry_run:
        logger.info("[DRY RUN] Would generate %d LSF jobs", len(ok_targets))
        return {
            "batch": batch_name, "construct": construct.name,
            "n_targets": len(ok_targets), "walltime_map": walltime_map,
            "total_gpu_h": total_gpu_h,
        }

    # 4. Generate LSF jobs with per-target walltimes
    hpc_batch_root = ""
    if hpc_root:
        hpc_batch_root = f"{hpc_root}/results/03_structural_prediction/{batch_name}"

    # generate_lsf_scripts uses compute_walltime() per target internally,
    # so each job gets the correct walltime without monkey-patching.
    scripts = _mod03.generate_lsf_scripts(
        candidates_csv, fasta_dir, jobs_dir,
        project_account=project_account,
        gpu_type=gpu_type,
        hpc_root=hpc_batch_root if hpc_batch_root else "",
        sflt1_length=construct.length,
    )

    # 6. Write batch summary
    summary = {
        "batch": batch_name,
        "construct": construct.name,
        "construct_residues": f"{construct.start}-{construct.end} ({construct.length} aa)",
        "n_targets": len(ok_targets),
        "n_scripts": len(scripts),
        "n_failed_fetches": len(failed),
        "failed_targets": list(failed.keys()),
        "total_gpu_h": total_gpu_h,
        "walltime_distribution": {wt: wt_counts[wt] for wt in sorted(wt_counts)},
        "walltime_map": walltime_map,
        "hpc_root": hpc_batch_root,
    }
    summary_path = batch_dir / "batch_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Batch summary: %s", summary_path)

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Generate Tier 1 full-screen AF2 batch (any construct)"
    )
    parser.add_argument(
        "--hpc-root",
        type=str,
        default="/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold",
        help="Absolute project root on HPC",
    )
    parser.add_argument("--project-account", type=str, default="acc_vascbrain")
    parser.add_argument("--gpu-type", type=str, default="",
                        help="GPU resource constraint. Defaults to module default (v100 || a100).")
    parser.add_argument("--construct", choices=["d1d3", "d1d6", "d1d7"], default="d1d3",
                        help="sFLT1 construct variant (default: d1d3)")
    parser.add_argument("--include-existing", action="store_true",
                        help="Include targets already in the original curated batch (re-run all)")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    project_root = SCRIPT_DIR.parent.parent
    results_dir = project_root / "results" / "03_structural_prediction"
    consensus_path = project_root / "results" / "01_cross_cohort_overlap" / "step01_consensus_proteins_pos.csv"
    existing_candidates = results_dir / args.construct / "candidates.csv"

    if not consensus_path.exists():
        logger.error("Consensus CSV not found: %s", consensus_path)
        sys.exit(1)

    # Select targets
    targets_df = select_tier1_targets(
        consensus_path,
        existing_candidates_path=existing_candidates,
        include_existing=args.include_existing,
    )

    if len(targets_df) == 0:
        logger.error("No targets selected. Check consensus CSV and existing batch.")
        sys.exit(1)

    # Generate batch
    summary = generate_fullscreen_batch(
        targets_df,
        base_dir=results_dir,
        construct_key=args.construct,
        hpc_root=args.hpc_root,
        project_account=args.project_account,
        gpu_type=args.gpu_type,
        dry_run=args.dry_run,
    )

    # Print next steps
    batch_name = summary["batch"]
    logger.info("")
    logger.info("=" * 60)
    logger.info("FULLSCREEN BATCH GENERATION COMPLETE (%s)", args.construct.upper())
    logger.info("=" * 60)
    logger.info("  Targets: %d", summary["n_targets"])
    logger.info("  GPU-hours: %d", summary["total_gpu_h"])
    logger.info("  At 8 A100s: ~%d days", summary["total_gpu_h"] // (8 * 24))
    logger.info("")
    logger.info("Next steps:")
    logger.info("  1. Validate: python analysis/checks/check_fullscreen_readiness.py")
    logger.info("  2. git add results/03_structural_prediction/%s/", batch_name)
    logger.info("  3. git commit -m 'feat(af2): fullscreen %s batch (%d Tier 1 targets)'",
                args.construct.upper(), summary["n_targets"])
    logger.info("  4. git push && ssh minerva 'cd <project> && git pull'")
    logger.info("  5. bash results/03_structural_prediction/%s/jobs/submit_all.sh --dry-run",
                batch_name)
    logger.info("  6. bash results/03_structural_prediction/%s/jobs/submit_all.sh",
                batch_name)


if __name__ == "__main__":
    main()
