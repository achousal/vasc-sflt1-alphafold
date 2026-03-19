#!/usr/bin/env python3
"""run_structural.py -- Step 3 entrypoint: Structural Prediction Pipeline.

Usage:
    python run_structural.py --step select    # Select candidates
    python run_structural.py --step fetch     # Fetch sequences from UniProt
    python run_structural.py --step generate  # Generate LSF job scripts
    python run_structural.py --step parse     # Parse AlphaFold results
    python run_structural.py --step plot      # Generate plots
    python run_structural.py --step all       # Run select + fetch + generate
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# Add script directory to path for imports
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import importlib

# Import numbered modules (filenames start with digits, need importlib)
_mod01 = importlib.import_module("01_select_candidates")
_mod02 = importlib.import_module("02_fetch_sequences")
_mod03 = importlib.import_module("03_generate_lsf_jobs")
_mod04 = importlib.import_module("04_parse_results")
_mod05 = importlib.import_module("05_plot_results")

select_candidates = _mod01.select_candidates
fetch_all_sequences = _mod02.fetch_all_sequences
generate_lsf_scripts = _mod03.generate_lsf_scripts
parse_all_results = _mod04.parse_all_results
generate_stats_report = _mod05.generate_stats_report
generate_summary_table = _mod05.generate_summary_table
plot_iptm_barplot = _mod05.plot_iptm_barplot
plot_pae_heatmap = _mod05.plot_pae_heatmap
load_pae_matrix = _mod04.load_pae_matrix


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Structural Prediction Pipeline"
    )
    parser.add_argument(
        "--step",
        choices=["select", "fetch", "generate", "parse", "plot", "all"],
        default="all",
        help="Which pipeline step to run (default: all = select+fetch+generate)",
    )
    parser.add_argument(
        "--consensus-dir",
        type=Path,
        default=Path("results/01_cross_cohort_overlap"),
        help="Step 1 consensus output directory",
    )
    parser.add_argument(
        "--enrichment-dir",
        type=Path,
        default=Path("results/02_pathway_enrichment"),
        help="Step 2 enrichment output directory",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/03_structural_prediction"),
        help="Output directory for Step 3",
    )
    parser.add_argument(
        "--af-output-dir",
        type=Path,
        default=None,
        help="AlphaFold output directory (for parse/plot steps)",
    )
    parser.add_argument(
        "--n-candidates",
        type=int,
        default=20,
        help="Max candidates to select (excluding controls)",
    )
    parser.add_argument(
        "--project-account",
        type=str,
        default="",
        help="LSF project account (e.g., acc_vascbrain). Required on Minerva.",
    )
    parser.add_argument(
        "--gpu-type",
        type=str,
        default="",
        help="GPU resource constraint (e.g., a100, v100, h100nvl).",
    )
    parser.add_argument(
        "--hpc-root",
        type=str,
        default="",
        help="Absolute path to project root on HPC. Sets working directory in job scripts.",
    )
    args = parser.parse_args()

    results_dir = args.results_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    fasta_dir = results_dir / "fasta"
    jobs_dir = results_dir / "jobs"
    candidates_path = results_dir / "step03_candidates.csv"
    scores_path = results_dir / "step03_interaction_scores.csv"

    # Resolve af_dir once so both parse and plot steps can use it
    af_dir = args.af_output_dir or results_dir / "af_output"

    steps = (
        ["select", "fetch", "generate"] if args.step == "all"
        else [args.step]
    )

    for step in steps:
        if step == "select":
            logger.info("=== Step 3a: Selecting candidates ===")
            consensus_path = args.consensus_dir / "step01_consensus_proteins_pos.csv"
            if not consensus_path.exists():
                logger.error("Consensus file not found: %s", consensus_path)
                sys.exit(1)

            candidates = select_candidates(
                consensus_path, args.enrichment_dir, args.n_candidates
            )
            candidates.to_csv(candidates_path, index=False)
            logger.info("Saved: %s (%d candidates)", candidates_path.name, len(candidates))

            # Print summary
            for _, row in candidates.iterrows():
                logger.info(
                    "  %2d. %-25s %s  %s",
                    row["rank"], row["target"], row["uniprot"], row["rationale"],
                )

        elif step == "fetch":
            logger.info("=== Step 3b: Fetching sequences ===")
            if not candidates_path.exists():
                logger.error("Run --step select first")
                sys.exit(1)

            results = fetch_all_sequences(candidates_path, fasta_dir)

            # Report
            n_ok = sum(1 for v in results.values() if v["status"] == "ok")
            n_fail = sum(1 for v in results.values() if v["status"] != "ok")
            logger.info("Fetched %d/%d sequences (%d failures)", n_ok, len(results), n_fail)

        elif step == "generate":
            logger.info("=== Step 3c: Generating LSF jobs ===")
            if not candidates_path.exists():
                logger.error("Run --step select first")
                sys.exit(1)

            scripts = generate_lsf_scripts(
                candidates_path, fasta_dir, jobs_dir,
                project_account=args.project_account,
                gpu_type=args.gpu_type,
                hpc_root=args.hpc_root,
            )
            logger.info("Generated %d LSF scripts in %s", len(scripts), jobs_dir)
            logger.info(
                "To submit (orchestrated parallel): bash %s/submit_orchestrator.sh",
                jobs_dir,
            )
            logger.info(
                "Dry run: bash %s/submit_orchestrator.sh --dry-run", jobs_dir
            )

        elif step == "parse":
            logger.info("=== Step 3d: Parsing AlphaFold results ===")
            if not candidates_path.exists():
                logger.error("Run --step select first")
                sys.exit(1)

            if not af_dir.exists():
                logger.error("AlphaFold output directory not found: %s", af_dir)
                logger.info("Run AlphaFold jobs on HPC first, then point --af-output-dir here")
                sys.exit(1)

            scores = parse_all_results(candidates_path, af_dir)
            scores.to_csv(scores_path, index=False)
            logger.info("Saved: %s", scores_path.name)

        elif step == "plot":
            logger.info("=== Step 3e: Generating plots ===")
            if not scores_path.exists():
                logger.error("Run --step parse first")
                sys.exit(1)

            import pandas as pd
            scores = pd.read_csv(scores_path)
            candidates = pd.read_csv(candidates_path)

            plot_iptm_barplot(scores, results_dir / "step03_iptm_barplot.pdf")
            generate_summary_table(scores, results_dir / "step03_summary_table.csv")
            generate_stats_report(scores, candidates, results_dir / "step03_stats_report.txt")

            import json
            parsed = scores[scores["n_models_parsed"] > 0]
            logger.info("Generating PAE heatmaps for %d parsed candidates", len(parsed))
            for _, row in parsed.iterrows():
                target = row["target"]
                safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")
                result_dir = af_dir / safe_name / f"sflt1_vs_{safe_name}"
                if not result_dir.exists():
                    result_dir = af_dir / safe_name  # fallback
                ranking_path = result_dir / "ranking_debug.json"
                if not ranking_path.exists():
                    logger.warning("No ranking_debug.json for %s, skipping PAE heatmap", target)
                    continue
                try:
                    with open(ranking_path) as _f:
                        ranking_data = json.load(_f)
                    order_map = ranking_data.get("order", [])
                    best_model = order_map[0] if order_map else None
                except Exception as _e:
                    logger.warning("Could not read ranking_debug.json for %s: %s", target, _e)
                    continue
                if best_model is None:
                    logger.warning("Empty model order for %s, skipping PAE heatmap", target)
                    continue
                pae_candidates = [
                    result_dir / f"pae_{best_model}.json",
                    result_dir / f"result_{best_model}.pkl",
                ]
                pae_matrix = None
                for pae_path in pae_candidates:
                    if pae_path.exists():
                        pae_matrix = load_pae_matrix(pae_path)
                        if pae_matrix is not None:
                            break
                if pae_matrix is None:
                    logger.warning("No PAE data found for %s, skipping heatmap", target)
                    continue
                out_path = results_dir / f"step03_pae_heatmap_{safe_name}.pdf"
                plot_pae_heatmap(pae_matrix, sflt1_length=338, target_name=target, output_path=out_path)
                logger.info("Saved PAE heatmap: %s", out_path.name)

    logger.info("=== Step 3 complete ===")


if __name__ == "__main__":
    main()
