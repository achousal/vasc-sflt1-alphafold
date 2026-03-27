"""08_cleanup_pkl.py -- Delete AlphaFold pkl files after parse + render are complete.

PKL files are the largest AF2 outputs (~100-500MB each, 5 per target).
After step 04 (parse) and step 06 (render), all useful data has been
extracted to CSV and JPEG. This step reclaims disk space on HPC.

Safety: only deletes pkl files for targets that have both:
  1. A row in interaction_scores.csv with n_models_parsed > 0
  2. A render entry (structure JPEG exists) OR --skip-render-check flag
"""

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def find_pkl_files(af_output_dir: Path, target: str) -> list[Path]:
    """Find all result pkl files for a given target.

    Parameters
    ----------
    af_output_dir : Path
        Root AF2 output directory.
    target : str
        Target name (will be sanitized for filesystem).

    Returns
    -------
    list[Path]
        Paths to result_model_*.pkl files.
    """
    safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")
    target_dir = af_output_dir / safe_name

    # AF2 nests in sflt1_vs_{name} subdir
    nested_dir = target_dir / f"sflt1_vs_{safe_name}"
    actual_dir = nested_dir if nested_dir.exists() else target_dir

    if not actual_dir.exists():
        return []

    return sorted(actual_dir.glob("result_model_*.pkl"))


def cleanup_pkl_files(
    scores_path: Path,
    af_output_dir: Path,
    render_dir: Path | None = None,
    dry_run: bool = False,
    skip_render_check: bool = False,
) -> dict[str, list[Path]]:
    """Delete pkl files for targets that have been fully parsed and rendered.

    Parameters
    ----------
    scores_path : Path
        Path to step03_interaction_scores.csv (from step 04).
    af_output_dir : Path
        Root AF2 output directory containing per-target subdirs.
    render_dir : Path or None
        Directory containing rendered JPEGs (from step 06).
        Required unless skip_render_check is True.
    dry_run : bool
        If True, log what would be deleted without actually deleting.
    skip_render_check : bool
        If True, skip checking for rendered JPEGs. Use when only parse
        results are needed and renders will not be re-generated.

    Returns
    -------
    dict[str, list[Path]]
        Mapping of outcome -> list of pkl paths.
        Keys: "deleted", "skipped_no_parse", "skipped_no_render", "skipped_missing".
    """
    scores = pd.read_csv(scores_path)

    result: dict[str, list[Path]] = {
        "deleted": [],
        "skipped_no_parse": [],
        "skipped_no_render": [],
        "skipped_missing": [],
    }

    for _, row in scores.iterrows():
        target = row["target"]
        safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")
        n_parsed = row.get("n_models_parsed", 0)

        # Gate 1: parse must have succeeded
        if pd.isna(n_parsed) or int(n_parsed) == 0:
            pkl_files = find_pkl_files(af_output_dir, target)
            result["skipped_no_parse"].extend(pkl_files)
            if pkl_files:
                logger.info("SKIP (no parse): %s (%d pkl files)", target, len(pkl_files))
            continue

        # Gate 2: render check (unless skipped)
        if not skip_render_check:
            if render_dir is None:
                logger.error("render_dir required when skip_render_check is False")
                return result

            struct_jpg = render_dir / f"{safe_name}_structure.jpg"
            pae_jpg = render_dir / f"{safe_name}_pae.jpg"

            if not struct_jpg.exists() and not pae_jpg.exists():
                pkl_files = find_pkl_files(af_output_dir, target)
                result["skipped_no_render"].extend(pkl_files)
                if pkl_files:
                    logger.info("SKIP (no render): %s (%d pkl files)", target, len(pkl_files))
                continue

        # Both gates passed -- delete pkl files
        pkl_files = find_pkl_files(af_output_dir, target)
        if not pkl_files:
            result["skipped_missing"].append(target)
            continue

        total_bytes = sum(f.stat().st_size for f in pkl_files)
        total_mb = total_bytes / (1024 * 1024)

        if dry_run:
            logger.info(
                "DRY RUN: would delete %d pkl files for %s (%.1f MB)",
                len(pkl_files), target, total_mb,
            )
        else:
            for pkl in pkl_files:
                pkl.unlink()
                logger.debug("Deleted: %s", pkl)
            logger.info(
                "Deleted %d pkl files for %s (%.1f MB freed)",
                len(pkl_files), target, total_mb,
            )

        result["deleted"].extend(pkl_files)

    # Summary
    n_deleted = len(result["deleted"])
    n_skip_parse = len(result["skipped_no_parse"])
    n_skip_render = len(result["skipped_no_render"])
    action = "Would delete" if dry_run else "Deleted"
    logger.info(
        "Cleanup summary: %s %d pkl files, skipped %d (no parse), %d (no render)",
        action, n_deleted, n_skip_parse, n_skip_render,
    )

    return result
