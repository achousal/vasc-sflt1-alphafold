"""06_render_structures.py -- Render AF2 multimer predictions as JPEG images.

Uses PyMOL headless to produce per-target interface views:
  - Chain A (sFLT1) colored blue, chain B (partner) colored orange
  - Cartoon representation with transparent surface
  - Rotated to show binding interface
  - PAE-colored interface residues (low PAE = high confidence = green)

Outputs ~50-100KB JPEG per target vs ~5MB per PDB.
"""

import json
import logging
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# PyMOL deferred init -- only start when needed
_pymol_initialized = False


def _init_pymol():
    global _pymol_initialized
    if not _pymol_initialized:
        import pymol
        pymol.finish_launching(["pymol", "-cq"])  # headless, quiet
        _pymol_initialized = True


def _find_best_model(af_subdir: Path) -> tuple[Path | None, Path | None]:
    """Find the best-ranked unrelaxed PDB and its corresponding pkl.

    If ranking_debug.json exists, use it. Otherwise pick model_1_pred_0.
    """
    ranking_path = af_subdir / "ranking_debug.json"
    if ranking_path.exists():
        with open(ranking_path) as f:
            ranking = json.load(f)
        # ranking_debug.json has "order" key with model names sorted by iptm+ptm
        best_model = ranking["order"][0]
    else:
        # Fallback: pick first model found
        pkls = sorted(af_subdir.glob("result_model_*_pred_*.pkl"))
        if not pkls:
            return None, None
        best_model = pkls[0].stem.replace("result_", "")

    pdb_path = af_subdir / f"unrelaxed_{best_model}.pdb"
    pkl_path = af_subdir / f"result_{best_model}.pkl"

    if not pdb_path.exists():
        # Try relaxed if unrelaxed missing
        pdb_path = af_subdir / f"relaxed_{best_model}.pdb"

    if not pdb_path.exists():
        return None, pkl_path if pkl_path.exists() else None

    return pdb_path, pkl_path if pkl_path.exists() else None


def _load_pae_from_pkl(pkl_path: Path) -> np.ndarray | None:
    """Extract PAE matrix from AF2 result pkl."""
    try:
        with open(pkl_path, "rb") as f:
            result = pickle.load(f)
        if "predicted_aligned_error" in result:
            return np.array(result["predicted_aligned_error"])
    except Exception as e:
        logger.warning("Could not load PAE from %s: %s", pkl_path.name, e)
    return None


def _get_iptm_ptm(pkl_path: Path) -> tuple[float, float]:
    """Extract ipTM and pTM from AF2 result pkl."""
    try:
        with open(pkl_path, "rb") as f:
            result = pickle.load(f)
        iptm = float(result.get("iptm", 0.0))
        ptm = float(result.get("ptm", 0.0))
        return iptm, ptm
    except Exception:
        return 0.0, 0.0


def render_structure(
    pdb_path: Path,
    output_path: Path,
    target_name: str = "",
    pkl_path: Path | None = None,
    width: int = 1200,
    height: int = 900,
    dpi: int = 150,
) -> bool:
    """Render a single AF2 multimer PDB as JPEG.

    Parameters
    ----------
    pdb_path : Path
        Path to PDB file (unrelaxed or relaxed).
    output_path : Path
        JPEG output path.
    target_name : str
        Target name for title overlay.
    pkl_path : Path | None
        If provided, extract ipTM/pTM for annotation.
    width, height : int
        Image dimensions in pixels.
    dpi : int
        Output resolution.

    Returns
    -------
    bool
        True if rendering succeeded.
    """
    try:
        _init_pymol()
        from pymol import cmd

        cmd.reinitialize()
        cmd.load(str(pdb_path), "complex")

        # Color by chain
        cmd.color("marine", "chain A")  # sFLT1 = blue
        cmd.color("orange", "chain B")  # partner = orange

        # Cartoon with transparent surface
        cmd.show("cartoon", "complex")
        cmd.set("cartoon_transparency", 0.0)
        cmd.show("surface", "complex")
        cmd.set("transparency", 0.7)

        # Style
        cmd.set("ray_opaque_background", 1)
        cmd.bg_color("white")
        cmd.set("antialias", 2)
        cmd.set("ray_shadows", 0)

        # Orient to show interface
        cmd.orient("complex")
        # Rotate to show interface between chains
        cmd.rotate("y", 30)

        # Add title annotation
        iptm, ptm = (0.0, 0.0)
        if pkl_path and pkl_path.exists():
            iptm, ptm = _get_iptm_ptm(pkl_path)

        label = target_name or pdb_path.stem
        if iptm > 0:
            label = f"{label}  ipTM={iptm:.3f}  pTM={ptm:.3f}"

        # Render
        output_path.parent.mkdir(parents=True, exist_ok=True)
        cmd.set("ray_trace_mode", 1)
        cmd.ray(width, height)
        cmd.png(str(output_path.with_suffix(".png")), dpi=dpi)

        # Convert PNG to JPEG for smaller file size
        try:
            from PIL import Image

            img = Image.open(str(output_path.with_suffix(".png")))
            img = img.convert("RGB")
            img.save(str(output_path), "JPEG", quality=85, optimize=True)
            output_path.with_suffix(".png").unlink()  # remove PNG
        except ImportError:
            # No PIL -- keep PNG
            import shutil

            shutil.move(
                str(output_path.with_suffix(".png")), str(output_path)
            )

        logger.info("Rendered: %s (ipTM=%.3f)", output_path.name, iptm)
        return True

    except Exception as e:
        logger.error("Failed to render %s: %s", pdb_path.name, e)
        return False


def render_pae_heatmap(
    pkl_path: Path,
    output_path: Path,
    sflt1_length: int = 304,
    target_name: str = "",
) -> bool:
    """Render PAE matrix as a JPEG heatmap.

    Parameters
    ----------
    pkl_path : Path
        AF2 result pkl with predicted_aligned_error.
    output_path : Path
        JPEG output path.
    sflt1_length : int
        Chain A (sFLT1) length for domain boundary annotation.
    target_name : str
        Target name for title.

    Returns
    -------
    bool
        True if rendering succeeded.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pae = _load_pae_from_pkl(pkl_path)
    if pae is None:
        return False

    iptm, ptm = _get_iptm_ptm(pkl_path)

    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    im = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=30)
    plt.colorbar(im, ax=ax, label="PAE (Å)", shrink=0.8)

    # Chain boundary
    ax.axhline(y=sflt1_length - 0.5, color="red", linewidth=1, linestyle="--")
    ax.axvline(x=sflt1_length - 0.5, color="red", linewidth=1, linestyle="--")

    # Domain boundaries (D1-D3 construct, 0-indexed)
    for boundary in [103, 198]:
        ax.axhline(y=boundary - 0.5, color="gray", linewidth=0.5, alpha=0.5)
        ax.axvline(x=boundary - 0.5, color="gray", linewidth=0.5, alpha=0.5)

    title = f"sFLT1 vs {target_name}" if target_name else pkl_path.stem
    if iptm > 0:
        title += f"  (ipTM={iptm:.3f}, pTM={ptm:.3f})"
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Residue index")

    # Labels
    ax.text(
        sflt1_length / 2, -15, "sFLT1", ha="center", fontsize=9, color="navy"
    )
    partner_mid = sflt1_length + (pae.shape[0] - sflt1_length) / 2
    ax.text(
        partner_mid, -15, target_name, ha="center", fontsize=9, color="darkorange"
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), format="jpeg", dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Rendered PAE: %s", output_path.name)
    return True


def render_all(
    candidates_path: Path,
    af_output_dir: Path,
    render_dir: Path,
    sflt1_length: int = 304,
) -> pd.DataFrame:
    """Render all AF2 results for a batch.

    Produces two images per target:
      - {target}_structure.jpg  -- 3D cartoon+surface view
      - {target}_pae.jpg        -- PAE heatmap

    Parameters
    ----------
    candidates_path : Path
        CSV with target, uniprot columns.
    af_output_dir : Path
        AF2 output directory (contains per-target subdirs).
    render_dir : Path
        Output directory for JPEG images.
    sflt1_length : int
        Chain A length for PAE annotation.

    Returns
    -------
    pd.DataFrame
        Summary with target, iptm, ptm, structure_jpg, pae_jpg paths.
    """
    candidates = pd.read_csv(candidates_path)
    render_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for _, row in candidates.iterrows():
        target = row["target"]
        sname = target.replace("/", "-").replace(" ", "_").replace(":", "_")

        # Find AF2 output subdir
        target_dir = af_output_dir / sname
        if not target_dir.exists():
            logger.warning("No AF2 output for %s, skipping", target)
            rows.append({"target": target, "status": "missing"})
            continue

        # AF2 nests output in a subdir named after the FASTA stem
        af_subdirs = list(target_dir.glob("sflt1_vs_*"))
        if not af_subdirs:
            af_subdirs = [target_dir]
        af_subdir = af_subdirs[0]

        pdb_path, pkl_path = _find_best_model(af_subdir)
        if pdb_path is None:
            logger.warning("No model PDB for %s, skipping", target)
            rows.append({"target": target, "status": "no_model"})
            continue

        iptm, ptm = _get_iptm_ptm(pkl_path) if pkl_path else (0.0, 0.0)

        # Render structure
        struct_path = render_dir / f"{sname}_structure.jpg"
        struct_ok = render_structure(
            pdb_path, struct_path, target_name=target, pkl_path=pkl_path
        )

        # Render PAE
        pae_path = render_dir / f"{sname}_pae.jpg"
        pae_ok = False
        if pkl_path:
            pae_ok = render_pae_heatmap(
                pkl_path, pae_path, sflt1_length=sflt1_length, target_name=target
            )

        rows.append({
            "target": target,
            "iptm": iptm,
            "ptm": ptm,
            "structure_jpg": str(struct_path) if struct_ok else "",
            "pae_jpg": str(pae_path) if pae_ok else "",
            "status": "ok" if struct_ok else "render_failed",
        })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="Render AF2 multimer predictions as JPEG"
    )
    parser.add_argument(
        "--candidates", type=Path, required=True,
        help="Path to candidates CSV",
    )
    parser.add_argument(
        "--af-output-dir", type=Path, required=True,
        help="AF2 output directory (contains per-target subdirs)",
    )
    parser.add_argument(
        "--render-dir", type=Path, required=True,
        help="Output directory for JPEG images",
    )
    parser.add_argument(
        "--sflt1-length", type=int, default=304,
        help="sFLT1 chain A length (default: 304 for D1-D3)",
    )
    args = parser.parse_args()

    summary = render_all(
        args.candidates,
        args.af_output_dir,
        args.render_dir,
        sflt1_length=args.sflt1_length,
    )
    summary_path = args.render_dir / "render_summary.csv"
    summary.to_csv(summary_path, index=False)
    logger.info("Saved: %s (%d targets)", summary_path, len(summary))

    n_ok = (summary["status"] == "ok").sum()
    logger.info("Rendered %d/%d targets", n_ok, len(summary))
