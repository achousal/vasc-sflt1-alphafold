#!/usr/bin/env python3
"""check_fullscreen_readiness.py -- Preflight for scaling AF2 from 24 to all consensus targets.

Reads the full consensus CSV (step01_consensus_proteins_pos.csv), deduplicates,
and checks readiness for a full-screen AF2 run:

  1. UniProt accession coverage (flag missing/invalid)
  2. Duplicate target names after safe_name normalization
  3. Tier and dual-somamer composition
  4. Residue distribution estimate (fetch from UniProt topology API)
  5. Walltime budget projection per construct
  6. GPU-hour and storage estimates
  7. Walltime scaling bin assignment

Run:
  python analysis/checks/check_fullscreen_readiness.py                    # fast, no network
  python analysis/checks/check_fullscreen_readiness.py --fetch-topology   # fetches UniProt topology (~6 min)
  python analysis/checks/check_fullscreen_readiness.py --construct d1d3   # single construct estimate
"""

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CONSENSUS_PATH = (
    PROJECT_ROOT / "results" / "01_cross_cohort_overlap" / "step01_consensus_proteins_pos.csv"
)
RESULTS_DIR = PROJECT_ROOT / "results" / "03_structural_prediction"

# sFLT1 construct lengths (mature protein, signal peptide removed)
CONSTRUCTS = {
    "d1d3": {"name": "D1-D3", "start": 27, "end": 330, "length": 304},
    "d1d6": {"name": "D1-D6", "start": 27, "end": 657, "length": 631},
    "d1d7": {"name": "D1-D7", "start": 27, "end": 747, "length": 721},
}

# Walltime bins (total residues -> walltime hours)
WALLTIME_BINS = [
    (700, 48),
    (1200, 72),
    (1600, 96),
    (float("inf"), 144),
]

# Storage per target (empirical average from d1d3 24-protein run)
STORAGE_PER_TARGET_MB = 50

# UniProt REST API
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
REQUEST_DELAY = 0.3
MAX_RETRIES = 2
RETRY_DELAY = 3

# Minimum ectodomain length to model
MIN_ECTODOMAIN_LENGTH = 50


def safe_name(target: str) -> str:
    return target.replace("/", "-").replace(" ", "_").replace(":", "_")


def walltime_for_residues(total_residues: int) -> int:
    """Return walltime in hours based on total residue count."""
    for threshold, hours in WALLTIME_BINS:
        if total_residues <= threshold:
            return hours
    return WALLTIME_BINS[-1][1]


def deduplicate_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate consensus table: keep best tier per target."""
    df = df.sort_values(["tier", "dual_somamer"], ascending=[True, False])
    df = df.drop_duplicates(subset="target", keep="first")
    return df.reset_index(drop=True)


def fetch_topology(accession: str) -> dict | None:
    """Fetch protein topology from UniProt JSON API. Returns parsed dict or None."""
    import urllib.request
    import urllib.error

    url = UNIPROT_JSON_URL.format(accession=accession)
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=15) as resp:
                data = json.loads(resp.read().decode())
                return _parse_topology(data)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY)
    return None


def _parse_topology(data: dict) -> dict:
    """Parse UniProt JSON into a lightweight topology dict."""
    seq_length = data.get("sequence", {}).get("length", 0)
    features = data.get("features", [])
    keywords = {kw.get("name", "") for kw in data.get("keywords", [])}

    signal_end = None
    tm_regions = []
    chain_start = None
    chain_end = None

    for feat in features:
        ftype = feat.get("type", "")
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue

        if ftype == "Signal":
            signal_end = end
        elif ftype == "Transmembrane":
            tm_regions.append((start, end))
        elif ftype == "Chain" and chain_start is None:
            chain_start = start
            chain_end = end

    # Classify
    if len(tm_regions) > 1:
        ptype = "multi_tm"
    elif len(tm_regions) == 1:
        ptype = "type_i_tm"
    elif "GPI-anchor" in keywords:
        ptype = "gpi_anchored"
    else:
        ptype = "soluble"

    # Compute modeled length
    if ptype in ("type_i_tm", "multi_tm") and tm_regions:
        ecto_start = (signal_end + 1) if signal_end else 1
        ecto_end = tm_regions[0][0] - 1
        if ecto_end >= ecto_start and (ecto_end - ecto_start + 1) >= MIN_ECTODOMAIN_LENGTH:
            modeled_length = ecto_end - ecto_start + 1
        else:
            modeled_length = seq_length  # fallback to full-length
    elif ptype == "gpi_anchored" and chain_start and chain_end:
        modeled_length = chain_end - chain_start + 1
    elif signal_end:
        modeled_length = seq_length - signal_end
    else:
        modeled_length = seq_length

    return {
        "seq_length": seq_length,
        "modeled_length": modeled_length,
        "protein_type": ptype,
        "signal_end": signal_end,
        "n_tm": len(tm_regions),
    }


def check_uniprot_coverage(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Check UniProt accession coverage. Returns (failures, warnings)."""
    failures = []
    warnings = []

    missing = df[df["uniprot"].isna() | (df["uniprot"].str.len() < 3)]
    if len(missing) > 0:
        failures.append(
            f"{len(missing)} targets with missing/invalid UniProt accessions:"
        )
        for _, row in missing.iterrows():
            failures.append(f"    {row['target']} -> uniprot='{row.get('uniprot', 'NaN')}'")

    valid = df[df["uniprot"].notna() & (df["uniprot"].str.len() >= 3)]
    warnings.append(f"{len(valid)}/{len(df)} targets have valid UniProt accessions")

    return failures, warnings


def check_name_collisions(df: pd.DataFrame) -> list[str]:
    """Check for target name collisions after safe_name normalization."""
    failures = []
    df = df.copy()
    df["safe"] = df["target"].apply(safe_name)
    dupes = df[df["safe"].duplicated(keep=False)]
    if len(dupes) > 0:
        for sname, group in dupes.groupby("safe"):
            originals = group["target"].tolist()
            failures.append(
                f"Name collision after normalization: {originals} -> '{sname}'"
            )
    return failures


def check_composition(df: pd.DataFrame) -> list[str]:
    """Report tier and dual-somamer composition. Returns info lines."""
    lines = []
    cross = pd.crosstab(df["tier"], df["dual_somamer"])
    lines.append("  Tier x Dual-somamer composition:")
    lines.append(f"  {'':>8s} {'single':>8s} {'dual':>8s} {'total':>8s}")
    for tier in sorted(cross.index):
        single = cross.loc[tier, False] if False in cross.columns else 0
        dual = cross.loc[tier, True] if True in cross.columns else 0
        total = single + dual
        lines.append(f"  Tier {tier:>2d} {single:>8d} {dual:>8d} {total:>8d}")
    lines.append(f"  {'Total':>8s} {cross[False].sum() if False in cross.columns else 0:>8d} "
                 f"{cross[True].sum() if True in cross.columns else 0:>8d} {len(df):>8d}")
    return lines


def estimate_budget_from_topology(
    df: pd.DataFrame,
    construct_key: str = "d1d3",
) -> tuple[list[str], pd.DataFrame]:
    """Fetch UniProt topology and estimate GPU budget. Returns (info_lines, enriched_df)."""
    construct = CONSTRUCTS[construct_key]
    bait_len = construct["length"]
    lines = []
    results = []

    valid = df[df["uniprot"].notna() & (df["uniprot"].str.len() >= 3)].copy()
    n_total = len(valid)

    lines.append(f"  Fetching topology for {n_total} targets (est. {n_total * 0.6:.0f}s)...")

    for i, (_, row) in enumerate(valid.iterrows()):
        topo = fetch_topology(row["uniprot"])
        time.sleep(REQUEST_DELAY)

        if topo is None:
            results.append({
                "target": row["target"],
                "uniprot": row["uniprot"],
                "seq_length": 0,
                "modeled_length": 0,
                "total_residues": 0,
                "protein_type": "fetch_failed",
                "walltime_h": 0,
            })
            continue

        total = bait_len + topo["modeled_length"]
        wt = walltime_for_residues(total)
        results.append({
            "target": row["target"],
            "uniprot": row["uniprot"],
            "seq_length": topo["seq_length"],
            "modeled_length": topo["modeled_length"],
            "total_residues": total,
            "protein_type": topo["protein_type"],
            "walltime_h": wt,
        })

        if (i + 1) % 100 == 0:
            lines.append(f"    ... {i + 1}/{n_total} fetched")

    rdf = pd.DataFrame(results)

    # Filter out fetch failures
    ok = rdf[rdf["protein_type"] != "fetch_failed"]
    failed = rdf[rdf["protein_type"] == "fetch_failed"]
    if len(failed) > 0:
        lines.append(f"  WARNING: {len(failed)} targets failed UniProt fetch")

    # Protein type distribution
    lines.append(f"\n  Protein type distribution ({construct_key}):")
    for ptype, count in ok["protein_type"].value_counts().items():
        lines.append(f"    {ptype:<20s} {count:>5d}")

    # Residue distribution
    lines.append(f"\n  Modeled length distribution (partner only):")
    lines.append(f"    min:    {ok['modeled_length'].min():>6.0f} aa")
    lines.append(f"    median: {ok['modeled_length'].median():>6.0f} aa")
    lines.append(f"    mean:   {ok['modeled_length'].mean():>6.1f} aa")
    lines.append(f"    max:    {ok['modeled_length'].max():>6.0f} aa")

    lines.append(f"\n  Total residues (bait {bait_len} + partner):")
    lines.append(f"    min:    {ok['total_residues'].min():>6.0f}")
    lines.append(f"    median: {ok['total_residues'].median():>6.0f}")
    lines.append(f"    mean:   {ok['total_residues'].mean():>6.1f}")
    lines.append(f"    max:    {ok['total_residues'].max():>6.0f}")

    # Walltime bins
    lines.append(f"\n  Walltime bin distribution ({construct_key}):")
    for threshold, hours in WALLTIME_BINS:
        count = len(ok[ok["walltime_h"] == hours])
        gpu_h = count * hours
        label = f"<={threshold:.0f}" if threshold < float("inf") else f">{WALLTIME_BINS[-2][0]:.0f}"
        lines.append(f"    {hours:>4d}h ({label:>6s} aa): {count:>5d} jobs = {gpu_h:>8,d} GPU-h")

    total_gpu_h = ok["walltime_h"].sum()
    lines.append(f"    {'Total':>21s}: {len(ok):>5d} jobs = {total_gpu_h:>8,} GPU-h")

    # Storage estimate
    storage_gb = len(ok) * STORAGE_PER_TARGET_MB / 1024
    lines.append(f"\n  Storage estimate: {storage_gb:,.1f} GB ({len(ok)} targets x {STORAGE_PER_TARGET_MB} MB)")

    # Parallel execution timeline
    lines.append(f"\n  Parallel execution timeline ({construct_key} only):")
    for n_gpu in [8, 24, 64]:
        days = total_gpu_h / (n_gpu * 24)
        lines.append(f"    {n_gpu:>3d} A100s: ~{days:,.0f} days")

    return lines, rdf


def estimate_budget_offline(
    df: pd.DataFrame,
    construct_key: str = "d1d3",
) -> list[str]:
    """Estimate GPU budget using the empirical residue distribution from existing 24-protein runs."""
    construct = CONSTRUCTS[construct_key]
    bait_len = construct["length"]
    lines = []

    valid = df[df["uniprot"].notna() & (df["uniprot"].str.len() >= 3)]
    n = len(valid)

    # Use empirical average from d1d3 24-protein run
    # Mean partner length ~640 aa, mean total ~944 aa, mean walltime ~73h
    avg_partner = 640
    avg_total = bait_len + avg_partner
    avg_wt = walltime_for_residues(avg_total)

    lines.append(f"  Offline estimate (using empirical avg from 24-protein run):")
    lines.append(f"    Targets with valid UniProt: {n}")
    lines.append(f"    Avg partner length (empirical): ~{avg_partner} aa")
    lines.append(f"    Avg total residues: ~{avg_total}")
    lines.append(f"    Avg walltime bin: {avg_wt}h")

    # Conservative: assume same distribution as 24-protein run
    # 37.5% at 48h, 37.5% at 72h, 16.7% at 96h, 8.3% at 144h
    bins = [(0.375, 48), (0.375, 72), (0.167, 96), (0.083, 144)]
    total_gpu_h = sum(int(frac * n) * h for frac, h in bins)

    lines.append(f"\n  Projected walltime distribution ({construct_key}):")
    for frac, hours in bins:
        count = int(frac * n)
        gpu_h = count * hours
        lines.append(f"    {hours:>4d}h: ~{count:>5d} jobs = ~{gpu_h:>8,d} GPU-h")
    lines.append(f"    {'Total':>6s}: ~{n:>5d} jobs = ~{total_gpu_h:>8,} GPU-h")

    storage_gb = n * STORAGE_PER_TARGET_MB / 1024
    lines.append(f"\n  Storage estimate: ~{storage_gb:,.1f} GB")

    lines.append(f"\n  Parallel execution timeline ({construct_key} only):")
    for n_gpu in [8, 24, 64]:
        days = total_gpu_h / (n_gpu * 24)
        lines.append(f"    {n_gpu:>3d} A100s: ~{days:,.0f} days")

    # All 3 constructs
    total_3x = total_gpu_h * 3
    lines.append(f"\n  All 3 constructs (d1d3 + d1d6 + d1d7): ~{total_3x:,} GPU-h")
    for n_gpu in [8, 24, 64]:
        days = total_3x / (n_gpu * 24)
        lines.append(f"    {n_gpu:>3d} A100s: ~{days:,.0f} days")

    return lines


def check_existing_overlap(df: pd.DataFrame) -> list[str]:
    """Check how many full-screen targets overlap with the existing 24-protein run."""
    lines = []
    candidates_path = RESULTS_DIR / "step03_candidates.csv"
    if not candidates_path.exists():
        lines.append("  (step03_candidates.csv not found -- skipping overlap check)")
        return lines

    existing = pd.read_csv(candidates_path)
    existing_targets = set(existing["target"])
    all_targets = set(df["target"])

    overlap = existing_targets & all_targets
    new = all_targets - existing_targets

    lines.append(f"  Existing 24-protein run overlap:")
    lines.append(f"    Already run:  {len(overlap):>5d}")
    lines.append(f"    New targets:  {len(new):>5d}")
    lines.append(f"    Total screen: {len(all_targets):>5d}")

    return lines


def main():
    parser = argparse.ArgumentParser(
        description="Preflight check for full-consensus AF2 screen"
    )
    parser.add_argument(
        "--fetch-topology",
        action="store_true",
        help="Fetch UniProt topology for accurate residue/walltime estimates (~6 min)",
    )
    parser.add_argument(
        "--construct",
        choices=list(CONSTRUCTS.keys()),
        default="d1d3",
        help="Construct for budget estimate (default: d1d3)",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Write enriched target table with topology to CSV",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Full-Screen AF2 Preflight Check")
    print("=" * 60)

    # Load and deduplicate
    if not CONSENSUS_PATH.exists():
        print(f"\nFAIL: consensus file not found: {CONSENSUS_PATH}")
        sys.exit(1)

    raw = pd.read_csv(CONSENSUS_PATH)
    df = deduplicate_consensus(raw)
    print(f"\n  Consensus CSV: {len(raw)} rows -> {len(df)} unique targets")

    all_failures = []

    # 1. UniProt coverage
    print("\n--- UniProt Coverage ---")
    failures, warnings = check_uniprot_coverage(df)
    all_failures.extend(failures)
    for w in warnings:
        print(f"  {w}")
    for f in failures:
        print(f"  FAIL: {f}")

    # 2. Name collisions
    print("\n--- Name Collision Check ---")
    collisions = check_name_collisions(df)
    if collisions:
        all_failures.extend(collisions)
        for c in collisions:
            print(f"  FAIL: {c}")
    else:
        print("  PASS: no name collisions after safe_name normalization")

    # 3. Composition
    print("\n--- Composition ---")
    for line in check_composition(df):
        print(line)

    # 4. Overlap with existing run
    print("\n--- Existing Run Overlap ---")
    for line in check_existing_overlap(df):
        print(line)

    # 5. Budget estimate
    print(f"\n--- GPU Budget Estimate (construct: {args.construct}) ---")
    if args.fetch_topology:
        budget_lines, enriched_df = estimate_budget_from_topology(df, args.construct)
        for line in budget_lines:
            print(line)
        if args.output_csv:
            enriched_df.to_csv(args.output_csv, index=False)
            print(f"\n  Enriched target table written to: {args.output_csv}")
    else:
        for line in estimate_budget_offline(df, args.construct):
            print(line)
        print("\n  (run with --fetch-topology for accurate per-target estimates)")

    # Summary
    print("\n" + "=" * 60)
    if all_failures:
        print(f"FAIL: {len(all_failures)} issues found\n")
        for f in all_failures:
            print(f"  {f}")
        sys.exit(1)
    else:
        n_valid = len(df[df["uniprot"].notna() & (df["uniprot"].str.len() >= 3)])
        print(f"PASS: {n_valid} targets ready for full-screen AF2 ({len(df) - n_valid} need UniProt resolution)")
        sys.exit(0)


if __name__ == "__main__":
    main()
