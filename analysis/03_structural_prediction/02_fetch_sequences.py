"""02_fetch_sequences.py -- Fetch FASTA sequences from UniProt REST API.

Produces two-chain FASTA files for AlphaFold Multimer:
  Chain A: sFLT1 construct (D1-D3, D1-D6, or D1-D7)
  Chain B: candidate protein (extracellular region only)

Construct definitions use structural boundaries from 5T89 (Markovic-Mueller
2017, VEGF + FLT1 D1-D6 crystal structure at 4.0A) rather than UniProt Ig-core
annotations. Signal peptide (residues 1-26) is always excluded -- all published
structural studies use the mature protein starting at Ser27.

Partner protein trimming (priority order):
  1. Annotated extracellular domain from UniProt "Topological domain" features
  2. Inferred ectodomain: after signal peptide, before first TM helix
  3. Chain annotation for GPI-anchored proteins
  4. Mature protein (signal peptide removed) for soluble/secreted

Numbering: all residue numbers refer to P17948 (FLT1_HUMAN) precursor.
"""

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# FLT1 / sFLT1 construct definitions (P17948 precursor numbering)
# ---------------------------------------------------------------------------
SFLT1_UNIPROT = "P17948"


@dataclass(frozen=True)
class FLT1Construct:
    """sFLT1/FLT1 ectodomain construct definition."""
    name: str
    start: int
    end: int
    description: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def fasta_header_tag(self) -> str:
        return f"sFLT1_{self.name}"


CONSTRUCT_D1D3 = FLT1Construct(
    name="D1-D3", start=27, end=330,
    description="Mature sFLT1 D1-D3, structural boundaries (5T89)",
)
CONSTRUCT_D1D6 = FLT1Construct(
    name="D1-D6", start=27, end=657,
    description="Mature sFLT1 D1-D6, structural boundaries (5T89)",
)
CONSTRUCT_D1D7 = FLT1Construct(
    name="D1-D7", start=27, end=747,
    description="FLT1 full ectodomain D1-D7, UniProt D7 end",
)

CONSTRUCTS = {"d1d3": CONSTRUCT_D1D3, "d1d6": CONSTRUCT_D1D6, "d1d7": CONSTRUCT_D1D7}

SFLT1_D1D3_START = CONSTRUCT_D1D3.start
SFLT1_D1D3_END = CONSTRUCT_D1D3.end

# UniProt REST API
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

REQUEST_DELAY = 0.5
MAX_RETRIES = 3
RETRY_DELAY = 5

MIN_ECTODOMAIN_LENGTH = 50

# Subcellular locations that can interact with extracellular sFLT1
ACCESSIBLE_KEYWORDS = {
    "Cell membrane", "Secreted", "Cell surface",
    "Extracellular space", "Extracellular matrix",
}

# Locations that are definitely intracellular
INTRACELLULAR_KEYWORDS = {
    "Cytoplasm", "Nucleus", "Mitochondrion", "Endoplasmic reticulum",
    "Golgi apparatus", "Lysosome", "Peroxisome",
}


# ---------------------------------------------------------------------------
# UniProt fetch helpers
# ---------------------------------------------------------------------------
def _uniprot_get(url: str, as_json: bool = False):
    """GET with retries and rate-limit handling."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                return resp.json() if as_json else resp.text
            if resp.status_code == 429:
                logger.warning("Rate limited, waiting %ds (attempt %d)", RETRY_DELAY, attempt)
                time.sleep(RETRY_DELAY)
                continue
            logger.warning("HTTP %d for %s (attempt %d)", resp.status_code, url, attempt)
        except requests.RequestException as e:
            logger.warning("Request failed: %s (attempt %d)", e, attempt)
        time.sleep(RETRY_DELAY)
    return None


def fetch_fasta_from_uniprot(accession: str) -> str | None:
    """Fetch FASTA text from UniProt."""
    return _uniprot_get(UNIPROT_API_URL.format(accession=accession))


def parse_fasta_sequence(fasta_text: str) -> str:
    """Extract amino acid sequence from FASTA text."""
    lines = fasta_text.strip().split("\n")
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def trim_sequence(sequence: str, start: int, end: int) -> str:
    """Trim sequence to residue range (1-indexed, inclusive)."""
    return sequence[start - 1 : end]


# ---------------------------------------------------------------------------
# Protein topology
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ProteinTopology:
    """Protein topology from UniProt annotations."""
    accession: str
    seq_length: int
    signal_peptide: tuple[int, int] | None
    transmembrane: tuple[tuple[int, int], ...] = field(default_factory=tuple)
    extracellular: tuple[tuple[int, int], ...] = field(default_factory=tuple)
    lumenal: tuple[tuple[int, int], ...] = field(default_factory=tuple)
    chain: tuple[int, int] | None = None
    gpi_anchor: int | None = None
    protein_type: str = "soluble"
    subcellular_keywords: tuple[str, ...] = field(default_factory=tuple)

    @property
    def is_transmembrane(self) -> bool:
        return len(self.transmembrane) > 0

    @property
    def is_accessible(self) -> bool:
        """Can this protein interact with extracellular sFLT1?"""
        if self.is_transmembrane or self.protein_type == "gpi_anchored":
            return True
        kw_set = set(self.subcellular_keywords)
        if kw_set & ACCESSIBLE_KEYWORDS:
            return True
        if kw_set and kw_set <= INTRACELLULAR_KEYWORDS:
            return False
        # No keywords -- assume accessible (conservative)
        return True

    @property
    def is_lumenal_only(self) -> bool:
        """TM protein with lumenal domains but no extracellular annotation.
        These face the Golgi/ER interior, not the extracellular space.
        See ADR-004."""
        return (self.is_transmembrane
                and len(self.lumenal) > 0
                and len(self.extracellular) == 0)

    @property
    def extracellular_range(self) -> tuple[int, int] | None:
        """Best extracellular region to model. Annotated ECD first, inferred second.

        Uses the largest single annotated extracellular segment, not the
        union span. For multi-pass TM, the span would include TM helices
        and cytoplasmic loops between extracellular segments. See ADR-004.

        Returns (start, end) 1-indexed, or None.
        """
        # Strategy A: largest annotated extracellular segment
        if self.extracellular:
            best = max(self.extracellular, key=lambda r: r[1] - r[0])
            start, end = best
            if (end - start + 1) >= MIN_ECTODOMAIN_LENGTH:
                return (start, end)

        # Strategy B: infer from signal peptide → first TM helix
        if self.is_transmembrane:
            ecto_start = (self.signal_peptide[1] + 1) if self.signal_peptide else 1
            ecto_end = self.transmembrane[0][0] - 1
            if ecto_end >= ecto_start and (ecto_end - ecto_start + 1) >= MIN_ECTODOMAIN_LENGTH:
                return (ecto_start, ecto_end)

        return None

    @property
    def mature_range(self) -> tuple[int, int] | None:
        """Mature protein range (signal peptide removed)."""
        if self.chain:
            return self.chain
        if self.signal_peptide:
            return (self.signal_peptide[1] + 1, self.seq_length)
        return None

    # Keep old property for backward compat with tests
    @property
    def ectodomain_range(self) -> tuple[int, int] | None:
        return self.extracellular_range

    @property
    def ectodomain_length(self) -> int | None:
        r = self.extracellular_range
        return (r[1] - r[0] + 1) if r else None


def fetch_topology_from_uniprot(accession: str) -> ProteinTopology | None:
    """Fetch protein topology from UniProt JSON API."""
    data = _uniprot_get(UNIPROT_JSON_URL.format(accession=accession), as_json=True)
    if data is None:
        return None
    return _parse_topology(accession, data)


def _parse_topology(accession: str, data: dict) -> ProteinTopology:
    """Parse UniProt JSON into ProteinTopology."""
    seq_length = data.get("sequence", {}).get("length", 0)
    features = data.get("features", [])

    signal_peptide = None
    transmembrane = []
    extracellular = []
    lumenal = []
    chain = None
    gpi_anchor = None

    # Subcellular location keywords
    keywords = {kw.get("name", "") for kw in data.get("keywords", [])}
    subcell_kw = tuple(sorted(keywords & (ACCESSIBLE_KEYWORDS | INTRACELLULAR_KEYWORDS)))
    is_gpi = "GPI-anchor" in keywords

    for feat in features:
        ftype = feat.get("type", "")
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        desc = feat.get("description", "")

        if start is None or end is None:
            continue

        if ftype == "Signal":
            signal_peptide = (start, end)
        elif ftype == "Transmembrane":
            transmembrane.append((start, end))
        elif ftype == "Topological domain" and "Extracellular" in desc:
            extracellular.append((start, end))
        elif ftype == "Topological domain" and "Lumenal" in desc:
            lumenal.append((start, end))
        elif ftype == "Chain" and chain is None:
            chain = (start, end)
        elif ftype == "Lipidation" and "GPI" in desc:
            gpi_anchor = start

    if len(transmembrane) > 1:
        protein_type = "multi_tm"
    elif len(transmembrane) == 1:
        protein_type = "type_i_tm"
    elif is_gpi:
        protein_type = "gpi_anchored"
    else:
        protein_type = "soluble"

    return ProteinTopology(
        accession=accession,
        seq_length=seq_length,
        signal_peptide=signal_peptide,
        transmembrane=tuple(transmembrane),
        extracellular=tuple(extracellular),
        lumenal=tuple(lumenal),
        chain=chain,
        gpi_anchor=gpi_anchor,
        protein_type=protein_type,
        subcellular_keywords=subcell_kw,
    )


# ---------------------------------------------------------------------------
# sFLT1 sequence
# ---------------------------------------------------------------------------
def get_sflt1_sequence(construct: FLT1Construct | None = None) -> str | None:
    """Fetch and trim the sFLT1/FLT1 ectodomain sequence."""
    if construct is None:
        construct = CONSTRUCT_D1D3
    fasta_text = fetch_fasta_from_uniprot(SFLT1_UNIPROT)
    if fasta_text is None:
        return None
    full_seq = parse_fasta_sequence(fasta_text)
    trimmed = trim_sequence(full_seq, construct.start, construct.end)
    logger.info("sFLT1 %s: %d aa (residues %d-%d)", construct.name, len(trimmed),
                construct.start, construct.end)
    return trimmed


# ---------------------------------------------------------------------------
# FASTA output
# ---------------------------------------------------------------------------
def write_two_chain_fasta(
    sflt1_seq: str, partner_seq: str, target_name: str,
    partner_uniprot: str, output_path: Path,
    construct: FLT1Construct | None = None, partner_region: str = "",
) -> None:
    """Write a two-chain FASTA for AlphaFold Multimer."""
    if construct is None:
        construct = CONSTRUCT_D1D3

    partner_header = f">{target_name}|{partner_uniprot}"
    if partner_region:
        partner_header += f"|{partner_region}"

    with open(output_path, "w") as f:
        f.write(f">{construct.fasta_header_tag}|{SFLT1_UNIPROT}"
                f"|residues_{construct.start}-{construct.end}\n")
        for i in range(0, len(sflt1_seq), 80):
            f.write(sflt1_seq[i : i + 80] + "\n")
        f.write(partner_header + "\n")
        for i in range(0, len(partner_seq), 80):
            f.write(partner_seq[i : i + 80] + "\n")


# ---------------------------------------------------------------------------
# Partner sequence selection
# ---------------------------------------------------------------------------
def select_partner_region(
    full_seq: str, topo: ProteinTopology, target: str,
) -> tuple[str, str, str]:
    """Pick the right region of a partner protein to model.

    Returns (sequence, region_tag, status).
    Region tag goes in the FASTA header for provenance.
    """
    # Intracellular-only → skip
    if not topo.is_accessible:
        logger.warning("  %s: intracellular-only (%s), skipping",
                       target, ", ".join(topo.subcellular_keywords))
        return ("", "", "intracellular")

    # Lumenal-only TM → skip (Golgi/ER interior, not extracellular; ADR-004)
    if topo.is_lumenal_only:
        logger.warning("  %s: lumenal-only TM (%d aa lumenal, no extracellular), skipping",
                       target, sum(r[1] - r[0] + 1 for r in topo.lumenal))
        return ("", "", "lumenal_only")

    # TM or multi-TM → extracellular region
    if topo.is_transmembrane:
        ecd = topo.extracellular_range
        if ecd is not None:
            seq = trim_sequence(full_seq, ecd[0], ecd[1])
            src = "annotated" if topo.extracellular else "inferred"
            logger.info("  %s: %s, %s ECD %d-%d (%d aa)",
                        target, topo.protein_type, src, ecd[0], ecd[1], len(seq))
            return (seq, f"ecd_{ecd[0]}-{ecd[1]}", "ok")
        else:
            logger.warning("  %s: TM but ECD < %d aa, skipping",
                           target, MIN_ECTODOMAIN_LENGTH)
            return ("", "", "ecd_too_short")

    # GPI-anchored → Chain annotation (excludes SP + GPI tail)
    if topo.protein_type == "gpi_anchored":
        if topo.chain is not None:
            seq = trim_sequence(full_seq, topo.chain[0], topo.chain[1])
            logger.info("  %s: GPI, chain %d-%d (%d aa)",
                        target, topo.chain[0], topo.chain[1], len(seq))
            return (seq, f"mature_{topo.chain[0]}-{topo.chain[1]}", "ok")
        elif topo.signal_peptide:
            start = topo.signal_peptide[1] + 1
            seq = trim_sequence(full_seq, start, topo.seq_length)
            logger.info("  %s: GPI, SP removed, %d-%d (%d aa)",
                        target, start, topo.seq_length, len(seq))
            return (seq, f"mature_{start}-{topo.seq_length}", "ok")

    # Soluble → remove signal peptide
    if topo.signal_peptide is not None:
        start = topo.signal_peptide[1] + 1
        seq = trim_sequence(full_seq, start, topo.seq_length)
        logger.info("  %s: soluble, SP 1-%d removed (%d aa)",
                     target, topo.signal_peptide[1], len(seq))
        return (seq, f"mature_{start}-{topo.seq_length}", "ok")

    # No annotations → full-length
    logger.info("  %s: soluble, no SP, full-length (%d aa)",
                target, len(full_seq))
    return (full_seq, "full_length", "ok")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def fetch_all_sequences(
    candidates_path: Path,
    fasta_dir: Path,
    construct: FLT1Construct | None = None,
    apply_ectodomain_filter: bool = True,
) -> dict[str, dict]:
    """Fetch sequences for all candidates and write two-chain FASTAs.

    Parameters
    ----------
    candidates_path : Path
        CSV with 'target' and 'uniprot' columns.
    fasta_dir : Path
        Output directory for FASTA files.
    construct : FLT1Construct, optional
        sFLT1 construct for Chain A. Defaults to D1-D3.
    apply_ectodomain_filter : bool
        Truncate TM proteins to ECD, remove signal peptides. Default True.

    Returns
    -------
    dict[str, dict]
        Per-target results with status, lengths, and topology info.
    """
    import pandas as pd

    if construct is None:
        construct = CONSTRUCT_D1D3

    candidates = pd.read_csv(candidates_path)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching sFLT1 %s sequence...", construct.name)
    sflt1_seq = get_sflt1_sequence(construct)
    if sflt1_seq is None:
        logger.error("Failed to fetch sFLT1 sequence. Aborting.")
        return {}

    # Save construct reference
    tag = construct.name.lower().replace("-", "")
    ref_path = fasta_dir / f"sflt1_{tag}_reference.fasta"
    with open(ref_path, "w") as f:
        f.write(f">{construct.fasta_header_tag}|{SFLT1_UNIPROT}"
                f"|residues_{construct.start}-{construct.end}\n")
        for i in range(0, len(sflt1_seq), 80):
            f.write(sflt1_seq[i : i + 80] + "\n")

    results = {}
    for _, row in candidates.iterrows():
        target = row["target"]
        uniprot = row["uniprot"]

        if not isinstance(uniprot, str) or len(uniprot) < 3:
            logger.warning("Skipping %s: invalid accession '%s'", target, uniprot)
            results[target] = {
                "uniprot": uniprot, "fasta_path": None,
                "partner_length": 0, "full_length": 0, "total_residues": 0,
                "status": "no_uniprot", "topology": "unknown", "region": "",
            }
            continue

        logger.info("Fetching %s (%s)...", target, uniprot)
        fasta_text = fetch_fasta_from_uniprot(uniprot)
        time.sleep(REQUEST_DELAY)

        if fasta_text is None:
            results[target] = {
                "uniprot": uniprot, "fasta_path": None,
                "partner_length": 0, "full_length": 0, "total_residues": 0,
                "status": "fetch_failed", "topology": "unknown", "region": "",
            }
            continue

        full_seq = parse_fasta_sequence(fasta_text)

        if apply_ectodomain_filter:
            topo = fetch_topology_from_uniprot(uniprot)
            time.sleep(REQUEST_DELAY)

            if topo is not None:
                partner_seq, partner_region, status = select_partner_region(
                    full_seq, topo, target
                )
                topology_type = topo.protein_type

                if status in ("intracellular", "ecd_too_short", "lumenal_only"):
                    results[target] = {
                        "uniprot": uniprot, "fasta_path": None,
                        "partner_length": 0, "full_length": len(full_seq),
                        "total_residues": 0, "status": status,
                        "topology": topology_type, "region": "",
                    }
                    continue
            else:
                partner_seq = full_seq
                partner_region = "full_length"
                topology_type = "unknown"
        else:
            partner_seq = full_seq
            partner_region = "full_length"
            topology_type = "unfiltered"

        safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")
        fasta_path = fasta_dir / f"sflt1_vs_{safe_name}.fasta"

        write_two_chain_fasta(
            sflt1_seq, partner_seq, target, uniprot, fasta_path,
            construct=construct, partner_region=partner_region,
        )

        total = construct.length + len(partner_seq)
        results[target] = {
            "uniprot": uniprot,
            "fasta_path": str(fasta_path),
            "partner_length": len(partner_seq),
            "full_length": len(full_seq),
            "total_residues": total,
            "status": "ok",
            "topology": topology_type,
            "region": partner_region,
        }

    # Summary
    ok = sum(1 for r in results.values() if r["status"] == "ok")
    skip_ic = sum(1 for r in results.values() if r["status"] == "intracellular")
    skip_short = sum(1 for r in results.values() if r["status"] == "ecd_too_short")
    skip_lum = sum(1 for r in results.values() if r["status"] == "lumenal_only")
    logger.info("Done: %d ok, %d lumenal_only, %d ecd_too_short, %d other",
                ok, skip_lum, skip_short,
                len(results) - ok - skip_lum - skip_short)

    return results
