"""02_fetch_sequences.py -- Fetch FASTA sequences from UniProt REST API.

Produces two-chain FASTA files for AlphaFold Multimer:
  Chain A: sFLT1 construct (D1-D3, D1-D6, or D1-D7)
  Chain B: candidate protein (ectodomain-only for transmembrane proteins)

Construct definitions use structural boundaries from 5T89 (Markovic-Mueller
2017, VEGF + FLT1 D1-D6 crystal structure at 4.0A) rather than UniProt Ig-core
annotations. Signal peptide (residues 1-26) is always excluded -- all published
structural studies use the mature protein starting at Ser27.

Transmembrane protein handling:
  sFLT1 is soluble/extracellular, so it can only interact with the extracellular
  portion of transmembrane partners. For TM proteins, we extract the full
  ectodomain (after signal peptide, before first TM helix) from UniProt topology
  annotations. Cytoplasmic domains are excluded to avoid biologically impossible
  interface predictions and reduce GPU cost.

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
# Signal peptide: residues 1-26 (cleaved, never included)
# Mature protein starts at Ser27
#
# Domain boundaries from 5T89 crystal structure (structural, includes linkers):
#   D1: 32-130    D2: 132-225    D3: 226-330
#   D4: 333-425   D5: 426-555    D6: 556-657
#   D7: 661-747   (UniProt; not resolved in 5T89)
#

SFLT1_UNIPROT = "P17948"


@dataclass(frozen=True)
class FLT1Construct:
    """sFLT1/FLT1 ectodomain construct definition."""
    name: str           # e.g. "D1-D3", "D1-D6", "D1-D7"
    start: int          # P17948 precursor residue (1-indexed, inclusive)
    end: int            # P17948 precursor residue (1-indexed, inclusive)
    description: str    # for FASTA header

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def fasta_header_tag(self) -> str:
        return f"sFLT1_{self.name}"


# Canonical constructs
CONSTRUCT_D1D3 = FLT1Construct(
    name="D1-D3",
    start=27,
    end=330,
    description="Mature sFLT1 D1-D3, structural boundaries (5T89)",
)

CONSTRUCT_D1D6 = FLT1Construct(
    name="D1-D6",
    start=27,
    end=657,
    description="Mature sFLT1 D1-D6, structural boundaries (5T89 D1-D6 + linker margin)",
)

CONSTRUCT_D1D7 = FLT1Construct(
    name="D1-D7",
    start=27,
    end=747,
    description="FLT1 full ectodomain D1-D7, UniProt D7 end",
)

CONSTRUCTS = {
    "d1d3": CONSTRUCT_D1D3,
    "d1d6": CONSTRUCT_D1D6,
    "d1d7": CONSTRUCT_D1D7,
}

# Backwards compat
SFLT1_D1D3_START = CONSTRUCT_D1D3.start
SFLT1_D1D3_END = CONSTRUCT_D1D3.end

# UniProt REST API
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"

# Rate limiting
REQUEST_DELAY = 0.5  # seconds between requests
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds


def fetch_fasta_from_uniprot(accession: str) -> str | None:
    """Fetch a protein FASTA sequence from UniProt REST API.

    Parameters
    ----------
    accession : str
        UniProt accession (e.g., "P17948").

    Returns
    -------
    str or None
        Raw FASTA text, or None on failure.
    """
    url = UNIPROT_API_URL.format(accession=accession)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response.text
            if response.status_code == 429:
                logger.warning(
                    "Rate limited for %s, waiting %ds (attempt %d/%d)",
                    accession, RETRY_DELAY, attempt, MAX_RETRIES,
                )
                time.sleep(RETRY_DELAY)
                continue
            logger.warning(
                "UniProt returned %d for %s (attempt %d/%d)",
                response.status_code, accession, attempt, MAX_RETRIES,
            )
        except requests.RequestException as e:
            logger.warning(
                "Request failed for %s: %s (attempt %d/%d)",
                accession, e, attempt, MAX_RETRIES,
            )
        time.sleep(RETRY_DELAY)

    logger.error("Failed to fetch %s after %d attempts", accession, MAX_RETRIES)
    return None


def parse_fasta_sequence(fasta_text: str) -> str:
    """Extract the sequence from FASTA text (strip header and whitespace).

    Parameters
    ----------
    fasta_text : str
        Raw FASTA text.

    Returns
    -------
    str
        Amino acid sequence.
    """
    lines = fasta_text.strip().split("\n")
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq_lines)


@dataclass(frozen=True)
class ProteinTopology:
    """Protein topology from UniProt annotations."""
    accession: str
    seq_length: int
    signal_peptide: tuple[int, int] | None  # (start, end) 1-indexed
    transmembrane: tuple[tuple[int, int], ...] = field(default_factory=tuple)
    protein_type: str = "soluble"  # soluble, type_i_tm, multi_tm, gpi_anchored

    @property
    def is_transmembrane(self) -> bool:
        return len(self.transmembrane) > 0

    @property
    def ectodomain_range(self) -> tuple[int, int] | None:
        """Return (start, end) 1-indexed of the full extracellular domain.

        For type I TM proteins: after signal peptide, before first TM helix.
        For multi-TM: after signal peptide, before first TM helix (N-terminal ecto).
        For soluble/secreted: None (use full mature protein).
        """
        if not self.is_transmembrane:
            return None
        ecto_start = (self.signal_peptide[1] + 1) if self.signal_peptide else 1
        ecto_end = self.transmembrane[0][0] - 1
        if ecto_end < ecto_start:
            return None
        return (ecto_start, ecto_end)

    @property
    def ectodomain_length(self) -> int | None:
        r = self.ectodomain_range
        if r is None:
            return None
        return r[1] - r[0] + 1


# UniProt JSON API for topology
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

# Minimum ectodomain length to model (skip if shorter)
MIN_ECTODOMAIN_LENGTH = 50


def fetch_topology_from_uniprot(accession: str) -> ProteinTopology | None:
    """Fetch protein topology annotations from UniProt JSON API.

    Parameters
    ----------
    accession : str
        UniProt accession (e.g., "O14786").

    Returns
    -------
    ProteinTopology or None
        Topology data, or None on failure.
    """
    url = UNIPROT_JSON_URL.format(accession=accession)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                data = response.json()
                return _parse_topology(accession, data)
            if response.status_code == 429:
                logger.warning(
                    "Rate limited for %s topology, waiting %ds (attempt %d/%d)",
                    accession, RETRY_DELAY, attempt, MAX_RETRIES,
                )
                time.sleep(RETRY_DELAY)
                continue
            logger.warning(
                "UniProt JSON returned %d for %s (attempt %d/%d)",
                response.status_code, accession, attempt, MAX_RETRIES,
            )
        except requests.RequestException as e:
            logger.warning(
                "Request failed for %s topology: %s (attempt %d/%d)",
                accession, e, attempt, MAX_RETRIES,
            )
        time.sleep(RETRY_DELAY)

    logger.error("Failed to fetch topology for %s after %d attempts", accession, MAX_RETRIES)
    return None


def _parse_topology(accession: str, data: dict) -> ProteinTopology:
    """Parse UniProt JSON into ProteinTopology."""
    seq_length = data.get("sequence", {}).get("length", 0)
    features = data.get("features", [])

    signal_peptide = None
    transmembrane = []

    for feat in features:
        ftype = feat.get("type", "")
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")

        if start is None or end is None:
            continue

        if ftype == "Signal":
            signal_peptide = (start, end)
        elif ftype == "Transmembrane":
            transmembrane.append((start, end))

    # Classify protein type
    if len(transmembrane) > 1:
        protein_type = "multi_tm"
    elif len(transmembrane) == 1:
        protein_type = "type_i_tm"
    else:
        protein_type = "soluble"

    return ProteinTopology(
        accession=accession,
        seq_length=seq_length,
        signal_peptide=signal_peptide,
        transmembrane=tuple(transmembrane),
        protein_type=protein_type,
    )


def trim_sequence(sequence: str, start: int, end: int) -> str:
    """Trim a protein sequence to a specific residue range (1-indexed).

    Parameters
    ----------
    sequence : str
        Full amino acid sequence.
    start : int
        Start residue (1-indexed, inclusive).
    end : int
        End residue (1-indexed, inclusive).

    Returns
    -------
    str
        Trimmed sequence.
    """
    return sequence[start - 1 : end]


def get_sflt1_sequence(
    construct: FLT1Construct | None = None,
) -> str | None:
    """Fetch and trim the sFLT1/FLT1 ectodomain sequence.

    Parameters
    ----------
    construct : FLT1Construct, optional
        Which construct to extract. Defaults to CONSTRUCT_D1D3.

    Returns
    -------
    str or None
        Trimmed sequence, or None on failure.
    """
    if construct is None:
        construct = CONSTRUCT_D1D3

    fasta_text = fetch_fasta_from_uniprot(SFLT1_UNIPROT)
    if fasta_text is None:
        return None

    full_seq = parse_fasta_sequence(fasta_text)
    trimmed = trim_sequence(full_seq, construct.start, construct.end)
    logger.info(
        "sFLT1 %s: %d residues (trimmed from %d, residues %d-%d)",
        construct.name, len(trimmed), len(full_seq),
        construct.start, construct.end,
    )
    return trimmed


def write_two_chain_fasta(
    sflt1_seq: str,
    partner_seq: str,
    target_name: str,
    partner_uniprot: str,
    output_path: Path,
    construct: FLT1Construct | None = None,
    partner_region: str = "",
) -> None:
    """Write a two-chain FASTA file for AlphaFold Multimer.

    Parameters
    ----------
    sflt1_seq : str
        sFLT1/FLT1 ectodomain amino acid sequence.
    partner_seq : str
        Partner protein amino acid sequence.
    target_name : str
        Short name for the partner protein.
    partner_uniprot : str
        UniProt accession of the partner.
    output_path : Path
        Output FASTA file path.
    construct : FLT1Construct, optional
        Construct used for Chain A header. Defaults to CONSTRUCT_D1D3.
    partner_region : str
        Description of the partner region (e.g., "ectodomain_22-856").
        Appended to the partner FASTA header for provenance.
    """
    if construct is None:
        construct = CONSTRUCT_D1D3

    partner_header = f">{target_name}|{partner_uniprot}"
    if partner_region:
        partner_header += f"|{partner_region}"

    with open(output_path, "w") as f:
        f.write(
            f">{construct.fasta_header_tag}|{SFLT1_UNIPROT}"
            f"|residues_{construct.start}-{construct.end}\n"
        )
        for i in range(0, len(sflt1_seq), 80):
            f.write(sflt1_seq[i : i + 80] + "\n")
        f.write(partner_header + "\n")
        for i in range(0, len(partner_seq), 80):
            f.write(partner_seq[i : i + 80] + "\n")


def fetch_all_sequences(
    candidates_path: Path,
    fasta_dir: Path,
    construct: FLT1Construct | None = None,
    apply_ectodomain_filter: bool = True,
) -> dict[str, dict]:
    """Fetch sequences for all candidates and write two-chain FASTAs.

    For transmembrane proteins, fetches UniProt topology and truncates to
    the extracellular domain (after signal peptide, before first TM helix).
    This is required because sFLT1 is extracellular and cannot interact with
    cytoplasmic domains. See CLAUDE.md "Protein localization filter" guardrail.

    Parameters
    ----------
    candidates_path : Path
        Path to candidates CSV (must have 'target' and 'uniprot' columns).
    fasta_dir : Path
        Output directory for FASTA files.
    construct : FLT1Construct, optional
        Which sFLT1/FLT1 construct to use for Chain A. Defaults to D1-D3.
    apply_ectodomain_filter : bool
        If True (default), truncate transmembrane proteins to their
        extracellular domain. If False, use full-length sequences.

    Returns
    -------
    dict[str, dict]
        Mapping of target name -> {"uniprot": str, "fasta_path": str,
        "partner_length": int, "total_residues": int, "status": str,
        "topology": str, "region": str}.
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
    sflt1_ref_path = fasta_dir / f"sflt1_{tag}_reference.fasta"
    with open(sflt1_ref_path, "w") as f:
        f.write(
            f">{construct.fasta_header_tag}|{SFLT1_UNIPROT}"
            f"|residues_{construct.start}-{construct.end}\n"
        )
        for i in range(0, len(sflt1_seq), 80):
            f.write(sflt1_seq[i : i + 80] + "\n")
    logger.info(
        "Saved sFLT1 reference: %s (%d aa)", sflt1_ref_path.name, len(sflt1_seq)
    )

    results = {}
    for _, row in candidates.iterrows():
        target = row["target"]
        uniprot = row["uniprot"]

        if not isinstance(uniprot, str) or len(uniprot) < 3:
            logger.warning("Skipping %s: invalid UniProt accession '%s'", target, uniprot)
            results[target] = {
                "uniprot": uniprot, "fasta_path": None,
                "partner_length": 0, "total_residues": 0,
                "status": "no_uniprot", "topology": "unknown", "region": "",
            }
            continue

        logger.info("Fetching %s (%s)...", target, uniprot)
        fasta_text = fetch_fasta_from_uniprot(uniprot)
        time.sleep(REQUEST_DELAY)

        if fasta_text is None:
            results[target] = {
                "uniprot": uniprot, "fasta_path": None,
                "partner_length": 0, "total_residues": 0,
                "status": "fetch_failed", "topology": "unknown", "region": "",
            }
            continue

        full_seq = parse_fasta_sequence(fasta_text)
        partner_seq = full_seq
        partner_region = ""
        topology_type = "soluble"

        # Apply ectodomain filter for transmembrane proteins
        if apply_ectodomain_filter:
            topo = fetch_topology_from_uniprot(uniprot)
            time.sleep(REQUEST_DELAY)

            if topo is not None and topo.is_transmembrane:
                topology_type = topo.protein_type
                ecto_range = topo.ectodomain_range

                if ecto_range is not None:
                    ecto_len = ecto_range[1] - ecto_range[0] + 1
                    if ecto_len < MIN_ECTODOMAIN_LENGTH:
                        logger.warning(
                            "  %s: ectodomain too short (%d aa < %d min), "
                            "using full-length",
                            target, ecto_len, MIN_ECTODOMAIN_LENGTH,
                        )
                    else:
                        partner_seq = trim_sequence(
                            full_seq, ecto_range[0], ecto_range[1]
                        )
                        partner_region = (
                            f"ectodomain_{ecto_range[0]}-{ecto_range[1]}"
                        )
                        logger.info(
                            "  %s: %s, ectodomain %d-%d (%d aa, "
                            "trimmed from %d aa, saved %d aa)",
                            target, topology_type,
                            ecto_range[0], ecto_range[1], len(partner_seq),
                            len(full_seq), len(full_seq) - len(partner_seq),
                        )
                else:
                    logger.warning(
                        "  %s: TM protein but no ectodomain range found, "
                        "using full-length",
                        target,
                    )
            elif topo is not None:
                topology_type = topo.protein_type

        safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")
        fasta_path = fasta_dir / f"sflt1_vs_{safe_name}.fasta"

        write_two_chain_fasta(
            sflt1_seq, partner_seq, target, uniprot, fasta_path,
            construct=construct,
            partner_region=partner_region,
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
            "region": partner_region if partner_region else "full_length",
        }
        if partner_region:
            logger.info(
                "  -> %s: %d aa ectodomain (of %d full), %d total, saved %s",
                target, len(partner_seq), len(full_seq), total, fasta_path.name,
            )
        else:
            logger.info(
                "  -> %s: %d aa (full-length), %d total, saved %s",
                target, len(partner_seq), total, fasta_path.name,
            )

    return results
