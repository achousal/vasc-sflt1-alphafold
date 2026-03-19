"""02_fetch_sequences.py -- Fetch FASTA sequences from UniProt REST API.

Produces two-chain FASTA files for AlphaFold Multimer:
  Chain A: sFLT1 construct (D1-D3, D1-D6, or D1-D7)
  Chain B: candidate protein

Construct definitions use structural boundaries from 5T89 (Markovic-Mueller
2017, VEGF + FLT1 D1-D6 crystal structure at 4.0A) rather than UniProt Ig-core
annotations. Signal peptide (residues 1-26) is always excluded -- all published
structural studies use the mature protein starting at Ser27.

Numbering: all residue numbers refer to P17948 (FLT1_HUMAN) precursor.
"""

import logging
import time
from dataclasses import dataclass
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
# Previous (incorrect) construct used residues 1-338, which included the
# 26-aa signal peptide and 8 residues of D4 linker. Fixed 2026-03-19.

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
    """
    if construct is None:
        construct = CONSTRUCT_D1D3

    with open(output_path, "w") as f:
        f.write(
            f">{construct.fasta_header_tag}|{SFLT1_UNIPROT}"
            f"|residues_{construct.start}-{construct.end}\n"
        )
        for i in range(0, len(sflt1_seq), 80):
            f.write(sflt1_seq[i : i + 80] + "\n")
        f.write(f">{target_name}|{partner_uniprot}\n")
        for i in range(0, len(partner_seq), 80):
            f.write(partner_seq[i : i + 80] + "\n")


def fetch_all_sequences(
    candidates_path: Path,
    fasta_dir: Path,
    construct: FLT1Construct | None = None,
) -> dict[str, dict]:
    """Fetch sequences for all candidates and write two-chain FASTAs.

    Parameters
    ----------
    candidates_path : Path
        Path to candidates CSV (must have 'target' and 'uniprot' columns).
    fasta_dir : Path
        Output directory for FASTA files.
    construct : FLT1Construct, optional
        Which sFLT1/FLT1 construct to use for Chain A. Defaults to D1-D3.

    Returns
    -------
    dict[str, dict]
        Mapping of target name -> {"uniprot": str, "fasta_path": str,
        "partner_length": int, "total_residues": int, "status": str}.
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
                "partner_length": 0, "total_residues": 0, "status": "no_uniprot",
            }
            continue

        logger.info("Fetching %s (%s)...", target, uniprot)
        fasta_text = fetch_fasta_from_uniprot(uniprot)
        time.sleep(REQUEST_DELAY)

        if fasta_text is None:
            results[target] = {
                "uniprot": uniprot, "fasta_path": None,
                "partner_length": 0, "total_residues": 0, "status": "fetch_failed",
            }
            continue

        partner_seq = parse_fasta_sequence(fasta_text)
        safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")
        fasta_path = fasta_dir / f"sflt1_vs_{safe_name}.fasta"

        write_two_chain_fasta(
            sflt1_seq, partner_seq, target, uniprot, fasta_path,
            construct=construct,
        )

        total = construct.length + len(partner_seq)
        results[target] = {
            "uniprot": uniprot,
            "fasta_path": str(fasta_path),
            "partner_length": len(partner_seq),
            "total_residues": total,
            "status": "ok",
        }
        logger.info(
            "  -> %s: %d aa partner, %d total residues, saved %s",
            target, len(partner_seq), total, fasta_path.name,
        )

    return results
