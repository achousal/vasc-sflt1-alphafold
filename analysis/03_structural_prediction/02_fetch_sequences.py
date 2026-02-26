"""02_fetch_sequences.py -- Fetch FASTA sequences from UniProt REST API.

Produces two-chain FASTA files for AlphaFold Multimer:
  Chain A: sFLT1 D1-D3 (residues 1-338 of FLT1/P17948)
  Chain B: candidate protein
"""

import logging
import time
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

# sFLT1 sequence: FLT1 (P17948) Ig-like domains 1-3, residues 1-338
# This is the soluble form (sVEGFR1/sFLT1) extracellular domain D1-D3.
SFLT1_UNIPROT = "P17948"
SFLT1_D1D3_START = 1
SFLT1_D1D3_END = 338

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


def get_sflt1_sequence() -> str | None:
    """Fetch and trim the sFLT1 D1-D3 sequence.

    Returns
    -------
    str or None
        sFLT1 D1-D3 sequence (338 residues), or None on failure.
    """
    fasta_text = fetch_fasta_from_uniprot(SFLT1_UNIPROT)
    if fasta_text is None:
        return None

    full_seq = parse_fasta_sequence(fasta_text)
    trimmed = trim_sequence(full_seq, SFLT1_D1D3_START, SFLT1_D1D3_END)
    logger.info(
        "sFLT1 D1-D3: %d residues (trimmed from %d)", len(trimmed), len(full_seq)
    )
    return trimmed


def write_two_chain_fasta(
    sflt1_seq: str,
    partner_seq: str,
    target_name: str,
    partner_uniprot: str,
    output_path: Path,
) -> None:
    """Write a two-chain FASTA file for AlphaFold Multimer.

    Parameters
    ----------
    sflt1_seq : str
        sFLT1 D1-D3 amino acid sequence.
    partner_seq : str
        Partner protein amino acid sequence.
    target_name : str
        Short name for the partner protein.
    partner_uniprot : str
        UniProt accession of the partner.
    output_path : Path
        Output FASTA file path.
    """
    with open(output_path, "w") as f:
        f.write(f">sFLT1_D1-D3|{SFLT1_UNIPROT}|residues_{SFLT1_D1D3_START}-{SFLT1_D1D3_END}\n")
        # Write sequence in 80-char lines
        for i in range(0, len(sflt1_seq), 80):
            f.write(sflt1_seq[i : i + 80] + "\n")
        f.write(f">{target_name}|{partner_uniprot}\n")
        for i in range(0, len(partner_seq), 80):
            f.write(partner_seq[i : i + 80] + "\n")


def fetch_all_sequences(
    candidates_path: Path,
    fasta_dir: Path,
) -> dict[str, dict]:
    """Fetch sequences for all candidates and write two-chain FASTAs.

    Parameters
    ----------
    candidates_path : Path
        Path to step03_candidates.csv.
    fasta_dir : Path
        Output directory for FASTA files.

    Returns
    -------
    dict[str, dict]
        Mapping of target name -> {"uniprot": str, "fasta_path": str,
        "partner_length": int, "status": str}.
    """
    import pandas as pd

    candidates = pd.read_csv(candidates_path)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    # Fetch sFLT1 first
    logger.info("Fetching sFLT1 D1-D3 sequence...")
    sflt1_seq = get_sflt1_sequence()
    if sflt1_seq is None:
        logger.error("Failed to fetch sFLT1 sequence. Aborting.")
        return {}

    # Save sFLT1 reference
    sflt1_ref_path = fasta_dir / "sflt1_d1d3_reference.fasta"
    with open(sflt1_ref_path, "w") as f:
        f.write(
            f">sFLT1_D1-D3|{SFLT1_UNIPROT}|residues_{SFLT1_D1D3_START}-{SFLT1_D1D3_END}\n"
        )
        for i in range(0, len(sflt1_seq), 80):
            f.write(sflt1_seq[i : i + 80] + "\n")
    logger.info("Saved sFLT1 reference: %s (%d aa)", sflt1_ref_path.name, len(sflt1_seq))

    results = {}
    for _, row in candidates.iterrows():
        target = row["target"]
        uniprot = row["uniprot"]

        if not isinstance(uniprot, str) or len(uniprot) < 3:
            logger.warning("Skipping %s: invalid UniProt accession '%s'", target, uniprot)
            results[target] = {
                "uniprot": uniprot, "fasta_path": None,
                "partner_length": 0, "status": "no_uniprot",
            }
            continue

        logger.info("Fetching %s (%s)...", target, uniprot)
        fasta_text = fetch_fasta_from_uniprot(uniprot)
        time.sleep(REQUEST_DELAY)

        if fasta_text is None:
            results[target] = {
                "uniprot": uniprot, "fasta_path": None,
                "partner_length": 0, "status": "fetch_failed",
            }
            continue

        partner_seq = parse_fasta_sequence(fasta_text)
        # Sanitize filename
        safe_name = target.replace("/", "-").replace(" ", "_").replace(":", "_")
        fasta_path = fasta_dir / f"sflt1_vs_{safe_name}.fasta"

        write_two_chain_fasta(sflt1_seq, partner_seq, target, uniprot, fasta_path)

        results[target] = {
            "uniprot": uniprot,
            "fasta_path": str(fasta_path),
            "partner_length": len(partner_seq),
            "status": "ok",
        }
        logger.info("  -> %s: %d aa, saved %s", target, len(partner_seq), fasta_path.name)

    return results
