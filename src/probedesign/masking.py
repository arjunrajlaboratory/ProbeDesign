"""Sequence masking using bowtie alignment.

This module provides functions to mask regions of a sequence that align to:
- Pseudogene databases (to avoid cross-hybridization)
- Genome databases (to avoid repetitive regions)

The masking logic mirrors the MATLAB findprobesLocal.m implementation.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

# Default index directory (can be overridden)
DEFAULT_INDEX_DIR = Path(__file__).parent.parent.parent.parent / "bowtie_indexes"


def find_bowtie() -> str:
    """Find the bowtie executable.

    Searches in order:
    1. conda/mamba environment (if active)
    2. /usr/bin/bowtie
    3. BOWTIEHOME environment variable

    Returns:
        Path to bowtie executable

    Raises:
        FileNotFoundError: If bowtie cannot be found
    """
    # Check if in PATH (e.g., conda environment)
    try:
        result = subprocess.run(
            ["which", "bowtie"],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        pass

    # Check common locations
    locations = [
        "/usr/bin/bowtie",
        os.path.expandvars("$BOWTIEHOME/bowtie"),
        os.path.expanduser("~/bowtie/bowtie"),
    ]

    for loc in locations:
        if os.path.exists(loc):
            return loc

    raise FileNotFoundError(
        "Could not find bowtie executable. "
        "Please install bowtie and ensure it is in your PATH, "
        "or set the BOWTIEHOME environment variable."
    )


def sequence_to_nmers(seq: str, mer_length: int) -> str:
    """Convert a sequence to FASTA-formatted n-mer substrings.

    Each n-mer is given a numeric ID corresponding to its start position.

    Args:
        seq: Input sequence (lowercase)
        mer_length: Length of each substring

    Returns:
        FASTA-formatted string with each n-mer as a separate entry
    """
    lines = []
    for i in range(len(seq) - mer_length + 1):
        nmer = seq[i:i + mer_length]
        lines.append(f">{i}")
        lines.append(nmer)
    return "\n".join(lines)


def run_bowtie(
    seq: str,
    mer_length: int,
    database: str,
    index_dir: Optional[Path] = None,
    mismatches: int = 0,
    max_hits: int = 1,
) -> List[int]:
    """Run bowtie alignment and return hit counts per position.

    Args:
        seq: Input sequence (lowercase)
        mer_length: Length of n-mers to align
        database: Name of the bowtie index (e.g., 'humanPseudo', 'human')
        index_dir: Directory containing bowtie indexes
        mismatches: Number of allowed mismatches (default 0)
        max_hits: Maximum alignments to report per read (default 1)

    Returns:
        List of hit counts, one per position in the sequence
    """
    if index_dir is None:
        index_dir = DEFAULT_INDEX_DIR

    # Initialize hit counts
    num_positions = len(seq) - mer_length + 1
    hits = [0] * len(seq)

    if num_positions <= 0:
        return hits

    # Generate n-mer FASTA
    nmers_fasta = sequence_to_nmers(seq, mer_length)

    # Find bowtie
    bowtie_path = find_bowtie()

    # Build index path
    index_path = index_dir / database

    # Run bowtie
    # -f: input is FASTA
    # -v: allowed mismatches
    # -k: max alignments to report
    # --quiet: suppress verbose output
    cmd = [
        bowtie_path,
        "-f",
        "-v", str(mismatches),
        "-k", str(max_hits),
        "--quiet",
        str(index_path),
        "-"  # read from stdin
    ]

    try:
        result = subprocess.run(
            cmd,
            input=nmers_fasta,
            capture_output=True,
            text=True,
            check=False  # Don't raise on non-zero exit (no alignments is OK)
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Bowtie not found: {e}")

    # Parse bowtie output
    # Default output format: read_name, strand, ref_name, offset, seq, qual, num_other_alignments, mismatches
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) >= 7:
            try:
                position = int(cols[0])  # Read name is the position
                num_other = int(cols[6])  # Number of other alignments
                hits[position] = num_other + 1  # Total hits = other + 1
            except (ValueError, IndexError):
                continue

    return hits


def hits_to_mask(
    hits: List[int],
    mer_length: int,
    threshold: int = 0,
) -> List[int]:
    """Convert hit counts to a binary mask.

    Positions where hits > threshold are marked as 1, and the mask is
    extended to cover the full mer_length from each hit position.

    Args:
        hits: Hit counts per position
        mer_length: Length of the n-mers used for alignment
        threshold: Minimum hits to trigger masking (default 0 = any hit)

    Returns:
        Binary mask (0 or 1 for each position)
    """
    mask = [0] * len(hits)

    for i, hit_count in enumerate(hits):
        if hit_count > threshold:
            # Extend mask to cover the full mer
            for j in range(i, min(i + mer_length, len(mask))):
                mask[j] = 1

    return mask


def remove_short_runs(
    mask: List[int],
    min_length: int = 20,
    tolerance: int = 2,
) -> List[int]:
    """Remove short masked runs from a mask.

    Keeps only runs of masked positions that are at least (min_length - tolerance)
    long. This helps avoid masking very short regions that are likely spurious.

    Args:
        mask: Binary mask
        min_length: Minimum run length to keep (default 20)
        tolerance: Allow this many gaps in the run (default 2)

    Returns:
        Filtered mask
    """
    out = [0] * len(mask)

    for i in range(len(mask)):
        if mask[i] == 1:
            # Check if this starts a sufficiently long run
            end = min(i + min_length, len(mask))
            window = mask[i:end]
            if sum(window) >= min_length - tolerance:
                # Keep this run
                for j in range(i, end):
                    if j < len(mask):
                        out[j] = mask[j]

    return out


def pseudogene_mask(
    seq: str,
    species: str = "human",
    index_dir: Optional[Path] = None,
) -> List[int]:
    """Create a mask for regions that align to pseudogenes.

    Uses 16-mer alignment with threshold of 0 (any hit triggers masking).
    Short runs (<18bp) are removed to avoid spurious masking.

    Args:
        seq: Input sequence (lowercase)
        species: Species for pseudogene database ('human', 'mouse', etc.)
        index_dir: Directory containing bowtie indexes

    Returns:
        Binary mask (1 = masked, 0 = unmasked)
    """
    # Map species to database name
    db_map = {
        "human": "humanPseudo",
        "mouse": "mousePseudo",
        "elegans": "celegansPseudo",
        "drosophila": "drosophilaPseudo",
        "rat": "ratPseudo",
    }

    db = db_map.get(species.lower())
    if not db:
        raise ValueError(f"Unknown species: {species}. Valid options: {list(db_map.keys())}")

    # Run bowtie with 16-mers
    hits = run_bowtie(seq, mer_length=16, database=db, index_dir=index_dir)

    # Convert to mask (threshold=0 means any hit)
    mask = hits_to_mask(hits, mer_length=16, threshold=0)

    # Remove short runs
    mask = remove_short_runs(mask, min_length=20, tolerance=2)

    return mask


def genome_mask(
    seq: str,
    species: str = "human",
    index_dir: Optional[Path] = None,
) -> List[int]:
    """Create a mask for regions that align to multiple genomic locations.

    Uses multiple mer lengths with different thresholds:
    - 12-mer with threshold 4000 (very repetitive)
    - 14-mer with threshold 500
    - 16-mer with threshold 20

    Args:
        seq: Input sequence (lowercase)
        species: Species for genome database ('human', 'mouse', etc.)
        index_dir: Directory containing bowtie indexes

    Returns:
        Binary mask (1 = masked, 0 = unmasked)
    """
    # Map species to database name
    db_map = {
        "human": "GCA_000001405.15_GRCh38_no_alt_analysis_set",
        "mouse": "mm10",
        "elegans": "celegans",
        "drosophila": "drosophila",
        "rat": "rat",
        "cow": "cow",
    }

    db = db_map.get(species.lower())
    if not db:
        raise ValueError(f"Unknown species: {species}. Valid options: {list(db_map.keys())}")

    # Run bowtie with different parameters
    # Note: We need -k to be high enough to count hits above threshold
    hits_12 = run_bowtie(seq, mer_length=12, database=db, index_dir=index_dir, max_hits=5000)
    hits_14 = run_bowtie(seq, mer_length=14, database=db, index_dir=index_dir, max_hits=1000)
    hits_16 = run_bowtie(seq, mer_length=16, database=db, index_dir=index_dir, max_hits=100)

    # Convert to masks with respective thresholds
    mask_12 = hits_to_mask(hits_12, mer_length=12, threshold=4000)
    mask_14 = hits_to_mask(hits_14, mer_length=14, threshold=500)
    mask_16 = hits_to_mask(hits_16, mer_length=16, threshold=20)

    # Combine masks (OR)
    combined = [0] * len(seq)
    for i in range(len(seq)):
        if mask_12[i] or mask_14[i] or mask_16[i]:
            combined[i] = 1

    return combined


def mask_to_badness(
    mask: List[int],
    oligo_length: int,
) -> List[float]:
    """Convert a mask to badness scores for probe design.

    An oligo that overlaps any masked position gets infinite badness.

    Args:
        mask: Binary mask (1 = masked)
        oligo_length: Length of oligonucleotides

    Returns:
        List of badness scores (0 or inf) for each oligo start position
    """
    num_oligos = len(mask) - oligo_length + 1
    badness = [0.0] * num_oligos

    for i in range(num_oligos):
        # Check if any position in this oligo is masked
        if any(mask[i:i + oligo_length]):
            badness[i] = float('inf')

    return badness


def create_full_mask(
    seq: str,
    species: str = "human",
    pseudogene: bool = True,
    genome: bool = True,
    index_dir: Optional[Path] = None,
) -> Tuple[List[int], List[str]]:
    """Create a combined mask from all masking sources.

    Args:
        seq: Input sequence (lowercase)
        species: Species for databases
        pseudogene: Whether to use pseudogene masking
        genome: Whether to use genome masking
        index_dir: Directory containing bowtie indexes

    Returns:
        Tuple of (combined_mask, mask_strings) where mask_strings contains
        visualization strings for each mask type
    """
    combined = [0] * len(seq)
    mask_strings = []

    if pseudogene:
        try:
            pmask = pseudogene_mask(seq, species, index_dir)
            for i, v in enumerate(pmask):
                combined[i] += v
            # Create visualization string
            pstr = "".join("P" if pmask[i] else seq[i] for i in range(len(seq)))
            mask_strings.append(pstr)
        except (FileNotFoundError, subprocess.CalledProcessError) as e:
            print(f"Warning: Pseudogene masking failed: {e}")

    if genome:
        try:
            gmask = genome_mask(seq, species, index_dir)
            for i, v in enumerate(gmask):
                combined[i] += v
            # Create visualization string
            gstr = "".join("B" if gmask[i] else seq[i] for i in range(len(seq)))
            mask_strings.append(gstr)
        except (FileNotFoundError, subprocess.CalledProcessError) as e:
            print(f"Warning: Genome masking failed: {e}")

    return combined, mask_strings
