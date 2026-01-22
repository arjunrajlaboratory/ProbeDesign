"""Sequence manipulation utilities."""

import re
from typing import Set

# Valid nucleotide characters
VALID_BASES: Set[str] = {'a', 'c', 'g', 't', 'A', 'C', 'G', 'T'}

# Complement mapping
COMPLEMENT_MAP = str.maketrans('acgtACGT', 'tgcaTGCA')


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence (any case)

    Returns:
        Reverse complement sequence in lowercase
    """
    return seq.lower().translate(COMPLEMENT_MAP)[::-1]


def complement(seq: str) -> str:
    """Return the complement of a DNA sequence.

    Args:
        seq: DNA sequence (any case)

    Returns:
        Complement sequence in lowercase
    """
    return seq.lower().translate(COMPLEMENT_MAP)


def percent_gc(seq: str) -> float:
    """Calculate the GC content of a sequence.

    Args:
        seq: DNA sequence (any case)

    Returns:
        GC content as a fraction (0.0 to 1.0)
    """
    seq = seq.lower()
    gc_count = seq.count('g') + seq.count('c')
    return gc_count / len(seq) if seq else 0.0


def is_valid_sequence(seq: str) -> bool:
    """Check if sequence contains only valid ACGT bases.

    Args:
        seq: DNA sequence to check

    Returns:
        True if sequence contains only A, C, G, T (case insensitive)
    """
    return bool(re.fullmatch(r'[aAcCgGtT]+', seq))


def has_invalid_chars(seq: str) -> bool:
    """Check if sequence contains invalid characters (masked regions, etc.).

    Characters that indicate masked/invalid regions:
    - x, X: General mask
    - n, N: Unknown nucleotide
    - m, M, h, H, p, P, b, B: Various mask types used in MATLAB
    - >: Sequence boundary marker

    Args:
        seq: DNA sequence to check

    Returns:
        True if sequence contains any invalid/masked characters
    """
    invalid_pattern = re.compile(r'[xXmMhHpPbBnN>]')
    return bool(invalid_pattern.search(seq))


def clean_sequence(seq: str, keep_chars: str = 'actgxn>') -> str:
    """Clean a sequence by keeping only specified characters.

    Converts to lowercase and removes any character not in keep_chars.

    Args:
        seq: Input sequence
        keep_chars: String of characters to keep (lowercase)

    Returns:
        Cleaned sequence in lowercase
    """
    seq = seq.lower()
    return ''.join(c for c in seq if c in keep_chars)
