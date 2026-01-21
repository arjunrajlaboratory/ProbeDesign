"""FASTA file I/O utilities."""

from typing import List, Tuple
from pathlib import Path


def read_fasta(filepath: str) -> Tuple[List[str], List[str]]:
    """Read a FASTA file and return headers and sequences.

    Args:
        filepath: Path to FASTA file

    Returns:
        Tuple of (headers, sequences) where:
        - headers: List of header lines (without leading '>')
        - sequences: List of corresponding sequences
    """
    headers: List[str] = []
    sequences: List[str] = []
    current_seq: List[str] = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Save previous sequence if exists
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
                # Store header without '>'
                headers.append(line[1:])
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_seq:
            sequences.append(''.join(current_seq))

    return headers, sequences


def sequences_to_single_string(seqs: List[str], mark_junctions: bool = True) -> str:
    """Concatenate multiple sequences into a single string.

    Used when a FASTA file has multiple entries (e.g., exons)
    that should be treated as one continuous sequence.

    Args:
        seqs: List of sequences
        mark_junctions: If True and there are multiple sequences,
            prefix each with '>' to mark exon junctions (MATLAB behavior).
            If there's only one sequence, just prefix with '>' for the header.

    Returns:
        Single concatenated sequence (with '>' at junctions if mark_junctions=True)
    """
    if not seqs:
        return ''

    if mark_junctions:
        # Include '>' at the start and at each exon junction
        return '>' + '>'.join(seqs)
    else:
        return ''.join(seqs)


def write_fasta(filepath: str, headers: List[str], sequences: List[str],
                line_width: int = 60) -> None:
    """Write sequences to a FASTA file.

    Args:
        filepath: Output file path
        headers: List of sequence headers (without '>')
        sequences: List of sequences
        line_width: Number of characters per line (default 60)
    """
    with open(filepath, 'w') as f:
        for header, seq in zip(headers, sequences):
            f.write(f'>{header}\n')
            # Write sequence in chunks of line_width
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + '\n')


def get_basename(filepath: str) -> str:
    """Extract the base name from a file path (without directory and extension).

    Args:
        filepath: File path

    Returns:
        Base name without extension
    """
    return Path(filepath).stem
