"""Output file generation for probe designs."""

from typing import List, Optional
from .core import Probe, ProbeDesignResult
from .sequence import complement


def write_oligos_file(result: ProbeDesignResult, filepath: str) -> None:
    """Write probe information to a TSV file.

    Output format (tab-separated):
    index  GC%  Tm  GibbsFE  sequence  name

    Args:
        result: ProbeDesignResult from design_probes()
        filepath: Output file path
    """
    with open(filepath, 'w') as f:
        for probe in result.probes:
            f.write(
                f"{probe.index}\t"
                f"{probe.gc_percent}\t"
                f"{probe.tm}\t"
                f"{probe.gibbs_fe}\t"
                f"{probe.sequence}\t"
                f"{probe.name}\n"
            )


def write_seq_file(
    result: ProbeDesignResult,
    filepath: str,
    mask_seqs: Optional[List[str]] = None,
    line_width: int = 110,
) -> None:
    """Write sequence alignment visualization file.

    Creates a multi-line visualization showing:
    - Original sequence (with '>' at exon junctions where present)
    - Masked regions (if any) - each mask shows sequence with mask chars (P, B, R, F)
    - Probe alignments with complementary sequences
    - Probe labels

    Format matches MATLAB output:
    - No '>' prefix added to lines - the '>' only appears where it exists in the sequence
    - Mask lines show the original sequence with masked positions replaced by mask characters

    Args:
        result: ProbeDesignResult from design_probes()
        filepath: Output file path
        mask_seqs: List of mask strings (same length as sequence, with mask chars at masked positions)
        line_width: Characters per line for wrapping
    """
    seq = result.input_sequence
    oligo_len = len(result.probes[0].sequence) if result.probes else 20

    # Build probe alignment string
    probe_align = [' '] * len(seq)
    probe_labels = [' '] * len(seq)

    for probe in result.probes:
        pos = probe.position
        # Get the actual sequence at this position (skip '>' characters)
        actual_seq = ''
        seq_pos = pos
        chars_collected = 0
        while chars_collected < oligo_len and seq_pos < len(seq):
            if seq[seq_pos] != '>':
                actual_seq += seq[seq_pos]
                chars_collected += 1
            seq_pos += 1

        comp_seq = complement(actual_seq)

        # Place complementary sequence (accounting for '>' in original)
        comp_idx = 0
        for i in range(oligo_len + 10):  # Allow for some '>' chars
            if pos + i >= len(probe_align):
                break
            if seq[pos + i] == '>':
                continue  # Skip '>' positions
            if comp_idx < len(comp_seq):
                probe_align[pos + i] = comp_seq[comp_idx]
                comp_idx += 1

        # Place probe label
        label = f"Prb# {probe.index},FE {probe.gibbs_fe},GC {probe.gc_percent}"
        for i, c in enumerate(label):
            if pos + i < len(probe_labels):
                probe_labels[pos + i] = c

    probe_align_str = ''.join(probe_align)
    probe_labels_str = ''.join(probe_labels)

    # Build output lines
    with open(filepath, 'w') as f:
        for start in range(0, len(seq), line_width):
            end = min(start + line_width, len(seq))

            # Original sequence (no prefix - '>' is already in sequence if present)
            f.write(f"{seq[start:end]}\n")

            # Mask sequences if provided
            if mask_seqs:
                for mask in mask_seqs:
                    if mask:
                        f.write(f"{mask[start:end]}\n")

            # Probe alignment
            f.write(f"{probe_align_str[start:end]}\n")

            # Probe labels
            f.write(f"{probe_labels_str[start:end]}\n")

            f.write("\n")


def format_probes_table(result: ProbeDesignResult) -> str:
    """Format probes as a printable table.

    Args:
        result: ProbeDesignResult from design_probes()

    Returns:
        Formatted string table
    """
    if not result.probes:
        return "No probes found."

    lines = [
        "Index\tGC%\tTm\tGibbs\tSequence\tName",
        "-" * 70,
    ]

    for probe in result.probes:
        lines.append(
            f"{probe.index}\t"
            f"{probe.gc_percent}\t"
            f"{probe.tm}\t"
            f"{probe.gibbs_fe}\t"
            f"{probe.sequence}\t"
            f"{probe.name}"
        )

    return "\n".join(lines)


def write_output_files(
    result: ProbeDesignResult,
    output_prefix: str,
    mask_seqs: Optional[List[str]] = None,
) -> None:
    """Write both oligos and seq output files.

    Args:
        result: ProbeDesignResult from design_probes()
        output_prefix: Prefix for output files (will add _oligos.txt and _seq.txt)
        mask_seqs: Optional mask sequences for seq file
    """
    write_oligos_file(result, f"{output_prefix}_oligos.txt")
    write_seq_file(result, f"{output_prefix}_seq.txt", mask_seqs)
