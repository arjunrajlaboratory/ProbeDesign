"""Core probe design algorithms.

This module implements the main probe design logic:
1. Calculate "badness" scores for each position based on thermodynamic properties
2. Use dynamic programming to find optimal probe placements
3. Generate probe designs from the results
"""

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from .thermodynamics import gibbs_rna_dna, tm_rna_dna
from .sequence import (
    reverse_complement,
    percent_gc,
    has_invalid_chars,
    clean_sequence,
)
from .fasta import read_fasta, sequences_to_single_string, get_basename


@dataclass
class Probe:
    """Represents a designed oligonucleotide probe."""
    index: int           # 1-based probe number
    position: int        # 0-based start position in template
    sequence: str        # Probe sequence (reverse complement of template)
    gc_percent: float    # GC content as percentage (0-100)
    tm: float            # Melting temperature (°C)
    gibbs_fe: float      # Gibbs free energy (kcal/mol)
    name: str            # Probe identifier


@dataclass
class ProbeDesignResult:
    """Result of probe design."""
    probes: List[Probe]
    score: float
    input_sequence: str
    template_name: str
    mask: Optional[List[int]] = None
    mask_strings: List[str] = field(default_factory=list)


def calculate_badness(
    seq: str,
    oligo_len: int,
    target_gibbs: float,
    allowable_range: Tuple[float, float],
) -> List[float]:
    """Calculate badness score for each position in the sequence.

    The badness is the squared difference from the target Gibbs free energy.
    Positions with invalid characters or Gibbs FE outside the allowable range
    get infinite badness.

    Args:
        seq: Input sequence (lowercase)
        oligo_len: Length of oligonucleotides
        target_gibbs: Target Gibbs free energy (kcal/mol)
        allowable_range: Tuple of (min_gibbs, max_gibbs) for valid probes

    Returns:
        List of badness scores, one per position. Length = len(seq) - oligo_len + 1
    """
    min_gibbs, max_gibbs = sorted(allowable_range)
    goodlen = len(seq) - oligo_len + 1
    badness = []

    for i in range(goodlen):
        oligo = seq[i:i + oligo_len]

        # Check for invalid/masked characters
        if has_invalid_chars(oligo):
            badness.append(float('inf'))
            continue

        # Calculate Gibbs free energy
        try:
            gibbs = gibbs_rna_dna(oligo)
        except KeyError:
            # Invalid dinucleotide (shouldn't happen with clean sequence)
            badness.append(float('inf'))
            continue

        # Check if within allowable range
        if gibbs < min_gibbs or gibbs > max_gibbs:
            badness.append(float('inf'))
            continue

        # Badness is squared distance from target
        badness.append((gibbs - target_gibbs) ** 2)

    return badness


def _calc_score(old_score: float, k: int, new_badness: float) -> float:
    """Calculate running average score.

    This implements the scoring function from the MATLAB code:
    score = (old_score * k + new_badness) / (k + 1)

    Args:
        old_score: Previous average score for k probes
        k: Number of probes so far (0-indexed, so k=0 means first probe)
        new_badness: Badness of the new probe

    Returns:
        New average score
    """
    return (old_score * k + new_badness) / (k + 1)


def find_best_probes(
    badness: List[float],
    seq_len: int,
    oligo_len: int,
    spacer_len: int,
    n_probes: int,
) -> List[Tuple[float, List[int]]]:
    """Find optimal probe positions using dynamic programming.

    This implements the MATLAB find_best_matches algorithm.

    Args:
        badness: List of badness scores per position
        seq_len: Length of input sequence
        oligo_len: Length of oligonucleotides
        spacer_len: Minimum spacing between probes
        n_probes: Maximum number of probes to find

    Returns:
        List of (score, positions) tuples for 1..n_probes solutions.
        Each positions list contains 0-based start positions.
    """
    goodlen = len(badness)
    probe_spacer_len = oligo_len + spacer_len

    # Initialize DP tables
    # bmsf_pos[x][k] = position of probe k+1 when ending at x
    # bmsf_sco[x][k] = best score for k+1 probes ending at or before x
    bmsf_pos = [[None] * n_probes for _ in range(goodlen)]
    bmsf_sco = [[float('inf')] * n_probes for _ in range(goodlen)]

    # Initialize first position
    bmsf_pos[0][0] = 0
    bmsf_sco[0][0] = badness[0]

    # Fill DP tables
    for x in range(1, goodlen):
        # Copy previous best solutions
        for k in range(n_probes):
            bmsf_pos[x][k] = bmsf_pos[x - 1][k]
            bmsf_sco[x][k] = bmsf_sco[x - 1][k]

        # Try placing probe at position x
        for k in range(n_probes):
            potential_score = float('inf')

            if k == 0:
                # First probe: just use this position's badness
                potential_score = badness[x]
            else:
                # Later probes: need to check if we can place after previous
                prev_x = x - probe_spacer_len
                if prev_x >= 0 and bmsf_pos[prev_x][k - 1] is not None:
                    potential_score = _calc_score(
                        bmsf_sco[prev_x][k - 1], k, badness[x]
                    )

            if potential_score < bmsf_sco[x][k]:
                bmsf_pos[x][k] = x
                bmsf_sco[x][k] = potential_score

    # Backtrack to find probe positions
    results = []
    for k in range(n_probes):
        x = goodlen - 1
        curr_k = k

        if bmsf_pos[x][curr_k] is None:
            continue

        score = bmsf_sco[x][curr_k]

        # Only include solutions with reasonable scores
        if score >= 1_000_000:
            continue

        positions = []
        for _ in range(k + 1):
            pos = bmsf_pos[x][curr_k]
            positions.append(pos)
            x = pos - probe_spacer_len
            curr_k -= 1

        positions.reverse()
        results.append((score, positions))

    return results


def design_probes(
    input_file: str,
    n_probes: int = 48,
    oligo_length: int = 20,
    spacer_length: int = 2,
    target_gibbs: float = -23.0,
    allowable_gibbs: Tuple[float, float] = (-26.0, -20.0),
    output_name: Optional[str] = None,
    species: str = "human",
    pseudogene_mask: bool = False,
    genome_mask: bool = False,
    index_dir: Optional[str] = None,
    repeatmask_file: Optional[str] = None,
) -> ProbeDesignResult:
    """Design oligonucleotide probes for a target sequence.

    Main entry point for probe design. Reads a FASTA file and designs
    optimal probes based on thermodynamic properties.

    Args:
        input_file: Path to input FASTA file
        n_probes: Number of probes to design (default 48)
        oligo_length: Length of each oligonucleotide (default 20)
        spacer_length: Minimum gap between probes (default 2)
        target_gibbs: Target Gibbs free energy in kcal/mol (default -23)
        allowable_gibbs: (min, max) Gibbs FE range (default (-26, -20))
        output_name: Base name for output files (default: derived from input)
        species: Species for masking databases (default: "human")
        pseudogene_mask: Whether to mask pseudogene alignments (default: False)
        genome_mask: Whether to mask repetitive genomic regions (default: False)
        index_dir: Directory containing bowtie indexes (default: auto-detect)

    Returns:
        ProbeDesignResult containing the designed probes
    """
    # Read input sequence
    headers, seqs = read_fasta(input_file)

    # Concatenate multi-entry FASTA into single sequence (with '>' marking junctions)
    full_seq = sequences_to_single_string(seqs, mark_junctions=True)

    # Check for N's in the sequence (indicates pre-repeatmasked input)
    has_n_masking = 'n' in full_seq.lower()

    # Clean sequence (lowercase, keep only actgn>)
    seq = clean_sequence(full_seq)

    # Determine output name
    if output_name is None:
        output_name = get_basename(input_file)

    # Calculate badness scores (thermodynamic filtering)
    badness = calculate_badness(
        seq, oligo_length, target_gibbs, allowable_gibbs
    )

    # Create F mask from initial badness==inf (thermodynamic filtering)
    # This is done BEFORE adding masking, matching MATLAB behavior
    # F shows positions where you CANNOT start a probe (badness==inf or past end)
    # This matches MATLAB mask_string(inseq, badness==inf, 'F')
    goodlen = len(badness)

    # Apply masking if requested
    full_mask = [0] * len(seq)
    mask_strings = []

    # Handle repeat masking from separate file or from N's in input
    if repeatmask_file:
        # Read repeatmasked sequence from separate file
        _, rm_seqs = read_fasta(repeatmask_file)
        rm_full_seq = sequences_to_single_string(rm_seqs, mark_junctions=True)
        rm_seq = clean_sequence(rm_full_seq)

        # Create repeat mask from N positions in the repeatmask file
        rmask = [1 if rm_seq[i].lower() == 'n' else 0 for i in range(len(rm_seq))]
        rstr = "".join("R" if rmask[i] else seq[i] for i in range(len(seq)))
        mask_strings.append(rstr)
        for i, v in enumerate(rmask):
            full_mask[i] += v
        print(f"Repeat masking (from file): {sum(rmask)} positions masked")
    elif has_n_masking:
        # Create repeat mask from N positions in input file
        rmask = [1 if seq[i].lower() == 'n' else 0 for i in range(len(seq))]
        rstr = "".join("R" if rmask[i] else seq[i] for i in range(len(seq)))
        mask_strings.append(rstr)
        for i, v in enumerate(rmask):
            full_mask[i] += v
        print(f"Repeat masking (from N's): {sum(rmask)} positions masked")

    if pseudogene_mask or genome_mask:
        from .masking import (
            pseudogene_mask as get_pseudogene_mask,
            genome_mask as get_genome_mask,
            mask_to_badness,
        )

        idx_dir = Path(index_dir) if index_dir else None

        if pseudogene_mask:
            try:
                pmask = get_pseudogene_mask(seq, species, idx_dir)
                for i, v in enumerate(pmask):
                    full_mask[i] += v
                # Create visualization string
                pstr = "".join("P" if pmask[i] else seq[i] for i in range(len(seq)))
                mask_strings.append(pstr)
                print(f"Pseudogene masking: {sum(pmask)} positions masked")
            except Exception as e:
                print(f"Warning: Pseudogene masking failed: {e}")

        if genome_mask:
            try:
                gmask = get_genome_mask(seq, species, idx_dir)
                for i, v in enumerate(gmask):
                    full_mask[i] += v
                # Create visualization string
                gstr = "".join("B" if gmask[i] else seq[i] for i in range(len(seq)))
                mask_strings.append(gstr)
                print(f"Genome masking: {sum(gmask)} positions masked")
            except Exception as e:
                print(f"Warning: Genome masking failed: {e}")

    # Add F mask string (thermodynamic filtering - from badness==inf BEFORE masking)
    # This matches MATLAB line 300: maskseqs = [maskseqs mask_string(inseq,badness==inf,'F')];
    # Show 'F' at position i if:
    #   - badness[i] == inf (can't start probe due to invalid chars or out-of-range Gibbs)
    #   - i >= goodlen (past the last valid probe start position)
    # Show sequence character at position i if badness[i] is finite (valid probe start)
    fstr_parts = []
    for i in range(len(seq)):
        if i < goodlen and badness[i] != float('inf'):
            fstr_parts.append(seq[i])
        else:
            fstr_parts.append('F')
    fstr = "".join(fstr_parts)
    mask_strings.append(fstr)

    # Add mask to badness (after F mask is created)
    if any(full_mask):
        from .masking import mask_to_badness
        mask_badness = mask_to_badness(full_mask, oligo_length)
        for i in range(len(badness)):
            if mask_badness[i] == float('inf'):
                badness[i] = float('inf')

    # Find optimal probe positions
    solutions = find_best_probes(
        badness, len(seq), oligo_length, spacer_length, n_probes
    )

    if not solutions:
        return ProbeDesignResult(
            probes=[],
            score=float('inf'),
            input_sequence=seq,
            template_name=output_name,
            mask=full_mask,
            mask_strings=mask_strings,
        )

    # Use the solution with the most probes (last in list)
    best_score, positions = solutions[-1]

    # Generate probe objects
    probes = []
    for i, pos in enumerate(positions):
        # Extract template region, skipping '>' characters
        template_region = ""
        j = pos
        while len(template_region) < oligo_length and j < len(seq):
            if seq[j] != '>':
                template_region += seq[j]
            j += 1

        probe_seq = reverse_complement(template_region)

        probe = Probe(
            index=i + 1,
            position=pos,
            sequence=probe_seq,
            gc_percent=round(percent_gc(probe_seq) * 100),
            tm=round(tm_rna_dna(template_region), 1),
            gibbs_fe=round(gibbs_rna_dna(template_region), 1),
            name=f"{output_name}_{i + 1}",
        )
        probes.append(probe)

    return ProbeDesignResult(
        probes=probes,
        score=best_score,
        input_sequence=seq,
        template_name=output_name,
        mask=full_mask,
        mask_strings=mask_strings,
    )
