"""Core probe design algorithms including thermodynamics and optimization."""

import math
import re
from .seq_utils import sequence_okay, reverse_complement, percent_gc


def calc_score(in_score, k, new_sd):
    """Calculate running average score for dynamic programming."""
    return (in_score * k + new_sd) / (k + 1)


def tm_rna_dna(sequence):
    """Calculate Tm of a sequence using RNA-DNA energetics (SantaLucia 98 parameters).
    
    This gives the free energy (dG) at 37°C for RNA-DNA hybridization.
    
    Args:
        sequence: DNA sequence string (lowercase)
        
    Returns:
        Free energy in kcal/mol
    """
    # SantaLucia 98 parameters for RNA-DNA
    del_h = {
        'aa': -7.8, 'ac': -5.9, 'ag': -9.1, 'at': -8.3,
        'ca': -9.0, 'cc': -9.3, 'cg': -16.3, 'ct': -7.0,
        'ga': -5.5, 'gc': -8.0, 'gg': -12.8, 'gt': -7.8,
        'ta': -7.8, 'tc': -8.6, 'tg': -10.4, 'tt': -11.5
    }
            
    del_s = {
        'aa': -21.9, 'ac': -12.3, 'ag': -23.5, 'at': -23.9,
        'ca': -26.1, 'cc': -23.2, 'cg': -47.1, 'ct': -19.7,
        'ga': -13.5, 'gc': -17.1, 'gg': -31.9, 'gt': -21.6,
        'ta': -23.2, 'tc': -22.9, 'tg': -28.4, 'tt': -36.4
    }
    
    if len(sequence) < 2:
        return float('inf')
    
    try:
        dh = sum([del_h[sequence[i:i+2]] for i in range(len(sequence)-1)])
        ds = sum([del_s[sequence[i:i+2]] for i in range(len(sequence)-1)])
    except KeyError:
        return float('inf')
    
    # Terminal corrections
    dh += 1.9
    ds += -3.9
    
    # Calculate free energy at 37°C (310.15 K)
    dg = dh * 1000.0 - (37.0 + 273.15) * ds
    dg = dg / 1000
        
    return dg


def tm_score_rna_dna(sequence, target_tm, target_range):
    """Calculate Tm score for a sequence with target range check.
    
    Args:
        sequence: DNA sequence string
        target_tm: Target melting temperature
        target_range: [min, max] allowable range
        
    Returns:
        Score (distance squared from target, or inf if outside range)
    """
    if sequence_okay(sequence):
        tm = tm_rna_dna(sequence.lower())
        score = (tm - target_tm) ** 2
        if tm < target_range[0] or tm > target_range[1]:
            score = float('inf')
    else:
        score = float('inf')
        
    return score


def gc_score(sequence):
    """Calculate GC content score (distance from 45% GC)."""
    if sequence_okay(sequence):
        score = (percent_gc(sequence) - 0.45) ** 2
    else:
        score = float('inf')

    return score


def find_goodness_rna_dna(inseq, oligo_len, block_len, score_fun, target_tm, target_range):
    """Find goodness scores for all candidate oligos in sequence.
    
    The output is the "goodness", which should be minimized (it's really badness).
    
    Args:
        inseq: Input sequence string
        oligo_len: Length of oligos to design
        block_len: Minimum spacer length between oligos
        score_fun: Scoring function to use
        target_tm: Target melting temperature
        target_range: [min, max] allowable Tm range
        
    Returns:
        List of goodness scores for each possible oligo position
    """
    goodness = []
    for i in range(len(inseq) - oligo_len + 1):
        goodness.append(score_fun(inseq[i:i + oligo_len], target_tm, target_range))
    return goodness


def find_oligos(n_oligos, inseq, goodness, oligo_len, block_len):
    """Main function to find optimal oligos using dynamic programming.
    
    Args:
        n_oligos: Number of oligos to design
        inseq: Input sequence string
        goodness: List of goodness scores for each position
        oligo_len: Length of each oligo
        block_len: Minimum spacer length between oligos
        
    Returns:
        List of solutions, where each solution contains [score, matches, oligos]
    """
    # Define constants
    nan = float('nan')
    inf = float('inf')
    
    # Length of oligo plus the "blocked" area
    short_len = oligo_len + block_len
    good_len = len(inseq) - oligo_len + 1
    
    # Initialize the 2D arrays for dynamic programming
    bmsf_pos = [[]]
    bmsf_sco = [[]]
    
    # Initialize assuming the first oligo is the best match at position 0
    bmsf_pos[0] += [0]
    bmsf_sco[0] += [goodness[0]]
    
    for i in range(n_oligos - 1):
        bmsf_pos[0] += [nan]
        bmsf_sco[0] += [inf]

    # Dynamic programming to find optimal oligo positions
    for x in range(1, len(goodness)):
        # Copy the best solution from the previous iteration
        bmsf_pos += [list(bmsf_pos[x-1])]
        bmsf_sco += [list(bmsf_sco[x-1])]
        
        for k in range(n_oligos):
            potential_score = inf
            if k == 0:
                potential_score = goodness[x]
            else:
                if x >= short_len:
                    if not math.isnan(bmsf_pos[x - short_len][k-1]):
                        potential_score = calc_score(bmsf_sco[x - short_len][k-1], k, goodness[x])

            if potential_score < bmsf_sco[x][k]:
                bmsf_pos[x][k] = x
                bmsf_sco[x][k] = potential_score
        
    # Read out the oligos from the dynamic programming table
    all_matches = []
    all_scores = []
    
    for k in range(n_oligos):
        matches = []
        x = good_len - 1
        curr_k = k
        
        if not math.isnan(bmsf_pos[x][curr_k]):
            all_scores += [bmsf_sco[x][curr_k]]
            for counter in range(k + 1):
                prev_match = bmsf_pos[x][curr_k]
                matches += [int(prev_match)]
                x = int(prev_match - short_len)
                curr_k -= 1
            matches.reverse()
            all_matches += [matches]

    # Find the maximum number of valid oligos (score < 100)
    max_pos = 0
    for i in range(len(all_scores)):
        if all_scores[i] < 100:
            max_pos = i
            
    all_scores = all_scores[0:max_pos + 1]
    all_matches = all_matches[0:max_pos + 1]
    n_oligos = max_pos + 1

    # Generate the actual oligo sequences
    all_oligos = []
    for i in range(n_oligos):
        curr_matches = all_matches[i]
        all_oligos += [[reverse_complement(inseq[j:j + oligo_len]) for j in curr_matches]]
        
    # Package output
    output = []
    for i in range(n_oligos):
        output += [[all_scores[i], all_matches[i], all_oligos[i]]]

    return output


def mask_runs(inseq, the_char, run_length, mismatches):
    """Mask runs of a specific character in sequence.
    
    Args:
        inseq: Input sequence
        the_char: Character to look for runs of
        run_length: Minimum length of run to mask
        mismatches: Number of allowed mismatches in run
        
    Returns:
        Binary mask array (1 = masked, 0 = unmasked)
    """
    outseq = [0] * len(inseq)
    for i in range(len(inseq) - run_length + 1):
        count = 0
        for j in range(i, i + run_length):
            if inseq[j] == the_char:
                count += 1
        if count >= run_length - mismatches:
            outseq[i:i + run_length] = [1] * run_length
            
    return outseq


def mask_to_badness(mask, mer_length):
    """Convert a binary mask to badness scores for oligos.
    
    Args:
        mask: Binary mask array
        mer_length: Length of oligos
        
    Returns:
        Badness array where masked regions get infinite badness
    """
    inf = float('inf')
    out = list(mask)
    for i in range(len(mask)):
        if mask[i] > 0:
            for j in range(max(0, i - mer_length + 1), i + 1):
                if j < len(out):
                    out[j] = inf
    return out


def mask_oligos_with_runs(inseq, the_char, run_length, mismatches, oligo_len):
    """Mask oligos that contain runs of a specific character.
    
    Args:
        inseq: Input sequence
        the_char: Character to check for runs
        run_length: Minimum run length to trigger masking
        mismatches: Allowed mismatches in run
        oligo_len: Length of oligos
        
    Returns:
        Binary mask for oligo positions
    """
    outseq = [0] * (len(inseq) - oligo_len + 1)
    
    for i in range(len(inseq) - oligo_len + 1):
        tmp = inseq[i:i + oligo_len]
        rn = mask_runs(tmp, the_char, run_length, mismatches)
        if sum(rn) > 0:
            outseq[i] = 1
        
    return outseq


def get_acgt_content(inseq):
    """Get ACGT content fractions for a sequence."""
    length = len(inseq)
    if length == 0:
        return [0, 0, 0, 0]
        
    a = float(inseq.count('a')) / length
    c = float(inseq.count('c')) / length
    g = float(inseq.count('g')) / length
    t = float(inseq.count('t')) / length
    
    return [a, c, g, t]


def bad_acgt(inseq):
    """Check if sequence has bad ACGT composition (too much or too little C/G)."""
    acgt = get_acgt_content(inseq)
    # Flag if C or G content is >= 40% or <= 10%
    bad = acgt[1] >= 0.40 or acgt[2] >= 0.40 or acgt[1] <= 0.10 or acgt[2] <= 0.10
    return 1 if bad else 0


def gc_badness(inseq, oligo_len):
    """Calculate GC badness for all oligo positions in sequence."""
    badness = [0] * (len(inseq) - oligo_len + 1)
    
    for i in range(len(inseq) - oligo_len + 1):
        tmp = inseq[i:i + oligo_len]
        badness[i] = bad_acgt(tmp)
        
    return badness


def remove_short_runs(in_mask, n, tolerance):
    """Remove short runs from a binary mask.
    
    From an array of 1s and 0s, removes all runs of 1s <= n in length.
    For n=20, tolerance=2, removes all runs <18 in length.
    
    Args:
        in_mask: Input binary mask
        n: Minimum run length to keep
        tolerance: Allowed gaps in runs
        
    Returns:
        Filtered binary mask
    """
    out = [0] * len(in_mask)
    for i in range(len(in_mask) - n + 1):
        if in_mask[i]:
            temp = in_mask[i:i + n]
            if sum(temp) >= n - tolerance:
                for j in range(n):
                    if i + j < len(out):
                        out[i + j] = in_mask[i + j]
            
    return out