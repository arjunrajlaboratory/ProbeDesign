"""Main probe designer interface that orchestrates the probe design process."""

import sys
from .fasta import Fasta
from .bowtie_alignment import align_for_hits
from .probe_design import (
    find_goodness_rna_dna, find_oligos, tm_score_rna_dna,
    mask_to_badness, mask_oligos_with_runs, gc_badness, remove_short_runs
)


def design_probes(input_file, n_oligos=48, 
                 oligo_length=20, spacer_length=2,
                 target_tm=-23.0, target_range=[-26.0, -20.0],
                 species='human',
                 pseudogene_mask=True, genome_mask=True, 
                 gc_run_mask=False, gc_mask=False,  # OFF by default like MATLAB
                 repeat_mask=True,  # ON by default like MATLAB
                 output_file=None, verbose=True):
    """Design oligo probes for RNA FISH.
    
    This is the main function that implements the probe design algorithm,
    matching the behavior of the MATLAB findprobesHD function.
    
    Args:
        input_file: Path to FASTA file or Fasta object
        n_oligos: Number of oligos to design (default: 48)
        oligo_length: Length of each oligo (default: 20)
        spacer_length: Minimum spacer between oligos (default: 2)
        target_tm: Target free energy in kcal/mol (default: -23.0)
        target_range: [min, max] allowable free energy range (default: [-26, -20])
        species: Species for masking databases ('human', 'mouse', etc.)
        pseudogene_mask: Enable pseudogene masking (default: True)
        genome_mask: Enable genome masking (default: True)
        gc_run_mask: Enable GC run masking (default: False)
        gc_mask: Enable GC composition masking (default: False)
        repeat_mask: Enable repeat masking (default: True, not implemented)
        output_file: Optional output file prefix
        verbose: Print progress messages (default: True)
        
    Returns:
        Dictionary with keys:
            'oligos': List of designed oligo sequences
            'positions': List of oligo positions in input sequence
            'scores': List of oligo scores
            'masked_sequence': Input sequence with masking applied
            'results': Full results from optimization algorithm
    """
    
    # Load input sequence
    if isinstance(input_file, str):
        fasta = Fasta(input_file)
    elif isinstance(input_file, Fasta):
        fasta = input_file
    else:
        raise ValueError("input_file must be a file path or Fasta object")
    
    inseq = fasta.one_line()
    
    if verbose:
        print(f"Loaded sequence of length {len(inseq)}")
    
    # Initialize masking
    full_mask = [0] * len(inseq)
    inf = float('inf')
    
    # Apply pseudogene masking
    if pseudogene_mask:
        if verbose:
            print("Applying pseudogene masking...")
        try:
            pseudo_db = species + 'Pseudo'
            hits_pseudo = _hits_to_mask(fasta, pseudo_db, 16, 0)
            hits_pseudo = remove_short_runs(hits_pseudo, 20, 2)
            
            for i in range(len(hits_pseudo)):
                if i < len(full_mask):
                    full_mask[i] += hits_pseudo[i]
                    
            if verbose:
                masked_bases = sum(1 for x in hits_pseudo if x > 0)
                print(f"  Pseudogene masking: {masked_bases} bases masked")
        except Exception as e:
            if verbose:
                print(f"  Warning: Pseudogene masking failed (bowtie not available?): {e}")
    
    # Apply genome masking (BLAST-style)
    if genome_mask:
        if verbose:
            print("Applying genome masking...")
        try:
            # Multiple seed lengths with different thresholds (like BLAST)
            hits_gen1 = _hits_to_mask(fasta, species, 12, 4000)
            hits_gen2 = _hits_to_mask(fasta, species, 14, 500) 
            hits_gen3 = _hits_to_mask(fasta, species, 16, 20)
            
            for i in range(len(inseq)):
                if i < len(hits_gen1):
                    full_mask[i] += hits_gen1[i]
                if i < len(hits_gen2):
                    full_mask[i] += hits_gen2[i]
                if i < len(hits_gen3):
                    full_mask[i] += hits_gen3[i]
                    
            if verbose:
                masked_bases = sum(1 for i in range(len(inseq)) 
                                 if (i < len(hits_gen1) and hits_gen1[i] > 0) or
                                    (i < len(hits_gen2) and hits_gen2[i] > 0) or  
                                    (i < len(hits_gen3) and hits_gen3[i] > 0))
                print(f"  Genome masking: {masked_bases} bases masked")
        except Exception as e:
            if verbose:
                print(f"  Warning: Genome masking failed (bowtie not available?): {e}")
    
    # Apply repeat masking (placeholder)
    if repeat_mask:
        if verbose:
            print("  Warning: Repeat masking not yet implemented")
    
    # Calculate thermodynamic goodness scores
    if verbose:
        print("Calculating thermodynamic scores...")
    goodness = find_goodness_rna_dna(inseq, oligo_length, spacer_length, 
                                   tm_score_rna_dna, target_tm, target_range)
    
    # Convert sequence masks to badness scores for oligos
    badness_mask = mask_to_badness(full_mask, oligo_length)
    
    # Apply GC run masking
    if gc_run_mask:
        if verbose:
            print("Applying GC run masking...")
        c_run_mask = mask_oligos_with_runs(inseq, 'c', 7, 2, oligo_length)
        g_run_mask = mask_oligos_with_runs(inseq, 'g', 7, 2, oligo_length)
        
        for i in range(len(goodness)):
            if i < len(c_run_mask) and c_run_mask[i] > 0:
                goodness[i] = inf
            if i < len(g_run_mask) and g_run_mask[i] > 0:
                goodness[i] = inf
    
    # Apply GC composition masking  
    if gc_mask:
        if verbose:
            print("Applying GC composition masking...")
        gc_mask_scores = gc_badness(inseq, oligo_length)
        
        for i in range(len(goodness)):
            if i < len(gc_mask_scores) and gc_mask_scores[i] > 0:
                goodness[i] = inf
    
    # Add masking badness to thermodynamic scores
    for i in range(len(goodness)):
        if i < len(badness_mask):
            goodness[i] += badness_mask[i]
    
    # Find optimal oligos using dynamic programming
    if verbose:
        print("Finding optimal oligo positions...")
    results = find_oligos(n_oligos, inseq, goodness, oligo_length, spacer_length)
    
    if not results:
        if verbose:
            print("No suitable oligos found!")
        return {
            'oligos': [],
            'positions': [], 
            'scores': [],
            'masked_sequence': inseq,
            'results': []
        }
    
    # Extract results for the maximum number of oligos found
    best_result = results[-1]  # Last result has the most oligos
    oligos = best_result[2]
    positions = best_result[1]
    score = best_result[0]
    
    if verbose:
        print(f"Successfully designed {len(oligos)} oligos with average score {score:.3f}")
    
    # Save output if requested
    if output_file:
        # Gather masking information for output (matching MATLAB symbols)
        masking_info = {
            'repeat_mask': [0] * len(inseq),  # R - not implemented yet
            'pseudogene_mask': _hits_to_mask(fasta, species + 'Pseudo', 16, 0) if pseudogene_mask else [0] * len(inseq),  # P
            'genome_mask': full_mask,  # B - BLAST/genome masking
            'gc_run_mask': [0] * (len(inseq) - oligo_length + 1) if not gc_run_mask else None,  # X
            'gc_composition_mask': [0] * (len(inseq) - oligo_length + 1) if not gc_mask else gc_badness(inseq, oligo_length),  # Z
            'thermodynamic_mask': goodness,  # F - free energy masking
            'oligo_length': oligo_length
        }
        
        # Add GC run masking if enabled
        if gc_run_mask:
            c_runs = mask_oligos_with_runs(inseq, 'c', 7, 2, oligo_length)
            g_runs = mask_oligos_with_runs(inseq, 'g', 7, 2, oligo_length) 
            masking_info['gc_run_mask'] = [c_runs[i] + g_runs[i] if i < len(c_runs) and i < len(g_runs) else 0 
                                         for i in range(len(inseq) - oligo_length + 1)]
        else:
            masking_info['gc_run_mask'] = [0] * (len(inseq) - oligo_length + 1)
            
        _save_output(oligos, positions, inseq, output_file, oligo_length, masking_info)
    
    return {
        'oligos': oligos,
        'positions': positions,
        'scores': [score] * len(oligos),  # Individual scores not tracked in current algorithm
        'masked_sequence': inseq,
        'results': results
    }


def _hits_to_mask(fasta, database, mer_length, threshold):
    """Convert bowtie hits to binary mask.
    
    Args:
        fasta: Fasta object with sequence
        database: Bowtie database name
        mer_length: Length of k-mers to align  
        threshold: Hit count threshold for masking
        
    Returns:
        Binary mask array (1 = masked, 0 = unmasked)
    """
    try:
        hits = align_for_hits(fasta, mer_length, database)
        mask = [0] * len(hits)
        
        for i in range(len(hits)):
            if hits[i] > threshold:
                for j in range(i, min(i + mer_length, len(mask))):
                    mask[j] = 1
                    
        return mask
    except Exception as e:
        # If alignment fails, return no masking
        print(f"Warning: Alignment against {database} failed: {e}")
        return [0] * len(fasta.one_line())


def _save_output(oligos, positions, sequence, output_file, oligo_length, masking_info=None):
    """Save probe design output to files in the classic MATLAB format.
    
    Args:
        oligos: List of oligo sequences
        positions: List of oligo positions
        sequence: Original sequence
        output_file: Output file prefix
        oligo_length: Length of oligos
        masking_info: Dictionary with masking information
    """
    from .seq_utils import percent_gc
    from .probe_design import tm_rna_dna
    
    # Save oligos in simple format
    with open(f"{output_file}_oligos.txt", 'w') as f:
        for i, oligo in enumerate(oligos):
            f.write(f"Probe_{i+1:02d}\t{oligo}\n")
    
    # Save alignment view in classic MATLAB format
    with open(f"{output_file}_alignment.txt", 'w') as f:
        
        # Process sequence in chunks of ~100 characters for readability
        chunk_size = 100
        for chunk_start in range(0, len(sequence), chunk_size):
            chunk_end = min(chunk_start + chunk_size, len(sequence))
            chunk = sequence[chunk_start:chunk_end]
            
            # Line 1: Original sequence (as-is from FASTA)
            f.write(f"{chunk}\n")
            
            # Line 2: Original sequence (repeated)
            f.write(f"{chunk}\n")
            
            # Line 3: Sequence with pseudogene masking (P markers)
            pseudo_line = list(chunk)
            if masking_info and 'pseudogene_mask' in masking_info:
                pseudogene_mask = masking_info['pseudogene_mask']
                for j in range(len(chunk)):
                    abs_pos = chunk_start + j
                    if abs_pos < len(pseudogene_mask) and pseudogene_mask[abs_pos] > 0:
                        pseudo_line[j] = 'P'
            f.write(f"{''.join(pseudo_line)}\n")
            
            # Line 4: Sequence with genome masking (B markers)  
            genome_line = list(chunk)
            if masking_info and 'genome_mask' in masking_info:
                genome_mask = masking_info['genome_mask']
                for j in range(len(chunk)):
                    abs_pos = chunk_start + j
                    if abs_pos < len(genome_mask) and genome_mask[abs_pos] > 0:
                        genome_line[j] = 'B'
            f.write(f"{''.join(genome_line)}\n")
            
            # Line 5: Thermodynamic masking (F markers for free energy)
            # This shows positions where oligos have bad thermodynamics
            mask_chunk = chunk[1:] if chunk.startswith('>') else chunk
            mask_start_offset = 1 if chunk.startswith('>') else 0
            thermo_line = list(mask_chunk)
            
            if masking_info and 'thermodynamic_mask' in masking_info:
                thermo_mask = masking_info['thermodynamic_mask']
                oligo_len = masking_info.get('oligo_length', 20)
                
                for j in range(len(mask_chunk)):
                    abs_pos = chunk_start + j + mask_start_offset
                    
                    # Check if this position would be masked by thermodynamics
                    is_thermo_masked = False
                    if abs_pos < len(thermo_mask) and thermo_mask[abs_pos] == float('inf'):
                        is_thermo_masked = True
                    
                    # Also check other mask types for combined display
                    is_other_masked = False
                    if 'pseudogene_mask' in masking_info and abs_pos < len(masking_info['pseudogene_mask']):
                        is_other_masked = is_other_masked or masking_info['pseudogene_mask'][abs_pos] > 0
                    if 'genome_mask' in masking_info and abs_pos < len(masking_info['genome_mask']):
                        is_other_masked = is_other_masked or masking_info['genome_mask'][abs_pos] > 0
                    
                    # Check oligo-level masks
                    oligo_pos = abs_pos
                    if oligo_pos < len(sequence) - oligo_len + 1:
                        if 'gc_run_mask' in masking_info and oligo_pos < len(masking_info['gc_run_mask']):
                            is_other_masked = is_other_masked or masking_info['gc_run_mask'][oligo_pos] > 0
                        if 'gc_composition_mask' in masking_info and oligo_pos < len(masking_info['gc_composition_mask']):
                            is_other_masked = is_other_masked or masking_info['gc_composition_mask'][oligo_pos] > 0
                    
                    if is_thermo_masked or is_other_masked:
                        thermo_line[j] = 'F'
                    else:
                        thermo_line[j] = mask_chunk[j].lower()
            
            # Add the > prefix to match other lines            
            if chunk.startswith('>'):
                f.write(f">{''.join(thermo_line)}\n")
            else:
                f.write(f"{''.join(thermo_line)}\n")
            
            # Line 6: Probe sequences aligned to their positions  
            # Show what would base-pair with the target (complement of target region)
            from .seq_utils import reverse_complement
            probe_chunk_len = len(thermo_line)
            probe_line = [' '] * probe_chunk_len
            for i, pos in enumerate(positions):
                # Position is relative to the sequence without >
                # But we need to position relative to the current chunk (also without >)
                if chunk_start + mask_start_offset <= pos < chunk_end:
                    rel_pos = pos - (chunk_start + mask_start_offset)
                    # Extract the target sequence region at this position
                    target_start = pos
                    target_end = pos + oligo_length
                    if target_end <= len(sequence):
                        # Get the target sequence (skip > if present)
                        if sequence.startswith('>'):
                            target_region = sequence[target_start+1:target_end+1]
                        else:
                            target_region = sequence[target_start:target_end]
                        # Show reverse complement of target (what would base-pair)
                        display_seq = reverse_complement(target_region)
                        # Place display sequence at correct position
                        for j, base in enumerate(display_seq):
                            if 0 <= rel_pos + j < len(probe_line):
                                probe_line[rel_pos + j] = base
            f.write(f"{''.join(probe_line)}\n")
            
            # Line 7: Probe information (number, FE, GC content)  
            info_line = [' '] * probe_chunk_len
            for i, pos in enumerate(positions):
                # Position is relative to the sequence without >
                if chunk_start + mask_start_offset <= pos < chunk_end:
                    rel_pos = pos - (chunk_start + mask_start_offset)
                    oligo = oligos[i]
                    
                    # Calculate free energy and GC content
                    try:
                        fe = tm_rna_dna(oligo.lower())
                        gc = percent_gc(oligo) * 100
                        info_text = f"Prb# {i+1},FE {fe:.1f},GC {gc:.0f}"
                    except:
                        info_text = f"Prb# {i+1}"
                    
                    # Place info text at probe position
                    for j, char in enumerate(info_text):
                        if 0 <= rel_pos + j < len(info_line):
                            info_line[rel_pos + j] = char
            f.write(f"{''.join(info_line)}\n")
            
            # Empty line between chunks
            f.write("\n")
    
    print(f"Output saved to {output_file}_oligos.txt and {output_file}_alignment.txt")