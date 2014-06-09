'''
@authors: 
    Marshall J. Levesque
    Daniel Wei
'''

import string
import array

from thermo import * 
import seq

def insert_base(inseq,insert):
    """Insert 0, 1, or 2 bases at the center of probe oligo or evenly spaced 
    e.g. position 13 w/in 25bp or position 8,16 in 24 oligo. Oligo gets longer.
    """
    
    seqA = array.array('c') #initialize a character array
    seqlen = len(inseq)
    
    if len(insert) == 1: #inserting one base at center of probe oligo
        seqA.fromstring(inseq[:(seqlen/2)]) #positions 1-center
        seqA.fromstring(insert) #insert a base at center position
        seqA.fromstring(inseq[(seqlen/2):]) #positions center-end
    elif len(insert) == 2:
        seqA.fromstring(inseq[:(seqlen/3)]) #first third
        seqA.fromstring(insert[0]) #insert a base
        seqA.fromstring(inseq[(seqlen/3):(2*seqlen/3)]) #second third
        seqA.fromstring(insert[1]) #insert a base
        seqA.fromstring(inseq[(2*seqlen/3):]) #last third
    else:  # no insert, just return the input sequence
        seqA.fromstring(inseq[:])
    
    return seqA.tostring()

    
def create_mask_probe(inseq,insert,end,toelength):
    """Construct probe and mask oligo sequences for an RNA target sequence.
    
    Arguments:
    inseq -- RNA sequence (5'->3' [actgACTG], uracil->thymidine) to design a probe
    sequence against, *INCLUDING* two flanking bases for overhang stability.
    insert -- 0, 1, or 2 bases to add insertions to the probe
    end -- integer 3 or 5 indicating where the mask should bind on the probe
    toelength -- length of the toehold sequence in # bp
    """
    # sanitizes input, make sure only characters [actgACTG]
    inseq = string.lower(inseq)
    assert containsAny(inseq,'bdefhijklmnopqrsvwxyz') == 0, 'Sequence contains invalid characters!'
    assert containsAny(insert,'bdefhijklmnopqrsvwxyz') == 0, 'Insert contains invalid characters!'
    assert containsAny(inseq,'u') == 0, 'DNA sequence contains u!'
    assert containsAny(insert,'u') == 0, 'Insert nucleotides cannot be uracil!'
    assert end == 5 or end == 3, 'Please specify which end of probe to mask!'
            
    #-----------------------------------------------
    # Calculate the RNA-DNA hybrid Gibbs free energy
    #-----------------------------------------------
    
    # overhang RNA/DNA data doesn't exist, taking the average of DNA/DNA & RNA/RNA
    oh5 = (overhang_rna(inseq[0:2],5) + overhang_dna(inseq[0:2],5)) / 2  # kcal/mol
    oh3 = (overhang_rna(inseq[-2:],3) + overhang_dna(inseq[-2:],3)) / 2  # kcal/mol

    inseq = inseq[1:-1]  # cut off the two flanking overhang bases
    [dHs,dSs] = stacks_rna_dna(inseq) 
    [dHi,dSi] = init_rna_dna() 
    binding_energy = gibbs(dHs+dHi,dSs+dSi,temp=37)  # cal/mol
    binding_energy += 1000*(oh5 + oh3)   # add dangling end stabilization
    binding_energy = salt_adjust(binding_energy/1000,len(inseq),saltconc=0.33)  # kcal/mol
    
    #------------------------------------------------------------------
    # Reverse complement to the a DNA probe sequence and insert one or
    # two bases evenly distributed in the sequence
    #------------------------------------------------------------------
    probe = seq.reverseComplement(inseq) # RNA target -> DNA probe
    probe = insert_base(probe,insert)  # DNA probe with bases inserted

    #------------------------------------------------------------------
    # Determine masked portion of DNA probe, depending on toehold size 
    #------------------------------------------------------------------
    mask = ''
    if end == 5:
        mask = probe[:-toelength]
        toeholdseq = probe[-toelength:]
    elif end == 3:
        mask = probe[toelength:]
        toeholdseq = probe[:toelength]
    mask = seq.reverseComplement(mask)

    #-----------------------------------------------
    # Calculate the DNA-DNA hybrid Gibbs free energy
    #-----------------------------------------------
    [dHs,dSs] = stacks_dna_dna(mask) 
    [dHi,dSi] = init_dna_dna(mask) 
    delGDD = gibbs(dHs+dHi,dSs+dSi,temp=37)  # cal/mol

    hangseq = '' # Add contribution of overhang in DNA probe + DNA mask hybrid
    if end == 5:
        hangseq = probe[len(mask)-1:len(mask)+1]
    elif end == 3:
        hangseq = probe[toelength-1:toelength+1]
    delGDD += 1000*overhang_dna(hangseq,end)  # cal/mol

    delGDD = salt_adjust(delGDD/1000,len(mask),saltconc=0.33)  # kcal/mol

    #------------------------------------------------------------------------------
    # calculation for delta delta G - the energy between RNA/DNA (target/probe) 
    # oligos and the DNA/DNA (mask/probe) oligos. We add the insertion contribution 
    # specified in SantaLucia 2004 Eqn 12 simply as +3.8kcal/mol for 1 nucleotide
    #------------------------------------------------------------------------------
    deldelG = binding_energy - delGDD + 3.8*len(insert)  # kcal/mol

    #------------------------------------------------------------------------------
    # Calculate toehold binding energy to  predict rate of strand displacement rxn.
    # Define:
    #                   Etoe = Etransition - Emasked 
    # Etransistion = dna/dna initiation + rna/dna initation + dna/dna stacks +       
    #       rna/dna stacks + rna/dna overhang stabilization + insertion penalties
    # Emasked = variable delGDD calculated above
    #------------------------------------------------------------------------------
    [dHiDD,dSiDD] = init_dna_dna(mask) 
    [dHiRD,dSiRD] = init_rna_dna()
    [dHsDD,dSsDD] = stacks_dna_dna(mask) 
    [dHsRD,dSsRD] = stacks_rna_dna(seq.reverseComplement(toeholdseq)) 
    dH = dHiDD + dHiRD + dHsDD + dHsRD  # cal/(mol*Kelvin)
    dS = dSiDD + dSiRD + dSsDD + dSsRD  # cal/(mol*Kelvin)
    dG = gibbs(dH,dS,temp=37)  # cal/mol
    if end == 5:  # this is the probe end that mask is placed
        dG += 1000*oh5  # cal/mol
    else:
        dG += 1000*oh3  # cal/mol
    dG = salt_adjust(dG/1000,len(mask)+len(toeholdseq),0.33)  # kcal/mol
    toehold_energy = dG - delGDD
    
    return (probe, mask, binding_energy, toehold_energy, deldelG)

def print_mask_design(inseq, probe, mask, binding_energy, toehold_energy, deldelG):
    #-------------------------------------------------------------------------------
    # Print out thermo results and sequences to screen as summary, sequences to file 
    #-------------------------------------------------------------------------------
    print 'Target RNA: ' + inseq
    print 'Probe: %s; %dbp' % (probe,len(probe))
    print 'Mask: %s; %dbp' % (mask, len(mask))
    #print 'Toehold1 %s %.2f kcal/mol' % (toeholdseq,dGsimple)
    toeholdseq = probe[:-len(mask)]
    print 'Toehold %s %.2f kcal/mol' % (toeholdseq,toehold_energy)
    print 'delta delta G = %.2f kcal/mol' % (deldelG)
    print 'delta G RNA/DNA = %.2f kcal/mol' % (binding_energy)
    print '\n'
    

if __name__ == '__main__':

    #chars   123456789012345678901234567890
    inseq = 'ggacgtctagaaaacgggaaccaagcaag'

    tlength = 4
    for c1 in 'c':
        # for c2 in 't':
        c2 = ''
        create_mask_probe(inseq, c1+c2, 3, 6)
        create_mask_probe(inseq, c1+c2, 3, 5)
        create_mask_probe(inseq, c1+c2, 3, 4)
        create_mask_probe(inseq, c1+c2, 3, 3)
        create_mask_probe(inseq, c1+c2, 3, 2)
