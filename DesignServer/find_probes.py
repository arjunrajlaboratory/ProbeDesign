"""Performs masking and design of RNA FISH probe oligos for an input FASTA.

Authors:
    Arjun Raj 2011
    Marshall J. Levesque 2011

"""
 
import sys
import bowtie_search
import fasta 
import seq
import probe_design
import re
import repeat_masker

def percent_masked(masked):
    """Print percentage of nucleotides masked in 'masked' string."""
    percent = float(masked.count('X')) / float((len(masked)-masked.count('>')))
    if percent == 0:
        print "no sequences masked"
    else:
        print "%.2f%% masked" % (100*percent)

def mask_hits(inseq,mer_length,hits,threshold):
    """Set bases of inseq to 'X' when their hit count is above threshold."""
    unmasked = inseq.one_line()
    masked = ''
    numseqs = unmasked.count('>')
    if len(hits) != len(unmasked):
        raise Exception("There should be a hit value for each inseq base")

        
    lst = list(unmasked)
    for i in range(len(hits)):
        if hits[i] > threshold:
            lst[i:(i+mer_length)] = 'X'*mer_length
    masked = ''.join(lst)
        
    percent_masked(masked)
    return masked  # now it is masked with X
        
def bowtie_align_and_mask(inseq,database,mer_length,threshold):
    """bowtie alignment & check for >= threshold hits, return masked sequence.
    
    Input:
        inseq - <fasta.Fasta> object created from the nucleotide sequence 
                we are designing probes agaisnt
        database   - <str> name of the bowtie index to align against
        mer_length - <int> length of subsequences we check for alignments
        threshold  - <int> number of hits that causes a subseq to be masked

    """
    print "Masking for %s...\t" % database,
    sys.stdout.flush()  # get the progress to print before a "lengthy" alignment
    hits = bowtie_search.align_for_hits(inseq,mer_length,database)
    mask = mask_hits(inseq,mer_length,hits,threshold)
    return mask
    
def hits_to_mask(inseq,database,mer_length,threshold):
    print "Masking for %s...\t" % database,
    sys.stdout.flush()
    hits = bowtie_search.align_for_hits(inseq,mer_length,database)
    unmasked = inseq.one_line()
    mask = [0]*len(unmasked)
    for i in range(len(hits)):
        if hits[i] > threshold:
            mask[i:(i+mer_length)] = [1]*mer_length
    return mask

def design(inseq,noligos,oligo_length,spacer_length,maskflag,species):
    """ Generic interface for designing probe oligos against FASTA seq.
    
    Input:
        inseq - <fasta.Fasta> object created from the nucleotide sequence 
                we are designing probes agaisnt
        noligos - <int>
        oligo_length - <int> between 10-30, most common is 20
        spacer_length - <int> between 0-10, most common is 2
        maskFlag - <bool> whether to perform masking steps w/ bowtie
        species - <str> name of species database to use in masking steps

    """
    if maskflag:
        #--- Perform the masking steps -----
        pmask = bowtie_align_and_mask(inseq,species+'Pseudo',16,1) # PSEUDO
        mmask = bowtie_align_and_mask(inseq,species+'Mito',7,35)   # MITOCHONDRIA
        gmask = bowtie_align_and_mask(inseq,species,16,10)         # GENOME
        '''
        rmask = repeat_mask(inseq,16)
        '''
        #---- Combine all masks --------
        unmasked_seq = inseq.one_line()
        maskBool = [False] * len(unmasked_seq)
        maskseqs = [pmask,mmask,gmask]
        for maskseq in maskseqs:
            for m in re.finditer('X',maskseq):
                maskBool[m.start()] = True
        
        masked_seq = ''
        for i,tf in enumerate(maskBool):
            if tf: 
                masked_seq += 'X'
            else:
                masked_seq += unmasked_seq[i]
        
        print "Total masking: ",
        percent_masked(masked_seq)
    else:
        masked_seq = inseq.one_line()

    # perform the oligo prop optimimzation to find probes in unmasked regions
    goodns = probe_design.findGoodness(masked_seq,oligo_length,spacer_length,probe_design.GCScore)
    #goodns = probe_design.findGoodness(masked_seq,oligo_length,spacer_length,probe_design.TmScore)
    output = probe_design.findOligos(noligos,masked_seq,goodns,oligo_length,spacer_length)
    
    results = {}
    results['masked_seq'] = masked_seq
    results['output'] = output
    return results

def designV3(inseq,noligos,oligo_length,spacer_length,maskflag,species):
    
    Inf = float('Inf')
    
    # perform the oligo prop optimimzation to find probes in unmasked regions
    #goodns = probe_design.findGoodness(masked_seq,oligo_length,spacer_length,probe_design.GCScore)
    goodns = probe_design.findGoodness(inseq.one_line(),oligo_length,spacer_length,probe_design.TmScore)
    
    badness_mask = [0]*len(goodns)

    if maskflag:
        
        hitsPseudo = [0]*len(inseq.one_line())
        hitsGen1   = [0]*len(inseq.one_line())
        hitsGen2   = [0]*len(inseq.one_line())
        hitsGen3   = [0]*len(inseq.one_line())
        hitsGen    = [0]*len(inseq.one_line())
        RMmask     = [0]*len(inseq.one_line())
        Grunmask   = [0]*len(inseq.one_line())
        Crunmask   = [0]*len(inseq.one_line())
        GCrunmask  = [0]*len(inseq.one_line())
        GCmask     = [0]*len(inseq.one_line())
        fullmask   = [0]*len(inseq.one_line())
        
        if maskflag >= 5:
            hitsPseudo = hits_to_mask(inseq,species+'Pseudo',16,0)
            hitsPseudo = probe_design.remove_short_runs(hitsPseudo,20,2)
            #print probe_design.convert_mask_to_seq(inseq.one_line(),hitsPseudo,'P')
        
        if maskflag >= 4:
            hitsGen1 = hits_to_mask(inseq,species,12,4000)
            hitsGen2 = hits_to_mask(inseq,species,14,500)
            hitsGen3 = hits_to_mask(inseq,species,16,20)
            #testseq = probe_design.convert_mask_to_seq(inseq.one_line(),hitsGen1,'B')
            #testseq = probe_design.convert_mask_to_seq(testseq,hitsGen2,'B')
            #print probe_design.convert_mask_to_seq(testseq,hitsGen3,'B')
            
        if maskflag >= 3:
            RMmask = repeat_masker.run_repeat_masker(inseq,species)
            #print len(inseq.one_line())
            #print len(RMmask)
            
        for i in range(len(hitsPseudo)):
            fullmask[i] = hitsPseudo[i] + hitsGen1[i] + hitsGen2[i] + hitsGen3[i] + RMmask[i]
        
        badness_mask = probe_design.mask_to_badness(fullmask,oligo_length)
            
        
        # The following masks are not really "masks" but rather mark places
        # where the oligo is bad.  In that case, this should just be added to 
        # badness
        if maskflag >= 2:        
            GCmask = probe_design.GC_badness(inseq.one_line(),oligo_length)
            #print probe_design.convert_mask_to_seq(inseq.one_line(),GCmask,'Z')
            
      
        if maskflag >= 1:
            Crunmask = probe_design.mask_oligos_with_runs(inseq.one_line(),'c',7,2,oligo_length)
            Grunmask = probe_design.mask_oligos_with_runs(inseq.one_line(),'g',7,2,oligo_length)
            #testseq = probe_design.convert_mask_to_seq(inseq.one_line(),Crunmask,'X')
            #print probe_design.convert_mask_to_seq(testseq,Grunmask,'X')
            
  
        for i in range(len(Crunmask)):
            tmp = Crunmask[i] + Grunmask[i] + GCmask[i]
            if tmp > 0:
                badness_mask[i] = Inf
                
                
    for i in range(len(goodns)):
        goodns[i] += badness_mask[i]

        
    output = probe_design.findOligos(noligos,inseq.one_line(),goodns,oligo_length,spacer_length)
    
    results = {}
    results['masked_seq'] = inseq.one_line()
    results['output'] = output
    return results
        
def designV3_1(inseq,noligos,oligo_length,spacer_length,maskflag,species,targetTm):
    
    Inf = float('Inf')
    
    # perform the oligo prop optimimzation to find probes in unmasked regions
    #goodns = probe_design.findGoodness(masked_seq,oligo_length,spacer_length,probe_design.GCScore)
    goodns = probe_design.findGoodness_Tm(inseq.one_line(),oligo_length,spacer_length,probe_design.TmScore_Tm,targetTm)
    
    badness_mask = [0]*len(goodns)

    if maskflag:
        
        hitsPseudo = [0]*len(inseq.one_line())
        hitsGen1   = [0]*len(inseq.one_line())
        hitsGen2   = [0]*len(inseq.one_line())
        hitsGen3   = [0]*len(inseq.one_line())
        hitsGen    = [0]*len(inseq.one_line())
        RMmask     = [0]*len(inseq.one_line())
        Grunmask   = [0]*len(inseq.one_line())
        Crunmask   = [0]*len(inseq.one_line())
        GCrunmask  = [0]*len(inseq.one_line())
        GCmask     = [0]*len(inseq.one_line())
        fullmask   = [0]*len(inseq.one_line())
        
        if maskflag >= 5:
            hitsPseudo = hits_to_mask(inseq,species+'Pseudo',16,0)
            hitsPseudo = probe_design.remove_short_runs(hitsPseudo,20,2)
            #print probe_design.convert_mask_to_seq(inseq.one_line(),hitsPseudo,'P')
        
        if maskflag >= 4:
            hitsGen1 = hits_to_mask(inseq,species,12,4000)
            hitsGen2 = hits_to_mask(inseq,species,14,500)
            hitsGen3 = hits_to_mask(inseq,species,16,20)
            #testseq = probe_design.convert_mask_to_seq(inseq.one_line(),hitsGen1,'B')
            #testseq = probe_design.convert_mask_to_seq(testseq,hitsGen2,'B')
            #print probe_design.convert_mask_to_seq(testseq,hitsGen3,'B')
            
        if maskflag >= 3:
            RMmask = repeat_masker.run_repeat_masker(inseq,species)
            #print len(inseq.one_line())
            #print len(RMmask)
            
        for i in range(len(hitsPseudo)):
            fullmask[i] = hitsPseudo[i] + hitsGen1[i] + hitsGen2[i] + hitsGen3[i] + RMmask[i]
        
        badness_mask = probe_design.mask_to_badness(fullmask,oligo_length)
            
        
        # The following masks are not really "masks" but rather mark places
        # where the oligo is bad.  In that case, this should just be added to 
        # badness
        if maskflag >= 2:        
            GCmask = probe_design.GC_badness(inseq.one_line(),oligo_length)
            #print probe_design.convert_mask_to_seq(inseq.one_line(),GCmask,'Z')
            
      
        if maskflag >= 1:
            Crunmask = probe_design.mask_oligos_with_runs(inseq.one_line(),'c',7,2,oligo_length)
            Grunmask = probe_design.mask_oligos_with_runs(inseq.one_line(),'g',7,2,oligo_length)
            #testseq = probe_design.convert_mask_to_seq(inseq.one_line(),Crunmask,'X')
            #print probe_design.convert_mask_to_seq(testseq,Grunmask,'X')
            
  
        for i in range(len(Crunmask)):
            tmp = Crunmask[i] + Grunmask[i] + GCmask[i]
            if tmp > 0:
                badness_mask[i] = Inf
                
                
    for i in range(len(goodns)):
        goodns[i] += badness_mask[i]

        
    output = probe_design.findOligos(noligos,inseq.one_line(),goodns,oligo_length,spacer_length)
    
    results = {}
    results['masked_seq'] = inseq.one_line()
    results['output'] = output
    return results
        
        
def designV4(inseq, noligos, oligo_length, spacer_length, maskflag, species, targetTm, targetRange):
    
    Inf = float('Inf')
    

    if maskflag:
        
        hitsPseudo = [0]*len(inseq.one_line())
        hitsGen1   = [0]*len(inseq.one_line())
        hitsGen2   = [0]*len(inseq.one_line())
        hitsGen3   = [0]*len(inseq.one_line())
        hitsGen    = [0]*len(inseq.one_line())
        RMmask     = [0]*len(inseq.one_line())
        Grunmask   = [0]*len(inseq.one_line())
        Crunmask   = [0]*len(inseq.one_line())
        GCrunmask  = [0]*len(inseq.one_line())
        GCmask     = [0]*len(inseq.one_line())
        fullmask   = [0]*len(inseq.one_line())
        
        if maskflag >= 5:
            hitsPseudo = hits_to_mask(inseq,species+'Pseudo',16,0)
            hitsPseudo = probe_design.remove_short_runs(hitsPseudo,20,2)
            #print probe_design.convert_mask_to_seq(inseq.one_line(),hitsPseudo,'P')
        
        if maskflag >= 4:
            hitsGen1 = hits_to_mask(inseq,species,12,4000)
            hitsGen2 = hits_to_mask(inseq,species,14,500)
            hitsGen3 = hits_to_mask(inseq,species,16,20)
            #testseq = probe_design.convert_mask_to_seq(inseq.one_line(),hitsGen1,'B')
            #testseq = probe_design.convert_mask_to_seq(testseq,hitsGen2,'B')
            #print probe_design.convert_mask_to_seq(testseq,hitsGen3,'B')
            
        if maskflag >= 3:
            RMmask = repeat_masker.run_repeat_masker(inseq,species)
            #print len(inseq.one_line())
            #print len(RMmask)
            
        for i in range(len(hitsPseudo)):
            fullmask[i] = hitsPseudo[i] + hitsGen1[i] + hitsGen2[i] + hitsGen3[i] + RMmask[i]

        # We only run the following if the mask flag is greater than 1, so I don't explicitly
        # test that condition, but it's implicit
        # if maskflag >= 1:
        # perform the oligo prop optimimzation to find probes in unmasked regions
        # Sets oligos to "bad" if they are outside the target Tm range.
        goodns = probe_design.findGoodness_RNA_DNA(inseq.one_line(), oligo_length, spacer_length, probe_design.TmScore_RNA_DNA, targetTm, targetRange)
        badness_mask = [0]*len(goodns)    
        badness_mask = probe_design.mask_to_badness(fullmask,oligo_length)
        
        for i in range(len(goodns)):  # Badness now flags all the bad oligos, so add to goodness to avoid those regions.
            goodns[i] += badness_mask[i]
    
    else:  # If there is no mask flag, we'll just do goodns with a large target range
        targetRange = [-100.0, 100.0]
        goodns = probe_design.findGoodness_RNA_DNA(inseq.one_line(), oligo_length, spacer_length, probe_design.TmScore_RNA_DNA, targetTm, targetRange)                
                
    output = probe_design.findOligos(noligos,inseq.one_line(),goodns,oligo_length,spacer_length)
    
    results = {}
    results['masked_seq'] = inseq.one_line()
    results['output'] = output
    return results
        
        
        
