import string
import math
import re
import seq

def calcscore(inscore,k,newSD):
    return (inscore*k + newSD)/(k+1)
    

def sequenceOkay(strg, search=re.compile('[^aAcCgGtT]').search):
    """Return true if all characters are valid [aAcCgGtT] in sequence."""
    return not bool(search(strg))
    
def TmScore(sequence):
    # This gives the Tm score of a sequence with target Tm of 67.5
    if sequenceOkay(sequence):
        score = (Tm(sequence.lower())-67.5)**2
    else:
        score = float('Inf')
        
    return score

def TmScore_Tm(sequence, targetTm):
    # This gives the Tm score of a sequence with target melting temp of Tm
    if sequenceOkay(sequence):
        score = (Tm(sequence.lower())-targetTm)**2
    else:
        score = float('Inf')
        
    return score

def TmScore_RNA_DNA(sequence, targetTm, targetRange):
    # This function gives the Tm score of a sequence (i.e., distance to target Tm)    
    # It also will set the score to infinity if the Tm is outside of the targetRange
    if sequenceOkay(sequence):
        score = (Tm_RNA_DNA(sequence.lower())-targetTm)**2
        if (Tm_RNA_DNA(sequence.lower()) < targetRange[0]) or (Tm_RNA_DNA(sequence.lower()) > targetRange[1]):
            score = float('Inf')
    else:
        score = float('Inf')
        
    return score
    
def Tm_RNA_DNA(sequence):
    # This gives the Tm of a sequence using RNA-DNA energetics
    primerConc = 0.00005
    temp = 30.0
    salt = 0.33
    
    # SantaLucia 98 parameters
    delH = {'aa':-7.8, 'ac':-5.9, 'ag':-9.1, 'at':-8.3,
            'ca':-9.0, 'cc':-9.3, 'cg':-16.3,'ct':-7.0,
            'ga':-5.5, 'gc':-8.0, 'gg':-12.8,'gt':-7.8,
            'ta':-7.8, 'tc':-8.6, 'tg':-10.4,'tt':-11.5}
            
    delS = {'aa':-21.9, 'ac':-12.3, 'ag':-23.5, 'at':-23.9,
            'ca':-26.1, 'cc':-23.2, 'cg':-47.1, 'ct':-19.7,
            'ga':-13.5, 'gc':-17.1, 'gg':-31.9, 'gt':-21.6,
            'ta':-23.2, 'tc':-22.9, 'tg':-28.4, 'tt':-36.4}
    dH = 0
    dS = 0
    
    dH = sum([delH[sequence[i:i+2]] for i in range(len(sequence)-1)])
    dS = sum([delS[sequence[i:i+2]] for i in range(len(sequence)-1)])
    
    dH += 1.9
    dS += -3.9
    
    dG = dH*1000.0 - (37.0+273.15)*dS
    dG = dG/1000
        
    #ans = dH*1000/(dS + (1.9872 * math.log(primerConc/4))) + (16.6 * math.log10(salt)) - 273.15
    
#    print(delH)
#    return ans    
    return dG
    
def Tm(sequence):
    # This gives the Tm of a sequence
    primerConc = 0.00005
    temp = 30.0
    salt = 0.33
    
    # SantaLucia 98 parameters
    delH = {'aa':-7.9, 'ac':-8.4, 'ag':-7.8, 'at':-7.2,
            'ca':-8.5, 'cc':-8.0, 'cg':-10.6,'ct':-7.8,
            'ga':-8.2, 'gc':-9.8, 'gg':-8.0, 'gt':-8.4,
            'ta':-7.2, 'tc':-8.2, 'tg':-8.5, 'tt':-7.9}
            
    delS = {'aa':-22.2, 'ac':-22.4, 'ag':-21.0, 'at':-20.4,
            'ca':-22.7, 'cc':-19.9, 'cg':-27.2, 'ct':-21.0,
            'ga':-22.2, 'gc':-24.4, 'gg':-19.9, 'gt':-22.4,
            'ta':-21.3, 'tc':-22.2, 'tg':-22.7, 'tt':-22.2}
    dH = 0
    dS = 0
    
    dH = sum([delH[sequence[i:i+2]] for i in range(len(sequence)-1)])
    dS = sum([delS[sequence[i:i+2]] for i in range(len(sequence)-1)])
    
    if (sequence[0] == 'c') or (sequence[0] == 'g'):
        dH += 0.1
        dS += -2.8
    else:
        dH += 2.3
        dS += 4.1
        
    if (sequence[-1] == 'c') or (sequence[-1] == 'g'):
        dH += 0.1
        dS += -2.8
    else:
        dH += 2.3
        dS += 4.1
        
    ans = dH*1000/(dS + (1.9872 * math.log(primerConc/4))) + (16.6 * math.log10(salt)) - 273.15
    
#    print(delH)
    return ans
    
    
def GCScore(sequence):
    if sequenceOkay(sequence):
        score = (seq.percentGC(sequence)-.45)**2
    else:
        score = float('Inf')

    return score
    
def findGoodness(inseq, oligolen, blocklen, scoreFun):
    # This finds the goodness of all candidate oligos in inseq
    # The output is the "goodness", which should be minimized
    # (i.e., it's really the badness)
    # The scoreFun is the function you use to find the goodness
    # Could be TmScore or GCScore
    goodness = []
    for i in range(len(inseq)-(oligolen)):
        goodness.append(scoreFun(inseq[i:i+oligolen]))
    return goodness
    
def findGoodness_Tm(inseq, oligolen, blocklen, scoreFun, targetTm):
    # This finds the goodness of all candidate oligos in inseq
    # The output is the "goodness", which should be minimized
    # (i.e., it's really the badness)
    # The scoreFun is the function you use to find the goodness
    # Could be TmScore or GCScore
    # targetTm is the target Tm
    goodness = []
    for i in range(len(inseq)-(oligolen)):
        goodness.append(scoreFun(inseq[i:i+oligolen], targetTm))
    return goodness

def findGoodness_RNA_DNA(inseq, oligolen, blocklen, scoreFun, targetTm, targetRange):
    # This finds the goodness of all candidate oligos in inseq
    # The output is the "goodness", which should be minimized
    # (i.e., it's really the badness)
    # The scoreFun is the function you use to find the goodness
    # Could be TmScore or GCScore
    # targetTm is the target Tm
    goodness = []
    for i in range(len(inseq)-(oligolen)):
        goodness.append(scoreFun(inseq[i:i+oligolen], targetTm, targetRange))
    return goodness
                
    
    
def findOligos(noligos, inseq, goodness, oligolen, blocklen):
    # This is the main function to find oligos
    # First, you have to run findGoodness to get the oligo scores
    # The return value is a list of lists with 3 elements
    # Each element of this list is the solution for n=1,2,3...noligos
    # For each element of the list, the 3 elements are as follows:
    # First element are the scores, second element is the matches
    # third element are the oligos.
    
    
    # Define a few constants
    nan = float('Nan')
    inf = float('Inf')
    
    # Length of oligo plus the "blocked" area
    shortlen = oligolen + blocklen
    goodlen = len(inseq) - oligolen
    
    # Initialize the 2D arrays
    bmsfpos = [[]]
    bmsfsco = [[]]
    
    # Initialize the values assuming the first oligo is the best match at length 0
    bmsfpos[0] += [0]
    bmsfsco[0] += [goodness[0]]
    
    for i in range(noligos-1):
        bmsfpos[0] += [nan]
        bmsfsco[0] += [inf]

    
    for x in range(1,len(goodness)):
        # First, let's copy in the best solution from the previous iteration
        bmsfpos += [list(bmsfpos[x-1])]
        bmsfsco += [list(bmsfsco[x-1])]
        
        for k in range(noligos):
            potentialscore = inf
            if k==0:
                potentialscore = goodness[x]
            else:
                if (x >= shortlen):
                    if bmsfpos[x-shortlen][k-1] >= 0:
                        potentialscore = calcscore(bmsfsco[x-shortlen][k-1],k,goodness[x])

            
            if potentialscore < bmsfsco[x][k]:
                bmsfpos[x][k] = x
                bmsfsco[x][k] = potentialscore
        
    # Okay, now we have to read out the oligos
    
    allmatches = []
    allscores  = []
    for k in range(noligos):
        matches = []
        x = goodlen - 1
        currk = k
        if not(math.isnan(bmsfpos[x][currk])):
            allscores += [bmsfsco[x][currk]]
            for counter in range(k+1):
                prevmatch = bmsfpos[x][currk]
                matches += [prevmatch]
                x = prevmatch-shortlen
                currk -= 1
            matches.reverse()
            allmatches += list([matches])

    # Let's find the maximum number of oligos
    maxpos = 0
    for i in range(len(allscores)):
        if allscores[i] < 100:
            maxpos = i
            
    allscores = allscores[0:maxpos+1]
    allmatches = allmatches[0:maxpos+1]
    noligos = maxpos + 1

    alloligos = []
    for i in range(noligos):
        currmatches = allmatches[i]
        alloligos += [[seq.reverseComplement(inseq[j:(j+oligolen)]) for j in currmatches]]
        
    output = []
    for i in range(noligos):
        output += [ [allscores[i], allmatches[i], alloligos[i] ] ]

    return output
    
def alignOutput(inseq,matches,oligolen):
    """Uses matches from findOligos() to make a nice output w/ probes aligned to inseq

    Returns 3 elements in a list:
        [0] - the inseq
        [1] - the oligos
        [2] - probe # information
    """
    noligos = len(matches)
    compoligos  = [seq.complement(inseq[j:(j+oligolen)]) for j in matches]
    spaceoligos = [' '*(matches[i]-matches[i-1]-oligolen) for i in range(1,noligos)]
    spaceoligos = [' '*matches[0]] + spaceoligos
    probenum    = [ ('Probe # '+str(i+1)).ljust(oligolen) for i in range(noligos)]
    
    compseq  = ''.join([spaceoligos[i] + compoligos[i] for i in range(noligos)])
    probeseq = ''.join([spaceoligos[i] + probenum[i] for i in range(noligos)])
    
    compseq = compseq.ljust(len(inseq))
    probeseq = probeseq.ljust(len(inseq))
    
    return [inseq, compseq, probeseq]
    
def splitOutput(seqs, width):
    ls = [split_list(x,width) for x in seqs]
    tmp = []
    for i in range(len(ls[0])):
        for j in range(len(ls)):
            tmp += [ls[j][i]]
    return tmp
        
def split_list(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]
    
    
def probeNames(oligos,probeName):
    output = []
    for i in range(len(oligos)):
        output += [[seq.percentGC(oligos[i])*100, oligos[i], probeName + '_' + str(i+1)]]
    return output

    
def mask_runs(inseq,thechar,runlength,mismatches):
    outseq = [0]*len(inseq)
    for i in range(len(inseq)-runlength+1):
        count = 0
        for j in range(i,i+runlength):
            if inseq[j] == thechar:
                count += 1
        if count >= runlength-mismatches:
            outseq[i:i+runlength] = [1]*runlength;
            
    return outseq

def mask_to_badness(mask,mer_length):
    Inf = float('Inf')
    out = list(mask)
    for i in range(len(mask)):
        if mask[i] > 0:
            for j in range(max(0,i-mer_length+1),i+1):
                out[j] = Inf
    return out

def mask_oligos_with_runs(inseq,thechar,runlength,mismatches,oligolen):
    outseq = [0]*len(inseq)
    
    for i in range(len(inseq)-oligolen+1):
        tmp = inseq[i:i+oligolen]
        rn = mask_runs(tmp,thechar,runlength,mismatches)
        if sum(rn)>0:
            outseq[i] = 1
        
    return outseq
    

def getacgt(inseq):
    a = float(inseq.count('a'))/len(inseq)
    c = float(inseq.count('c'))/len(inseq)
    g = float(inseq.count('g'))/len(inseq)
    t = float(inseq.count('t'))/len(inseq)
    
    return [a,c,g,t]
    
    
def bad_acgt(inseq):
    tmp = getacgt(inseq)
    bd = tmp[1] >= 0.40 or tmp[2] >= 0.40 or tmp[1] <= 0.10 or tmp[2] <= 0.10
    if bd:
        out = 1
    else:
        out = 0
    return out

def GC_badness(inseq,oligolen):
    badness = [0]*len(inseq)
    
    for i in range(len(inseq)-oligolen+1):
        tmp = inseq[i:i+oligolen]
        badness[i] = bad_acgt(tmp)
        
    return badness


def remove_short_runs(inmask,n,tolerance):
# From an array of 1s and 0s, removes all runs of 1s <= n in length.
# For n=20, tolerance=2, you will remove all runs <18 in length,
# or all runs of 20 in length with 2 0s in the middle.

    out = [0]*len(inmask)
    for i in range(len(inmask)):
        if inmask[i]:
            temp = inmask[i:i+n]
            if sum(temp) >= n-tolerance:
                out[i:i+n] = inmask[i:i+n]
            
    return out

            
            
def convert_mask_to_seq(inseq,mask,chr):
# takes all places where mask > 0 and converts the character of inseq to chr
    output = list(inseq)
    
    for i in range(len(mask)):
        if mask[i]>0:
            output[i] = chr
    
    return ''.join(output)
