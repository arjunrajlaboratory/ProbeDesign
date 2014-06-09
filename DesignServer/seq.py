import re
import string

def stripExtraneousChars(instring,goodLetters):
    p = re.compile("[^"+goodLetters+"]")
    return(p.sub("",instring))

def reverseComplement(sequence):
    complement = string.maketrans('atcgn', 'tagcn')
    return sequence.lower().translate(complement)[::-1]
    
def complement(sequence):
    complement = string.maketrans('atcgn', 'tagcn')
    return sequence.lower().translate(complement)
    
def percentGC(sequence):
    ll = len(sequence)
    gc = len(stripExtraneousChars(sequence.lower(),'gc'))
    return float(gc)/float(ll)
