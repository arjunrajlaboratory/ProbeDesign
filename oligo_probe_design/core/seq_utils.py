"""Sequence utility functions for oligo probe design."""

import re


def strip_extraneous_chars(instring, good_letters):
    """Remove characters not in good_letters from instring."""
    p = re.compile("[^" + good_letters + "]")
    return p.sub("", instring)


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement_map = str.maketrans('atcgn', 'tagcn')
    return sequence.lower().translate(complement_map)[::-1]
    

def complement(sequence):
    """Return the complement of a DNA sequence."""
    complement_map = str.maketrans('atcgn', 'tagcn')
    return sequence.lower().translate(complement_map)
    

def percent_gc(sequence):
    """Calculate the GC content of a sequence as a fraction (0-1)."""
    ll = len(sequence)
    gc = len(strip_extraneous_chars(sequence.lower(), 'gc'))
    return float(gc) / float(ll)


def sequence_okay(strg, search=re.compile('[^aAcCgGtT]').search):
    """Return True if all characters are valid [aAcCgGtT] in sequence."""
    return not bool(search(strg))