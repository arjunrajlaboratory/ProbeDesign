"""FASTA file handling utilities for oligo probe design."""

import os
import re


class Fasta:
    """Stores contents of a FASTA file with useful retrieval methods.
    
    Usage:
        # From file
        fa = Fasta('/path/to/yourfavgene.fasta')
        
        # From string
        seq = 'ACTATCTACTACTTTCATACTTATACTCTATC'
        fa = Fasta(seq, from_string=True)
    """

    def __init__(self, input_data, from_string=False):
        """Creates a Fasta object from provided FASTA filename or string.
    
        Args:
            input_data: Path to FASTA file or raw FASTA sequence string
            from_string: True if input_data is a FASTA string, False if filename
        """
        if from_string:
            self.raw = input_data
        else:
            self._read_fasta(input_data)
    
    def _read_fasta(self, filename):
        """Open provided filename and read file contents into self.raw."""
        curdirpath = os.path.join(os.getcwd(), filename)
        if os.path.exists(curdirpath):
            fastapath = curdirpath
        elif os.path.exists(filename):
            fastapath = filename
        else:
            raise FileNotFoundError(f"Could not locate input FASTA file {filename}")
            
        print(f"Opening file {fastapath}")
        with open(fastapath, 'r') as f:
            self.raw = f.read()
    
    def one_line(self):
        """Converts a raw FASTA string to a string without new-line characters.
        
        Returns a single-line string with no newline or other whitespace 
        characters, splitting multi-sequence FASTA with '>' characters. 
        This allows for simple indexing of nucleotide position.
        """
        one_line_seq = ''
        header_line = False
        
        for c in self.raw:
            if c == '>':
                header_line = True
                one_line_seq += c
            elif c in ('\r', '\n'):
                header_line = False
            elif not header_line:
                if c.lower() in ('a', 'c', 't', 'g', 'n'):
                    one_line_seq += c

        one_line_seq = one_line_seq.rstrip()  # remove any trailing whitespace
        return one_line_seq.lower()

    def validate(self):        
        """Returns True if input Fasta string is valid format.
    
        Allows 'n' or 'N' as a valid character for pre-masked sequences.
        """
        anybad = re.search('[^>aAcCtTgGnN]', self.one_line())
        if anybad:
            raise ValueError("Found invalid characters [^>aAcCtTgGnN] in the input sequence")
        
    def to_substrings(self, mer_length, delim_char='>'):
        """Convert one-line version of the FASTA into n-mer subsequences.

        From the FASTA sequence file, creates a string made up of 
        n-mer subsequences separated by user-specified character. This is
        useful for performing alignments of subsequences to a reference 
        sequence and masking for hits (e.g., BLAST or bowtie).
        
        Args:
            mer_length: int between 5-30, length of subsequences
            delim_char: delimiter character - 't' for tab, 'c' for comma, 
                       's' for space, '>' for valid FASTA sequence with 
                       numbered headers corresponding to character index
        
        Returns:
            String of FASTA subsequences
        """
        inseq = self.one_line()
        substrings = ''

        if mer_length not in range(5, 31):
            raise ValueError('Provided n-mer length must be between 5-30')

        if delim_char in 'tcsn>':
            if delim_char == 't':
               delim_char = '\t' 
            elif delim_char == 'c':
               delim_char = ',' 
            elif delim_char == 's':
               delim_char = ' ' 
            elif delim_char == 'n':
               delim_char = '\n' 
        else:
           delim_char = '>' 
            
        # iterate over each nucleotide char in the FASTA one-line str
        # from char[i] up to but not including i+n
        for i in range(0, len(inseq) - mer_length + 1):  
            subseq = inseq[i:i + mer_length]
            if '>' in subseq:  # skip the header characters
                pass
            else:
                if delim_char == '>':
                    # make a multi-seq FASTA with char-index as header 
                    substrings += f'>{i}\n{subseq}\n'
                else:
                    substrings += f'{subseq}{delim_char}'
        
        return substrings