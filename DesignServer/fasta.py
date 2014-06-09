import os
import re

class Fasta:
    """ Stores contents of a FASTA file with useful retreival methods 
        
    Usage: 
    >> fa = fasta.Fasta('/path/to/yourfavgene.fasta')
        or
    >> seq = 'ACTATCTACTACTTTCATACTTATACTCTATC'
    >> fa = fasta.Fasta(s,True)

    Author:
        Marshall J. Levesque 2011

    """

    raw = ''
    """ here we store the file contents. """

    def __init__(self,file,strflag=False):
        """ Creates a Fasta object from provided FASTA filename or string
    
        Input:
            file - /path/to a Fasta file or a raw Fasta sequence string
            strFlag  - must be True if filename argument is a Fasta string

        """
        if strflag:
            self.raw = file
        else:
            self.read_fasta(file)
            
        # The validation step (below) is for checking for "bad" characters
        # However, many people put in stuff with numbering, etc.  Easiest thing
        # to do is just to strip out those characters.  We'll do that in one_line
        #self.validate()
    
    def read_fasta(self,filename):
        """Open provided filename and read file contents into self.raw""" 
        curdirpath = os.path.join(os.getcwd(),filename)
        if os.path.exists(curdirpath):  # input was a filename in current dir
            fastapath = curdirpath
        elif os.path.exists(filename):  # input was a full path to a fasta
            fastapath = filename
        else:
            msg = "Could not locate input Fasta file %s" % filename
            raise Exception(msg)
            
        print "Opening file %s" % curdirpath
        f = open(curdirpath,'r')
        self.raw = f.read()
        f.close()
    
    def one_line(self):
        """Converts a raw FASTA string to a string w/o new-line characters
        
        We find that a useful format when dealing with FASTA sequences is 
        to treat it as single-line string (no newline or other whitespace 
        characters) and to split up a multi-sequence FASTA with '>'
        characters. This design allows for simple indexing of nucleotide
        position, if you remember to include '>'.
        
        TODO: A more 'sophisticated' approach may be to have structured data
        for FASTA headers and sequences. For now this works since current 
        methods don't care about this information, only sequence.
        """
        one_line_seq = ''
        headerLine = False
        for c in self.raw:
            if c == '>':
                headerLine = True
                one_line_seq += c
            elif c in ('\r','\n'):
                headerLine = False
            elif not headerLine:
                if c in ('a','A','c','C','t','T','g','G','n','N'):
                    one_line_seq += c

        one_line_seq = one_line_seq.rstrip()  # remove any trailing whitespace
        return one_line_seq.lower()

    def validate(self):        
        """Returns True if input Fasta string is valid format.
    
        We allow 'n' or 'N' as a valid character for pre-masked sequences.
        """
        anybad = re.search('[^>aAcCtTgGnN]',self.one_line())
        if anybad:
            msg = "Found invalid characters [^>aAcCtTgGnN] in the input sequence"
            raise Exception(msg)
        
    def to_substrings(self,mer_length,delimchar):
        ''' Convert one-line version of the FASTA into n-mer subsequences

        From the FASTA (multi) sequence file, creates a string made up of 
        n-mer subsequences separated by user-specified character. This is
        useful for performing alignments of subsequences to a reference 
        sequence and masking for hits (eg BLAST or bowtie).
        
        Arguments:
        mer_length -- int between 5-30
        delimchar  -- 't' for tab, 'c' for comma, 's' for space, '>' for
            a valid FASTA sequence with numbered headers corresponding to
            the character index when using the Fasta.one_line() function
        
        Returns the FASTA subsequence string 
        '''
        inseq = self.one_line()
        substrings = ''

        if mer_length not in range(5,31):
            raise Exception('Provided n-mer length must be between 5-30')

        if delimchar in 'tcsn>':
            if delimchar is 't':
               delimchar = '\t' 
            elif delimchar is 'c':
               delimchar = ',' 
            elif delimchar is 's':
               delimchar = ' ' 
            elif delimchar is 'n':
               delimchar = '\n' 
        else:
           delimchar = '>' 
            

        # iterate over each nucleotide char in the FASTA one-line str
        # from char[i] up to but not including i+n
        for i in range(0,len(inseq)-mer_length+1):  
            subseq = inseq[i:i+mer_length]
            if '>' in subseq:  # skip the header characters
                pass
            else:
                if delimchar == '>':
                    # make a multi-seq FASTA with char-index as header 
                    substrings += '>%d\n%s\n' % (i,subseq)
                else:
                    substrings += '%s%c' % (subseq,delimchar) # append
        
        return substrings


