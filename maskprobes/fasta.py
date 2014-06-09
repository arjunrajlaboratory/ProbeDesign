import os
import re
import string

class Fasta:
    """ Stores contents of a FASTA file with useful retreival methods 
        
    Usage: 
    >> fa = fasta.Fasta('/path/to/yourfavgene.fasta')
        or
    >> seq = 'ACTATCTACTACTTTCATACTTATACTCTATC'
    >> fa = fasta.Fasta(s,True)

    Author:
        Marshall J. Levesque 2011-2013

    """

    
    
    """ here we store the file contents. """

    def __init__(self,file,strflag=False):
        """ Creates a Fasta object from provided FASTA filename or string
    
        Input:
            file - /path/to a Fasta file or a raw Fasta sequence string
            strFlag  - must be True if filename argument is a Fasta string

        """

        self.path     = ''
        self.filename = ''
        self.headers  = []
        self.seqs     = []

        self.introns  = []
        self.exons    = []
        self.CDS      = []

        if strflag:
            self.parse_fasta(file)
        else:
            self.read_fasta_file(file)

            
        # The validation step (below) is for checking for "bad" characters
        # However, many people put in stuff with numbering, etc.  Easiest thing
        # to do is just to strip out those characters.  We'll do that in one_line
        #self.validate()
    
    def read_fasta_file(self,filename):
        """Open provided filename,read and parse file contents""" 
        curdirpath = os.path.join(os.getcwd(),filename)
        if os.path.exists(curdirpath):  # input was a filename in current dir
            self.path = os.getcwd()
            self.filename = filename
        elif os.path.exists(filename):  # input was a full path to a fasta
            (self.path,self.filename) = os.path.split(filename)
        else:
            msg = "Could not locate input Fasta file %s" % filename
            raise Exception(msg)
            
        print "Opening file %s" % curdirpath
        f = open(os.path.join(self.path,self.filename),'r')
        self.parse_fasta(f.read())
        self.parse_headers()
        f.close()

    def parse_fasta(self,raw):
        """Breaks up fasta file into .headers and .seqs"""

        seq = ''  # temporary sequence storage we build line-by-line
        lines = raw.split('\n')
        for line in lines: # for each line in the FASTA
            line.strip()  # remove leading and trailing whitespace
            if line is '':  # empty line
                continue

            if line[0] is '>': # we are at a header line

                if seq is not '':  # >1 entry in this fasta, append seqs
                    self.seqs.append(seq)
                    seq = ''

                self.headers.append(line)
            else:
                seq += line

        if seq is not '':  # add the last or only sequence
            self.seqs.append(seq)

        #TODO: Validation of the fasta sequences and parsing of header info"""

        if len(self.seqs) == 0: 
            msg = "No sequences in FASTA %s" % self.filename
            raise Exception(msg)

        self.mark_exon_intron_CDS()

    def mark_exon_intron_CDS(self):
        """Identify each sequence as EXON,CDS,INTRON

        Given the FASTA uses the following UCSC genomic sequence export style:

                1) EXONS in UPPERCASE, everything else lowercase
                2) One FASTA record per region (exon, intron, etc.)
                3) Split UTR and CDS parts of an exon into separate FASTA records

        Then we can label the sequences in the Fasta.seqs list as exon, intron, CDS
        and store this info as a list of indicies

        CDS sequences are a subset of the EXON sequences. So for example, we could
        have EXONS and CDS annotation like this:
            self.EXONS = [0 1 3 5 7 8]
            self.CDS   = [1 3 5 7]
        which indicates that there is a 5'UTR exon adjacent to the translation start
        codon and a 3'UTR just after the stop codon.

        This logic and tracking of parts of a gene is separate from FASTA parsing 
        so it can change when there are updates to UCSC or the way we want to handle
        this type of data."""

        START_CODONS = ['ATG']
        STOP_CODONS  = ['TGA','ACC']
        CDS_FLAG = False
        
        for (i,seq) in enumerate(self.seqs):
            if seq.islower():  # intron
                self.introns.append(i)
            else:  # EXON
                self.exons.append(i)
                
                if self.seqs[i][0:3] in START_CODONS:
                    CDS_FLAG = True

                if CDS_FLAG:
                    self.CDS.append(i)

                if self.seqs[i][-3:] in STOP_CODONS:  # 3' UTR
                    CDS_FLAG = False


    def parse_headers(self):
        """Gets the information out of UCSC fasta headers
        example: 
            >mm9_refGene_NM_010028_0 range=chrX:12858148-12858238 
        
        We are hoping to find these vaules in the header for each sequence:
            seqID: str identifying the sequence
            chr: str for chromosome (e.g. chr19 or chrX)
            start: int for starting genomic position
            end: int for starting genomic position
            strand: str (+/-)
            repeats: str indicating if there is repeat masking ('none' or 'N' or 'n')
        
        There are object instance attributes  with these same names
        """
        n = len(self.seqs) 
        self.seqID   = ['']*n
        self.chr     = ['']*n
        self.strand  = ['']*n
        self.start   = [-1]*n
        self.end     = [-1]*n
        self.repeats = ['']*n

        re_seqID = re.compile('>(\w+)')
        re_chr_and_range = re.compile('range=(\w+):(\d+)-(\d+)')
        re_strand  = re.compile('strand=(\W)')
        re_repeats = re.compile('repeatMasking=(\w+)')
        
        for (i,head) in enumerate(self.headers):
            attrs = head.split(' ')
            for a in attrs:
                """Go through all known header elements using string matching"""
                if re_seqID.match(a):
                    # TODO: We can do genomic build here
                    self.seqID[i] = re_seqID.match(a).group(1)
                elif re_chr_and_range.match(a): # range=chrX:12858148-12858238
                    results = re_chr_and_range.match(a)
                    self.chr[i] = results.group(1)
                    self.start[i] = int(results.group(2))
                    self.end[i] = int(results.group(3))
                elif re_strand.match(a):
                    self.strand[i] = re_strand.match(a).group(1)
                elif re_repeats.match(a):
                    self.repeats[i] = re_repeats.match(a).group(1)
        

    def seq_in_range(self,start,end,onlyexon=False):
        """Returns (str) sequence for the provided range (start to end on + strand)
        
        The returned sequence is inclusive for the position coordinates provided.
        
        Input:
            start - (int) 5' genomic position on + strand 
            end   - (int) 3' genomic position on + strand 
        """
        # check if we even have position information extracted from headers
        if -1 in self.start or -1 in self.end:
            raise Exception('Missing sequence range information in this FASTA')

        # check if we even have strand information extracted from headers
        if '' in self.strand:
            raise Exception('Missing DNA strand information in this FASTA')

        # check if we were given numbers

        # sanity checks
        if end < start:
            msg = 'Provided range [' + str(start) + ' - ' + str(end) + '] '
            msg += 'was not valid for + strand genomic coordinates'
            raise Exception(msg)

        # TODO: ensure that start/end values match sequence lengths
        # TODO: ensure that start/end values are exclusive and continuous

        # Find the sequences in the FASTA that contain the start and end positions
        # of the provided range.  Ensure they make sense
        seq_start = -1
        seq_end = -1
        for i in range(len(self.seqs)):
            if start >= self.start[i] and start <= self.end[i]:
                seq_start = i

            if end >= self.start[i] and end <= self.end[i]:
                seq_end = i
        
        if seq_start == -1 or seq_end == -1:
            msg  = 'Could not find sequences for the provided range: '
            msg += '[' + str(start) + ' - ' + str(end) + '] '
            raise Exception(msg)

        # For both the start and end genomic posiiton, we have the seq_index
        # We need to map genomic position to index of the characters in the seqs. 
        if self.strand[seq_start] == '+':
            index_start = start - self.start[seq_start]
        else:  # minus strand
            index_start = self.start[seq_start] - start

        if self.strand[seq_end] == '+':
            index_end =  end - self.end[seq_end]
        else:  # minus strand
            index_end = self.end[seq_end] - end

        # construct the eequence
        subseq = ''
        seqs = [seq_start,seq_end]
        seqs.sort()
        for i in range(seqs[0],seqs[1]+1):
    
            if onlyexon and i in self.introns:
                continue
            
            from_index = to_index = None

            if i == seq_start:
                if self.strand[seq_start] == '+':
                    from_index = index_start
                else:  # minus strand
                    to_index = index_start
            
            if i == seq_end:
                if self.strand[seq_end] == '+':
                    to_index = index_end
                else:  # minus strand
                    from_index = index_end
            if to_index == 0:   # handles [n:0] case
                to_index = None

            subseq += self.seqs[i][from_index:to_index] 
 
        return subseq
        
    def range_of_seq(self,seq):
        """A substring match but returns genomic coordinates [5',3']
        """
        oneseq = ''
        if self.strand[0] == '-':
            seq = seq.
        for i in range(
    
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

    # TODO Look up python static method styling

    def revcomp(self,seq=''):
        """The reverse complement string for the sequence provided. If no
            sequence is provided, we used the FASTA's entire sequence
        Input Args:
            seq - (str) 
        """

        # if no sequence is provided, use this FASTA's entire sequence
        if seq == '':
            for i in len(self.seqs):
                seq += self.seqs[i]

        complement = string.maketrans('ATCGNatcgn', 'ATCGNtagcn')
        return seq.lower().translate(complement)[::-1]
