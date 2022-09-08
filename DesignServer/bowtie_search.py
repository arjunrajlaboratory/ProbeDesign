##------------------- bowtie_search.py ---------------------##
# Description: sets up and performs a bowtie search for short
#   DNA oligos
#
# Usage: 
#   >>> hits = bowtie_search(inseq,mer_length,database)
#
# Author:
#   Marshall J. Levesque 2011

import os  # use functions here in path manipulation for platform portability
from subprocess import Popen,PIPE

def find_bowtie():
    """ Find bowtie executable and return path
    
    There are common places the user may have bowtie installed and we can
    check for these, or the user can specify the path to bowtie. Our order
    of preference in searching for bowtie is (updated 11/30/2021):
        1. /usr/bin/bowtie
        2. $BOWTIEHOME/bowtie (user sets envrionment variable in .bashrc etc)
        3. $HOME/bowtie/bowtie (try users home dir)
        4. $HOME/Dropbox (RajLab)/probeDesign/bowtie/bowtie (Rajlab Dropbox)
        5. $HOME/Downloads/bowtie/bowtie (try users Downloads dir)
        6  search the current working directory (not sure how useful this is)
    """
    if os.sep == '/':  # Unix, Linux, MacOS X
        default = '/usr/bin/bowtie'
        homeDir =  os.environ['HOME']
        if 'BOWTIEHOME' in os.environ:
            bowtiedir = os.environ['BOWTIEHOME'] + os.sep + 'bowtie'
        else:
            bowtiedir = os.path.join(homeDir,'bowtie','bowtie')
        bowtiedir2 = os.path.join(homeDir,'Dropbox (RajLab)','probeDesign','bowtie','bowtie')
        bowtiedir3 = os.path.join(homeDir,'Downloads','bowtie','bowtie')
        curdir = os.getcwd() + os.sep + 'bowtie'  # current working dir
        if os.path.exists(default): btpath = default
        elif os.path.exists(bowtiedir): btpath = bowtiedir
        elif os.path.exists(bowtiedir2): btpath = bowtiedir2
        elif os.path.exists(bowtiedir3): btpath = bowtiedir3
        elif os.path.exists(curdir): btpath = curdir
        else:
            msg =  "Could not find the bowtie executable!\n" 
            msg += "Try setting the $BOWTIEHOME environment variable"
            raise Exception(msg)
    else:
        msg = "We currently do not support Windows because its difficult"
        raise Exception(msg)

    return btpath

def setup_input(v=2,k=1,database='human'):
    """ Setup the bowtie input argments. Defaults set for simple usage.
    
    -f : input reads are FASTA formatted
    -v <int> : number of allowed mismatches in sequence 
    -k <int> : number of alignments to report in output
    --quiet  : print nothing but alignments
    database : prefix of the bowtie EBWT files in $BOWTIEHOME/indexes dir
    - : specify that the input reads will come from STDIN

    """

    known_db = ('human','mouse','celegans','drosophila','cow','rat',
                'humanPseudo','mousePseudo','celegansPseudo','drosophilaPseudo','ratPseudo',
                'humanMito','mouseMito','humanReference','humanMito_rRNA','hg19',
                'human_rDNA',
                'humanMito','mouseMito','celegansMito','drosophilaMito')

    if database not in known_db:
        msg = "Provided 'species' database name is not known"
        raise Exception(msg)

    bowtie_path = find_bowtie()
    argin = [bowtie_path,'-f','--quiet',
                         '-v',str(v),'-k',str(k),
                         database,'-']
    return argin

def output_to_hits(bowtie_default_output,numchars):
    """Create list of integer hit counts for each Fasta sub-sequence.

    Main objective is to parse the default bowtie output to obtain hit counts
    for each read
    """
    
    hits = [0] * numchars  # list is same length as the Fasta.one_line() string
    for line in bowtie_default_output.split('\n'):
        cols = line.split('\t')
        if len(cols) == 8:  # default output has 8 columns, final blank line
            seqID = int(cols[0])
            numhits  = int(cols[6]) + 1  # listed as zero for single hit alignment
            hits[seqID] = numhits  # convert to ZERO-indexing

    return hits
    
        

def align_for_hits(fasta,mer_length,database):
    """Performs bowtie alignment of Fasta n-mer sub-sequences against database."""

    if mer_length not in range(5,31):
        raise Exception('Provided n-mer length should be between 5-30')
    subseqs = fasta.to_substrings(mer_length,'>')
    args = setup_input(v=0,k=1,database=database)
    bwtoutput = runbowtie(args,subseqs)

    # hits list should have an entry for every char in fasta.one_line()
    numchars = len(fasta.one_line())
    hits = output_to_hits(bwtoutput,numchars)
    
    return hits

def default_align(reads,database='human'):
    """Perform bowtie alignment with fasta formatted reads using defaults."""
    args = setup_input(database=database)
    bwtoutput = runbowtie(args,reads)
    return bwtoutput

def runbowtie(args,reads):
    
    # create a subprocess and execute the bowtie search. Have bowtie print
    # to standard out and capture that to parse out the information we want
    bwtjob = Popen(args,stdin=PIPE,stdout=PIPE, universal_newlines=True)
    stdoe = bwtjob.communicate(reads)  # stdoe[0] = stdout  stdoe[1] = stderr
    bwtoutput = stdoe[0]

    return bwtoutput
