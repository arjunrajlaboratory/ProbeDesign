"""Bowtie alignment utilities for oligo probe design."""

import os
import subprocess
from subprocess import PIPE, Popen


def find_bowtie():
    """Find bowtie executable and return path.
    
    Search order:
        1. /usr/bin/bowtie
        2. $BOWTIEHOME/bowtie (user sets environment variable)
        3. $HOME/bowtie/bowtie (users home dir)
        4. $HOME/Dropbox (RajLab)/probeDesign/bowtie/bowtie (Rajlab Dropbox)
        5. $HOME/Downloads/bowtie/bowtie (users Downloads dir)
        6. ./indexes/bowtie/bowtie (local package indexes dir)
    """
    if os.sep == '/':  # Unix, Linux, MacOS X
        default = '/usr/bin/bowtie'
        home_dir = os.environ['HOME']
        package_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'indexes', 'bowtie', 'bowtie')
        
        if 'BOWTIEHOME' in os.environ:
            bowtiedir = os.path.join(os.environ['BOWTIEHOME'], 'bowtie')
        else:
            bowtiedir = os.path.join(home_dir, 'bowtie', 'bowtie')
            
        bowtiedir2 = os.path.join(home_dir, 'Dropbox (RajLab)', 'probeDesign', 'bowtie', 'bowtie')
        bowtiedir3 = os.path.join(home_dir, 'Downloads', 'bowtie', 'bowtie')
        
        search_paths = [default, bowtiedir, bowtiedir2, bowtiedir3, package_dir]
        
        for path in search_paths:
            if os.path.exists(path):
                return path
                
        raise FileNotFoundError(
            "Could not find the bowtie executable!\n"
            "Try setting the $BOWTIEHOME environment variable or see README.md for installation instructions"
        )
    else:
        raise OSError("Windows is not currently supported")


def setup_bowtie_args(v=0, k=1, database='human'):
    """Setup the bowtie input arguments.
    
    Args:
        v: number of allowed mismatches in sequence 
        k: number of alignments to report in output
        database: prefix of the bowtie EBWT files
    
    Returns:
        List of command arguments for bowtie
    """
    known_db = (
        'human', 'mouse', 'celegans', 'drosophila', 'cow', 'rat',
        'humanPseudo', 'mousePseudo', 'celegansPseudo', 'drosophilaPseudo', 'ratPseudo',
        'humanMito', 'mouseMito', 'humanReference', 'humanMito_rRNA', 'hg19',
        'human_rDNA', 'celegansMito', 'drosophilaMito'
    )

    if database not in known_db:
        raise ValueError(f"Provided 'species' database name '{database}' is not known")

    bowtie_path = find_bowtie()
    args = [
        bowtie_path, '-f', '--quiet',
        '-v', str(v), '-k', str(k),
        database, '-'
    ]
    return args


def output_to_hits(bowtie_output, num_chars):
    """Create list of integer hit counts for each FASTA sub-sequence.

    Parse the default bowtie output to obtain hit counts for each read.
    
    Args:
        bowtie_output: String output from bowtie alignment
        num_chars: Length of original sequence (for creating hits array)
        
    Returns:
        List of hit counts, same length as original sequence
    """
    hits = [0] * num_chars
    
    for line in bowtie_output.split('\n'):
        cols = line.split('\t')
        if len(cols) == 8:  # default output has 8 columns
            seq_id = int(cols[0])
            num_hits = int(cols[6]) + 1  # listed as zero for single hit alignment
            if seq_id < len(hits):
                hits[seq_id] = num_hits

    return hits
        

def run_bowtie(args, reads):
    """Execute bowtie alignment with given arguments and input reads."""
    bwt_job = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    stdout, stderr = bwt_job.communicate(reads)
    
    if bwt_job.returncode != 0 and stderr:
        raise RuntimeError(f"Bowtie alignment failed: {stderr}")
    
    return stdout


def align_for_hits(fasta, mer_length, database):
    """Perform bowtie alignment of FASTA n-mer sub-sequences against database.
    
    Args:
        fasta: Fasta object containing the sequence
        mer_length: Length of subsequences to align (5-30)
        database: Name of bowtie database to align against
        
    Returns:
        List of hit counts for each position in the original sequence
    """
    if mer_length not in range(5, 31):
        raise ValueError('Provided n-mer length should be between 5-30')
        
    subseqs = fasta.to_substrings(mer_length, '>')
    args = setup_bowtie_args(v=0, k=1, database=database)
    bwt_output = run_bowtie(args, subseqs)

    # hits list should have an entry for every char in fasta.one_line()
    num_chars = len(fasta.one_line())
    hits = output_to_hits(bwt_output, num_chars)
    
    return hits


def default_align(reads, database='human'):
    """Perform bowtie alignment with FASTA formatted reads using defaults."""
    args = setup_bowtie_args(database=database)
    bwt_output = run_bowtie(args, reads)
    return bwt_output