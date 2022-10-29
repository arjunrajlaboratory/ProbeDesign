"""Interface for findprobesLocal.m to screen fasta sequence using bowtie on local computer

Description:
    For users who cannot access the Raj lab probedesign server remotely, 
    and who would like to perform the default masking filters (pseudogenemask, genomemask)
    we created findprobesLocal.m to enable these maskings steps on your local computer. 
    Please see README_findprobesLocal.txt for additional instructions 
Authors:
    Benjamin L. Emert 2022

"""

import fasta
import bowtie_search
import os
import sys

from time import localtime, strftime
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inseq", type = str, help = "",) 
parser.add_argument("mer_length", type = int, help = "")
parser.add_argument("seq_database", type = str, help = "")

args = parser.parse_args()

def screen_seqence(inseq,mer_length, seq_database):
    # time_str = strftime("%a, %d %b %Y %H:%M:%S", localtime())
    # sys.stderr.write("Screen sequence request: {}\n".format(time_str))
    inseq_fasta = fasta.Fasta(inseq,strflag=True)

    inseq_fasta = fasta.Fasta(args.inseq,strflag=True)
    hts = bowtie_search.align_for_hits(inseq_fasta,mer_length,seq_database)
    
    # sys.stderr.write("Successful screen\n")
    return hts

print(screen_seqence(args.inseq, args.mer_length, args.seq_database))