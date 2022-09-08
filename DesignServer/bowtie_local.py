"""Interface for findprobesLocal.m to run bowtie on local computer

Description:
    For users who cannot access the Raj lab probedesign server remotely, 
    and who would like to perform the default masking filters (pseudogenemask, genomemask)
    we created findprobesLocal.m to enable these maskings steps on your local computer. 
    Please see README_findprobesLocal.txt for additional instructions 
Authors:
    Benjamin L. Emert 2021

"""

import fasta
import bowtie_search
import os

from time import localtime, strftime
from argparse import ArgumentParser

# parser = ArgumentParser()
# parser.add_argument("-s", "--inseq", help = "", type = str) 
# parser.add_argument("-n", "--mer_length", help = "", type = int)
# parser.add_argument("-d", "--seq_database", help = "", type = str)

# args = parser.parse_args()

def screen_seqence(inseq,mer_length, seq_database):
	time_str = strftime("%a, %d %b %Y %H:%M:%S", localtime())
	print("Screen sequence request: {}".format(time_str))

	inseq_fasta = fasta.Fasta(inseq,strflag=True)
	hts = bowtie_search.align_for_hits(inseq_fasta,mer_length,seq_database)
	
	print("Successful screen")
	return hts

def test_function(test_str, test_int):
	time_str = strftime("%a, %d %b %Y %H:%M:%S", localtime())
	print("Screen sequence request: {}".format(time_str))

	testSeq = fasta.Fasta(test_str,strflag=True)
	testSub = testSeq.to_substrings(test_int,'>').split('\n')
	return testSub