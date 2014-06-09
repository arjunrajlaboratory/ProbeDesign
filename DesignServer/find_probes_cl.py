""" Command-line version of find_probes 

Description:
    A module to process command line input and running a probe design job 
    using an input FASTA sequence. Logic in the main() function is to 
    process the input arguments and then pass them along to the generic 
    find_probes method in the find_probes.pl module

Usage:
    $ python find_probes_cl.py HOTAIR.fasta human 48
                               <input_seq> <species> <noligos>

Author:
    Marshall J. Levesque 2011

"""

import sys
import fasta
import find_probes
import probe_design


def main():
    #--------------- get and check input arguments------------#
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        raise Exception("No input filename was provided")
    
    if len(sys.argv) > 2:
        prefix = sys.argv[2]  
    else:
        prefix = 'ProbeOligo'

    noligos = 0
    if len(sys.argv) > 3:
        noligos = int(sys.argv[3])

    if (not noligos) or (noligos not in range(5,151)):
        print "No valid # of oligos was provided. Using 48"
        noligos = 48
    else:
        print "Designing for %d oligos" % noligos

    if len(sys.argv) > 4:
        oligo_length = int(sys.argv[4])
    else:
        oligo_length = 20

    if len(sys.argv) > 5:
        spacer_length = int(sys.argv[5])
    else:
        spacer_length = 2

    maskflag = False
    if len(sys.argv) > 6:
        maskflag = bool(sys.argv[6])

    if maskflag:
        print "Performing masking steps"
    else:
        print "Not performing masking steps"

    if len(sys.argv) > 7:
        species =  str.lower(sys.argv[7])
    else:
        species = ''

    if maskflag:
        if not species:
            print "No species argument was provided. Using 'human'"
            species = 'human'
        elif species not in ('human','mouse','celegans','drosophila','cow','rat'):
            raise Exception("No valid species argument was provided")
        else:
            print "Using provided species: %s" % species

    #--------- Done with processing input arguments -----------# 


    # read the provided FASTA filename and store as a Fasta object
    inseq = fasta.Fasta(filename)

    # perform the (masking and) design
    results = find_probes.design(inseq,noligos,oligo_length,spacer_length,maskflag,species)
    alignment = probe_design.alignOutput(results['masked_seq'],
                                         results['output'][-1][1],oligo_length)

    # print our results
    probeset = results['output'][-1]
    probeset = probe_design.probeNames(probeset[2],prefix)
    print "%GC,OligoSeq,OligoLabel"
    for i in range(len(probeset)):
        print "%.1f,%s,%s" % ( probeset[i][0],  # %GC
                               probeset[i][1],  # probe oligo seq
                               probeset[i][2])  # probe oligo label
    strList = probe_design.splitOutput(alignment,80)
    for x in strList:
        print x


main()
