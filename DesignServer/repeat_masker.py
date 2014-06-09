##------------------- bowtie_search.py ---------------------##
# Description: sets up and performs a bowtie search for short
#   DNA oligos
#
# Usage: 
#   >>> hits = bowtie_search(inseq,mer_length,database)
#
# Author:
#   Marshall J. Levesque 2011

import fasta
import shutil
import os  # use functions here in path manipulation for platform portability
from subprocess import Popen,PIPE,call

def find_RM():
    """ Find RepeatMasker executable and return path
    
    There are common places the user may have bowtie installed and we can
    check for these, or the user can specify the path to bowtie. Our order
    of preference in searching for bowtie is:
        1. /usr/bin/RepeatMasker
        2. $REPEATMASKERHOME/RepeatMasker (user sets environment variable in .bashrc etc)
        3. $HOME/repeatmasker/RepeatMasker (try users home dir)
        4  search the current working directory (not sure how useful this is)
    """
    
    if os.sep == '/':  # Unix, Linux, MacOS X
        default = '/usr/bin/RepeatMasker'
        if os.environ.has_key('REPEATMASKERHOME'):  
            repeatmaskerdir = os.environ['REPEATMASKERHOME']
        else:
            repeatmaskerdir = os.path.join(os.environ['HOME'],'repeatmasker','RepeatMasker')
        curdir = os.getcwd() + os.sep + 'RepeatMasker'  # current working dir
        if os.path.exists(default): rmpath = default
        elif os.path.exists(repeatmaskerdir): rmpath = repeatmaskerdir
        elif os.path.exists(curdir): rmpath = curdir
        else:
            msg =  "Could not find the RepeatMasker executable!\n" 
            msg += "Try setting the $REPEATMASKERHOME environment variable"
            raise Exception(msg)
    else:
        msg = "We currently do not support Windows because its difficult"
        raise Exception(msg)

    return rmpath
        
def run_repeat_masker(inseq,species):
    # This runs repeatmasker on inseq
    
    known_db = ('human','mouse','celegans','drosophila','rat','cow')

    if species not in known_db:
        msg = "Provided 'species' database name is not known"
        raise Exception(msg)
        
    if species == 'celegans':
        species = 'elegans'  # A little workaround, since repeatmasker uses "elegans" instead of "celegans"

    RM_path = find_RM()
    tmp_path = RM_path + os.sep + 'tmp'
    
    f = open(tmp_path + os.sep + 'tmp.txt','w')
    f.write(inseq.raw)
    f.close()
    
    call([RM_path+os.sep+'RepeatMasker','-s','-species',species,'-dir',tmp_path,tmp_path + os.sep + 'tmp.txt'])

    try:
        out_fasta = fasta.Fasta(tmp_path + os.sep + 'tmp.txt.masked')
    except:
        out_fasta = inseq
    
    #print inseq.one_line()[0:3]
    #print out_fasta.one_line()[0:3]
        
    for root, dirs, files in os.walk(tmp_path):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))
    
    # One issue is that if the original file is just raw sequence (not a fasta)
    # then the output will be a fasta and so the one_line() function will return a 
    # sequence that is longer by 1 because it starts with a '>'.  So we have to account
    # for that offset.
    if inseq.one_line()[0] != '>' and out_fasta.one_line()[0] == '>':
        offset = 1
    else:
        offset = 0
    
    tmp = out_fasta.one_line().lower()
    out = [0]*len(inseq.one_line())
    inf = float('Inf')
    for i in range(len(out)):
        if tmp[i+offset] == 'n':
            out[i] = inf
    
    return out
