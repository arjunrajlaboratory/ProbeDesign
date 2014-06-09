#==========================================================================
# FILENAME: /path/to/DesignServer/README.txt
#
# PURPOSE: Provide a summary of the contents of the DesignServer directory
#   and the steps necessary to setup the working environment
#
# Authors:
#   Arjun Raj 2011
#   Marshall J. Levesque 2011
#==========================================================================


#--------------
# REQUIREMENTS:
#--------------
- Python 2.6+ (not extensively tested for all versions, simply devel env)    
- bowtie 0.12.7 (http://bowtie-bio.sourceforge.net/index.shtml)
    installed as /usr/bin/bowtie symlink or $BOWTIEHOME dir environ variable
    indexed databases must be in $BOWTIEHOME/indexes (see "bowtie setup")

#--------------
# BOWTIE SETUP:
#--------------
After obtaining the pre-built binary from the bowtie sourceforge page or
building from source on your system, the easiest install is to set an 
environment variable in your 'bashrc' or 'bash_profile':
    export BOWTIEHOME=/path/to/bowtie-dir
Remember to open a new bash session after editing this file.

Databases supported in the code and listed in 'bowtie_search.py':
    * human (genome), humanPseudo, humanMito
    * mouse (genome), -Pseudo, -Mito
    * celegans (genome), -Pseudo, -Mito
    * drosophila (genome), -Pseudo, -Mito

These databases are either downloaded as pre-built indexes from bowtie's
sourceforge page or built using command:
    $ ./bowtie-build -a -r -f database.fa db_prefix
        (-a for auto memory management, -r for no paired-end indexing,
         -f for input being Fasta-formatted)

Shared Memory
--------------
$ ./bowtie --shmem ... # load the database into shared memory, then runs alignment
$ ./bowtie --mm ... # uses database in shared memory instead of usual I/O, runs alignment

$ ipcs -M   # see summary of shared memory settings on your system
 
To increase amount of allocatable shared memory, use the following commands:

Mac OS X:
$ sysctl -w kern.sysv.shmmax=1077936128 kern.sysv.shmall=263168

#--------------------------
# Default bowtie alignment
#--------------------------
In addition to using bowtie for masking sequences with hits to certain databases,
one can perform a bowtie alignment with default alignment parameters and one of 
the pre-built indexes used in the Raj Lab ProbeDesignServer. On the Python
interpreter prompt and having your input Fasta file in the current working dir:
    >>> import fasta
    >>> import bowtie_search as bt
    >>> hot = fasta.Fasta('HOTAIR.txt')  # Create a Fasta object instance
    Opening file HOTAIR.txt
    >>> bwtout = bt.default_align(hot.to_substrings(16,'>'),database='humanPseudo')
    >>> f = open('HOTAIR-pseudo.out', 'w')  # 16-mer reads screened against Pseudogene
    >>> f.write(bwtout)                     # with <str> output & then save to disk
    >>> f.close()


#------------------
# DESIGNING PROBES
#------------------
In order to design complimentary oligos against a Fasta-formatted input 
sequence, we use the 'find_probes_cl.py' module that runs as a script:

    $ python find_probes_cl.py input_mRNA.fa noligos maskFlag database

where the word 'database' would be replaced with your species of iterest 
(see "bowtie setup" for available databases) and 'noligos' is an <int> 
specified the desired number of oligos to design for the probe.



