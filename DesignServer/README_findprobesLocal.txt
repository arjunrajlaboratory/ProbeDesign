#==========================================================================
# Filename: /path/to/DesignServer/README_findprobesLocal.txt
#
# PURPOSE: The default parameters for findprobesHD masks probes that align to pseudogenes or to multiple genomic regions. For convenience, findprobesHD uses a central probe design server located in the Raj lab for performing alignments with bowtie. With the lab's move to the school of medicine in 2019, the probe design server became inaccessible to users outside of the network. To circumvent this issue, we've written this document to provide instructions for how to design probes using a local installation of bowtie and the new function finsprobesLocal.m. We attempted to make as few changes as possible to the probe design pipeline, even if additional changes would've improved the efficiency or functionality of the original code written nearly 10 years ago. In pilot testing, probes designed following these instructions were identical to the probes designed using findprobesHD run within the Raj lab. 

# Authors:
#   Benjamin Emert 2021
#==========================================================================

#--------------
# Requirements:
#--------------
- Python 2.6+ (not extensively tested for all versions, simply devel env)    
- bowtie 0.12.7 (this is the version run on the probe design Server. Alignment has not been tested using other versions of bowtie or bowtie2). 
 -- Legacy binaries can be accessed at https://sourceforge.net/projects/bowtie-bio/files/bowtie/
 -- Executable must be located at /usr/bin/bowtie, $HOME/bowtie, $HOME/Dropbox (RajLab)/probeDesign/bowtie, or $HOME/Downloads/bowtie. 
 -- indexed databases must be in the /indexes subdirectory.

#--------------
# Bowtie setup:
#--------------
Download the pre-built binary from sourceforge () or from the Raj lab Dropbox at https://www.dropbox.com/sh/6n354swv2t74lp2/AABKoCIPjeLytwasW61APoqxa?dl=0 (assuming the link has not rotted). If needed, unzip the download, change the name of the folder to 'bowtie' (e.g. bowtie-0.12.7 -> bowtie) and move the folder to one of these paths: /usr/bin/bowtie, $HOME/bowtie, $HOME/Dropbox (RajLab)/probeDesign/bowtie, or $HOME/Downloads/bowtie. On a Mac, you can rename & move a folder in the terminal by typing (for example): mv -R /original/path/bowtie-0.12.7 $HOME/Downloads/bowtie

In the terminal on Mac, you can check the installation by typing: $HOME/Downloads/bowtie -h 
You should see:
"Usage: 
  bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]" 
followed by a long list of options

The bowtie binary download from Dropbox (https://www.dropbox.com/sh/6n354swv2t74lp2/AABKoCIPjeLytwasW61APoqxa?dl=0) includes several pre-indexed genomes & sequence databases. These include all the genomes & sequence databases on the Rajlab probe design server. If you downloaded bowtie from sourceforge, you may copy the indexes from Dropbox (https://www.dropbox.com/sh/81dpffus03h0ctu/AABt-TgJgIfoaWxqV--IijkMa?dl=0) to the indexes/ subdirectory within your bowtie folder. To build a new index, use the bowtie-build function (see instructions at http://bowtie-bio.sourceforge.net/manual.shtml). Note, however, that by default the probe design software only accepts the names of genomes & sequences database installed on the probe design server (see the list of names assigned to 'okdnasource' in findprobesHD.m). If you'd like to mask probes that align to a different reference sequence, try using the 'masksequences' parameter.   

#--------------------------
# Running FindprobesLocal.m:
#---------------------------
FindprobesLocal.m can be run from your MATLAB command window using the same parameters as findprobesHD.m. For example, to design 32 RNA FISH probes against the template FASTA sequence in myinfile.fa, type the following:
>> findprobesLocal('/path/to/myinfile.fa', 32)

