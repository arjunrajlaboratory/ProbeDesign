#=============================================================================
# Where does this data come from???
#=============================================================================
Pseudogene database files come from pseudogene.org by navigating from their
homepage to "DATABASE->Eukaryotic/Prokaryotic Pseudogenes" and then finding
the set for your species of interest. Pseudogene.org has chosen to use some
common names vs scientific names so be sure to double check the Tax ID to
ensure it is the correct species. By clicking 'download', a TSV file is 
downloaded to the web browser's download folder. Each TSV file is named after
the species name and build number (eg. Human61.txt). 


#=============================================================================
# How to convert from the Pseudogene.org TSV format to usable FASTA
#=============================================================================
In this directory (/rajlab/sequenceanalysis/probedesign/pseudogeneDBs) we keep
the latest build of pseudogenes for each species in their TSV format, and the
usable FASTA format created by:

********************************
*****  pseudoTSVtoFASTA.m  *****
********************************

This MATLAB script reads in the pseudogene.org TSV and converts it to usable
FASTA format. Usage of this script is described in the comments of the m-file.


----- README author -----
Marshall J. Levesque 2011
-------------------------
