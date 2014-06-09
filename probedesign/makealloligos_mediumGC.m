function makealloligos_mediumGC(infile,outfile);

s = readsequence(infile);

numseq = translateseq(s);

makeoligos_mediumGC

outoligofile = [outfile '_csv.txt'];
outseqfile   = [outfile '.txt'];

writeoligos
writesequence