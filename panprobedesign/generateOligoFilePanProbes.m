function oligos = generateOligoFilePanProbes(alignedFasta,oligoLength,allMatches,probename,outoligofile)

for i = 1:length(allMatches)
    seqFromTarget = ...
        alignedFasta(allMatches(i).whichOligo).Sequence(allMatches(i).position:(allMatches(i).position+oligoLength-1));
    oligoSeq = revcomp(seqFromTarget);
    oligos(i).seq = oligoSeq;
    
    oligos(i).gc = getGC(oligoSeq);
    oligos(i).Tm = getTm_RNA_DNA(seqFromTarget);
    oligos(i).Gibbs = getGibbs_RNA_DNA(seqFromTarget);
    oligos(i).probename = [probename '_' num2str(i)];
    oligos(i).whichOligo = allMatches(i).whichOligo;
    
end;


% Write out in a TSV
fid = fopen([outoligofile],'w');
for i = 1:length(oligos)
  fprintf(fid,'%d\t%d\t%2.1f\t%2.1f\t%s\t%s\n',i,round(oligos(i).gc),oligos(i).Tm,oligos(i).Gibbs,oligos(i).seq,oligos(i).probename);
end;
fclose(fid);
