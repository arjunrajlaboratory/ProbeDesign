function generateAlignedFilePanProbe(maskedSequences,oligos,oligoLength,allMatches,probename,outseqfile)

probeSeq = blanks(length(maskedSequences(1).Sequence));
infoSeq  = blanks(length(maskedSequences(1).Sequence));

for i = 1:length(allMatches)
    probeSeq(allMatches(i).position:allMatches(i).position+oligoLength-1) = fliplr(oligos(i).seq);
    numstr = ['Prb# ' num2str(i) ',MS ' num2str(oligos(i).whichOligo) ',GC ' num2str(oligos(i).gc)];
    infoSeq(allMatches(i).position:allMatches(i).position + length(numstr)-1) = numstr;
end;

fid = fopen(outseqfile,'w');

linelength = 110;

done = 0;
i = 1;
l = length(maskedSequences(1).Sequence);
while ~done
    idx = i:min(i+linelength-1,l);
    
    for j = 1:length(maskedSequences)
        fprintf(fid,'%s\n',maskedSequences(j).Sequence(idx));
    end;
    fprintf(fid,'%s\n',probeSeq(idx));
    fprintf(fid,'%s\n',infoSeq(idx));
    fprintf(fid,'\n');
    
    if max(idx)==l
        done = 1;
    end
    i=i+linelength;
end
fclose(fid);
