function align_oligo_file(inseq,maskseqs,oligolen,currmatches,probename,outseqfile)

for i = 1:length(currmatches)
    seq = seqrcomplement(inseq(currmatches(i):(currmatches(i)+oligolen-1)));
    oligos(i).seq = seq;
    oligos(i).comp = seqcomplement(seqrcomplement(seq));
    oligos(i).gc = getGC(seq);
    oligos(i).Tm = getTm(seq);
    oligos(i).probename = [probename '_' num2str(i)];
    
end;


probeseq = blanks(length(inseq));
infoseq  = blanks(length(inseq));

for i = 1:length(currmatches)
  probeseq(currmatches(i):currmatches(i)+oligolen-1) = oligos(i).comp;
  numstr = ['Prb# ' num2str(i) ',Tm ' num2str(oligos(i).Tm,'%2.3g') ',GC ' num2str(oligos(i).gc)];
  infoseq(currmatches(i):currmatches(i) + length(numstr)-1) = numstr;
end;

fid = fopen(outseqfile,'w');

linelength = 110;

done = 0;
i = 1;
l = length(inseq);
while ~done
  idx = i:min(i+linelength-1,l);
  
  fprintf(fid,'%s\n',inseq(idx));
  for j = 1:length(maskseqs)
      fprintf(fid,'%s\n',maskseqs{j}(idx));
  end;
  fprintf(fid,'%s\n',probeseq(idx));
  fprintf(fid,'%s\n',infoseq(idx));
  fprintf(fid,'\n');
  
  if max(idx)==l
    done = 1;
  end
  i=i+linelength;
end
fclose(fid);
% 
% fprintf(fid,'%s\n',inseq);
% for i = 1:length(maskseqs);
%     fprintf(fid,'%s\n',maskseqs{i});
% end;
% fprintf(fid,'%s\n',probeseq);
% fprintf(fid,'%s\n',infoseq);
% fclose(fid);
