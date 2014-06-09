
%function writesequence(outfile);


offset = 5;

outseq1 = untranslateseq(numseq);
outseq1 = [blanks(offset) outseq1];

outseq2 = blanks(length(outseq1));

for i = 1:length(outoligos)
  currpos = outoligos{i}.pos;
  outseq2(currpos-1+offset+actualoligo) = fliplr(outoligos{i}.letterseqrc);
  num = num2str(i);
  l = length(num);
%  outseq2(currpos-1 +2+ (1:l)) = num;
  outseq2(currpos-1 +4-(l-1)+ (1:l)) = num;

end;


linelength = 75;

fid = fopen(outseqfile,'w');
done = 0;
i = 1;
l = length(outseq2);
while ~done
  idx = i:min(i+linelength-1,l);
  fprintf(fid,'%s\n%s\n\n',outseq1(idx),outseq2(idx));
  if max(idx)==l
    done = 1;
  end
  i=i+linelength;
end