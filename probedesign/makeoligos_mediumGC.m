
%numseq = ***;

seqlen = length(numseq);

noligos = 48;
oligolength = 20+4;
strstringency = 1;

actualoligo = 3:22;

for i = 1:(seqlen-oligolength)
  candidates{i}.seq = numseq(i:i+oligolength-1);
  candidates{i}.shortseq = candidates{i}.seq(actualoligo);
  candidates{i}.strscore = structurescore(candidates{i}.seq);
  candidates{i}.gc  = sum(candidates{i}.seq(actualoligo) > 2.5)/length(actualoligo);
end

clear gcs;

for i = 1:length(candidates)
  gcs(i) = candidates{i}.gc;
  strscores(i) = candidates{i}.strscore;
  gcend(i) = candidates{i}.shortseq(1);
  gcstub(i,:) = candidates{i}.seq(1:2)';
end

clear outoligos

tempgcs = gcs;

for i = 1:length(gcs)
  if strscores(i)>strstringency
%    tempgcs(i) = 0;
  end
  if gcend(i) > 2.5 %If there is a G or a C at the terminus
    tempgcs(i) = 0;
  end;
  if sum(gcstub(i,:)>2.5)
    tempgcs(i) = 0;
  end;
end


for i = 1:noligos
  %[y,j] = max(tempgcs);
  [y,j] = min(abs(tempgcs-.45));
  outoligos{i} = candidates{j};
  outoligos{i}.pos = j;
  
  seqtozero = max([j-oligolength,1]):min([j+oligolength,seqlen]);
  tempgcs(seqtozero) = 0;
  
  oligopositions(i) = j;
  oligogcs(i) = outoligos{i}.gc;
  oligostrscores(i) = outoligos{i}.strscore;
  
  outoligos{i}.rc = reversecomplement(outoligos{i}.seq);
  outoligos{i}.shortrc = reversecomplement(outoligos{i}.shortseq);
  
  outoligos{i}.letterseq = untranslateseq(outoligos{i}.shortseq);
  outoligos{i}.letterseqrc = untranslateseq(outoligos{i}.shortrc);
  
  
end

