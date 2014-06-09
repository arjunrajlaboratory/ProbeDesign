seqlen = length(numseq);

noligos = 36;

oligolength = 20;

actualoligo = 1:oligolength;  %This is just to make compatible with writesequence


positions = round(linspace(1,seqlen-oligolength,noligos));

for i = 1:noligos
  outoligos{i}.pos = positions(i);
  outoligos{i}.shortseq = numseq(positions(i):positions(i)+oligolength-1);
  outoligos{i}.strscore = structurescore(outoligos{i}.shortseq);
  outoligos{i}.gc = sum(outoligos{i}.shortseq>2.5)/oligolength;
  
%  outoligos{i}.rc = reversecomplement(outoligos{i}.seq);
  outoligos{i}.shortrc = reversecomplement(outoligos{i}.shortseq);
  outoligos{i}.letterseq = untranslateseq(outoligos{i}.shortseq);
  outoligos{i}.letterseqrc = untranslateseq(outoligos{i}.shortrc);
end;


%  
%  outoligos{i}.rc = reversecomplement(outoligos{i}.seq);
%  outoligos{i}.shortrc = reversecomplement(outoligos{i}.shortseq);
%  
%  outoligos{i}.letterseq = untranslateseq(outoligos{i}.shortseq);
%  outoligos{i}.letterseqrc = untranslateseq(outoligos{i}.shortrc);


%
%%numseq = ***;
%
%seqlen = length(numseq);
%
%noligos = 24;
%oligolength = 20+4;
%strstringency = 1;
%
%actualoligo = 3:22;
%
%for i = 1:(seqlen-oligolength)
%  candidates{i}.seq = numseq(i:i+oligolength-1);
%  candidates{i}.shortseq = candidates{i}.seq(actualoligo);
%  candidates{i}.strscore = structurescore(candidates{i}.seq);
%  candidates{i}.gc  = sum(candidates{i}.seq(actualoligo) > 2.5)/length(actualoligo);
%end
%
%clear gcs;
%
%for i = 1:length(candidates)
%  gcs(i) = candidates{i}.gc;
%  strscores(i) = candidates{i}.strscore;
%end
%
%clear outoligos
%
%tempgcs = gcs;
%
%for i = 1:length(gcs)
%  if strscores(i)>strstringency
%    tempgcs(i) = 0;
%  end
%end
%
%
%currpos = 1;
%
%done = 0;
%i = 1;
%
%while ~done
%
%%for i = 1:noligos
%  
%  j = currpos;
%  
%  %[y,j] = max(tempgcs);
%  outoligos{i} = candidates{j};
%  outoligos{i}.pos = j;
%  
%  seqtozero = max([j-oligolength,1]):min([j+oligolength,seqlen]);
%  tempgcs(seqtozero) = 0;
%  
%  oligopositions(i) = j;
%  oligogcs(i) = outoligos{i}.gc;
%  oligostrscores(i) = outoligos{i}.strscore;
%  
%  outoligos{i}.rc = reversecomplement(outoligos{i}.seq);
%  outoligos{i}.shortrc = reversecomplement(outoligos{i}.shortseq);
%  
%  outoligos{i}.letterseq = untranslateseq(outoligos{i}.shortseq);
%  outoligos{i}.letterseqrc = untranslateseq(outoligos{i}.shortrc);
%  
%  currpos = currpos + oligolength;
%  if currpos > seqlen - oligolength-1
%    done = 1;
%  end
%  i = i+1;
%  
%end
%
