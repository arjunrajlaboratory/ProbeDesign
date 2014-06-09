function outstructure = findprobes2seq(infile1,infile2,outfile,noligos)


%s = fileread('skn-1_2c.1_ORF.txt');
s = fileread(infile1);

inseq = s;

s = fileread(infile2);
inseq2 = s;

%noligos = 100;

inseq = lower(inseq);
inseq = inseq(ismember(inseq,'actg'));

inseq2 = lower(inseq2);
inseq2 = inseq2(ismember(inseq2,'actg'));

mismatches = inseq~=inseq2;

shortlen = 22;
blocklen = 2;
goodlen = length(inseq)-shortlen;
for i = 1:goodlen
  goodness(i) = sum(ismember(inseq(i:(i+shortlen-blocklen-1)),'gc'))/(shortlen-blocklen);
end;
percentgc = goodness;

mismatches2 = mismatches;
for i = 1:length(mismatches)
  mismatches2(i) = sum(mismatches(i:(min([i+shortlen-blocklen-1 length(mismatches)]))));
  
end;

mismatches2 = 1000*mismatches2(1:end-shortlen);

goodness = (goodness-.45).^2 + mismatches2;


clear bmsf;
bmsf(1,1).position = 1;
bmsf(1,1).score = goodness(1);
for i = 2:noligos
  bmsf(1,i).position = nan;
  bmsf(1,i).score = inf;
end;


for x = 2:length(goodness)
  bmsf = [bmsf ; bmsf(x-1,:)];
  
  for k = 1:noligos
    potentialscore = inf;
    if k == 1
      potentialscore = goodness(x);
    else
      if x-shortlen > 0
        if ~isnan(bmsf(x-shortlen,k-1).position)
          potentialscore = calcscore(bmsf(x-shortlen,k-1).score,k-1,goodness(x));
        end;
      end;
    end;
    if potentialscore < bmsf(x,k).score
      bmsf(x,k).position = x;
      bmsf(x,k).score = potentialscore;
    end;
  end;
end;

for k = 1:noligos
  matches = [];
  x = goodlen;
  currk = k;
  if ~isnan(bmsf(x,currk).position)
    pos{k}.score = bmsf(x,currk).score;
    for counter = 1:k
      prevmatch = bmsf(x,currk).position;
      matches = [matches prevmatch];
      x = prevmatch-shortlen;
      currk = currk - 1;
    end;
    pos{k}.matches = matches;
  end;
end;
  

outstructure = [pos{:}];
for i = 1:length(outstructure)
  outstructure(i).gcs = percentgc(outstructure(i).matches);
end

noligos = length(pos);
fprintf('Found a total of %d oligos...\n',noligos);

% Now let's output the files...

outoligofile = [outfile '_oligos.txt'];
outseqfile = [outfile '_seq.txt'];

%fid = fopen(outoligofile,'w');

oligopositions = fliplr(outstructure(noligos).matches);

for i = 1:noligos
  curroligo{i}.seq = inseq(oligopositions(i):(oligopositions(i)+shortlen-blocklen-1));
  curroligo{i}.rc = revcomp(curroligo{i}.seq);
  curroligo{i}.comp = comp(curroligo{i}.seq);
  curroligo{i}.position = oligopositions(i);
  curroligo{i}.gc = percentgc(oligopositions(i));
end;

probeseq = blanks(length(inseq));
infoseq  = blanks(length(inseq));

for i = 1:noligos
  probeseq(oligopositions(i):oligopositions(i)+shortlen-blocklen-1) = curroligo{i}.comp;
  numstr = ['Probe # ' num2str(i) ', ' num2str(round(curroligo{i}.gc*100)) '% GC'];
  infoseq(oligopositions(i):oligopositions(i) + length(numstr)-1) = numstr;
end;

%probeseq
%infoseq

fid = fopen(outoligofile,'w');
for i = 1:noligos
  fprintf(fid,'%d,%d,%d,%s\n',i,round(curroligo{i}.gc*100),curroligo{i}.position,curroligo{i}.rc);
end;
fclose(fid);



linelength = 75;

fid = fopen(outseqfile,'w');
fprintf(fid,'Processed from files %s, %s, oligos in file %s\n\n',infile1,infile2,outoligofile);
%fprintf(fid,'%d\n\n',length(oligos));

m = inseq~=inseq2;
m2 = double(m)*'x';
m3 = char(m2);  % This is the location of the mismatches;

m3 = inseq;
m3(~m) = ' ';

%fprintf('%s\n',m3);

done = 0;
i = 1;
l = length(inseq);
while ~done
  idx = i:min(i+linelength-1,l);
  i;
%  fprintf('%s\n%s\n%s\n%s\n%s\n%s\n\n',m3(idx),inseq(idx),inseq2(idx),probeseq(idx),infoseq(idx));
  
  fprintf(fid,'\n%s\n%s\n%s\n%s\n%s\n\n',m3(idx),inseq(idx),inseq2(idx),probeseq(idx),infoseq(idx));
  if max(idx)==l
    done = 1;
  end
  i=i+linelength;
end

fclose(fid);

