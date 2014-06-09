function [outstructure,oligos] = findprobes_introns_mjl(infile,outfile,dnasource,noligos,pseudogenemaskflag)
%
% pseudogenemaskflag = 1     mask pseudogenes
% pseudogenemaskflag = 0    do not mask pseudogenes
% 
% FINDPROBES_INTRONS_MJL  Create specific oligo probes against target mRNA sequence
%
%   Input Args:
%   infile    - the template FASTA sequence you are designing probes for
%   outfile   - output filename prefix, added to *_oligos.txt and *_seq.txt
%   dnasource - a string value of species to use as a repeat source (see
%     repeatmask.m) Now only use 'human' or 'mouse'
%   noligos   - the desired number of oligos for the probe
%
%   Output Args:
%   outstructure - stucture with fields: score, matches, gcs
%   oligos - structure with fields: seq, rc, comp, position, gc, Tm, seqprop
%
% Example Usage:
%----------------------------------------------------------------
% >> [os,ol] = findprobes_introns_mjl('MKNK2_int.fasta','MKNK2_int','human',16)
% OUTPUT OF REPEATMASK() FUNCTION...
% Searching for matches to human RNAs...
% Searching for matches to miRNAs...
% Searching for matches to human pseudo-genes...
% hit! (pseudogene) 251
% hit! (pseudogene) 252
% hit! (pseudogene) 333
% hit! (pseudogene) 346
% Done searching
% Found a total of 16 oligos...
% >> 
%----------------------------------------------------------------
% 
% 2009-06-01  Modified by Marshall J. Levesque
% This is a modified version of Arjun Raj's "findprobes" function 
% that incorporates searching for matches in the human pseudo-gene 
% database and creates probes for target introns
%=============================================================


%------------------------------------------------------------
% RepeatMask the input sequences using repeatmasker.org
%------------------------------------------------------------

[headers, seqs] = repeatmask(infile,dnasource); % see repeatmask.m

s = multiseq_singlestring(seqs);

RMinseq = s;

RMinseq = lower(RMinseq);  % convert the sequence to lowercase
RMinseq = RMinseq(ismember(RMinseq,'actgxn>'));  % keeps only characters that match a c t g x n or > 


%-------------------------------------------------------------------
% Read the FASTA input file containing the template strand sequences
%-------------------------------------------------------------------

[headers, seqs] = fastaread_blah(infile);

s = multiseq_singlestring(seqs);
    
inseq = s;

inseq = lower(inseq);  % convert the sequence to lowercase
inseq = inseq(ismember(inseq,'actgxn>'));  % keeps only characters that match a c t g x n or > 

%%%%%%%%%%%%%%
% Some basic parameters...
shortlen = 22;
blocklen = 2;

%----------------------------------------------------------------------
% Check for pre-filtered sequence w/ filename: [outfile '_inseq2.txt']
% and use it for designing probes
%----------------------------------------------------------------------

[inseq2,inseq2_old,merge_action] = readinseq2(inseq,outfile,shortlen,blocklen);

if strmatch(merge_action,{'new','append'})
%   ----------------------------
%   Start filtering sequences
%   ----------------------------
%   inseq2 is tested as potential oligo lengths one nucleotide at a time
%   against the filter databases. Any hits result in a masking of that
%   tested oligo length. We use an unedited copy of the inseq2 for creating
%   the test sequences so prior masking doesn't affect subsquent steps.


% Iteration over the sequence uses this 'goodlen' variable, needs to be set
% each time our sequence of interest is changing length. Here, after we
% examine any previously filtered sequence, inseq2 may be only a short
% additional sequence or the full length input. We reset the goodlen value
% below after a potential merge of two sequences.
goodlen = length(inseq2)-shortlen;

% store an unmasked version of the sequence for creating testseq oligos.
% This could be the appending sequence, so doesn't have to be exact copy of
% original input sequence.
inseq2_unedited = inseq2;

fprintf('Searching for matches to human RNAs...\n');
[human_headers, human_seqs] = fastaread_blah('human_filter_reference.fasta');
for i = 1:goodlen
  testseq = inseq2_unedited( (i+2):(i+shortlen-blocklen-1-2) );
  if ~isempty(cell2mat(strfind(human_seqs,testseq)))
    fprintf('hit on human RNAs! %d\n',i);
    inseq2( (i+2):(i+shortlen-blocklen-1-2) ) = 'H';
  end;
  
end;

fprintf('Searching for matches to miRNAs...\n');
mirseqs = load_mirdb('mature.fa','hsa');  % Put hsa for human, mmu for mouse
for i = 1:goodlen
  testseq = inseq2_unedited( (i+2):(i+shortlen-blocklen-1-2) );
  if ~isempty(cell2mat(strfind(mirseqs,testseq)))
    fprintf('\n hit (microRNA), %d!\n',i);
    inseq2( (i+2):(i+shortlen-blocklen-1-2) ) = 'M';
  end;
  
end;

% Indexed values for allowed DNASource map to array of pseudogene filenames
okdnasource = {'human','mouse'};
pseudoDBs   = {'Human57.fasta','Mouse56.fasta'};
tf = strcmpi(dnasource,okdnasource);
if any(tf)
    filename = pseudoDBs{tf};
else
    error('Could not find a pseudogene DB for the specified DNA source %s\n',dnasource);
end
fprintf('Searching for matches to %s pseudo-genes...\n',filename);
[pseudo_headers pseudo_seqs] = fastaread_blah(filename);
pseudo_seqs = lower(pseudo_seqs);
for i = 1:goodlen
  testseq = inseq2_unedited( (i+2):(i+shortlen-blocklen-1-2) );
  if ~isempty(cell2mat(strfind(pseudo_seqs,testseq)))
    fprintf('hit! (pseudogene) %d\n',i);
    % GN - can uncomment the code below to make it not do pseudogenes
    if pseudogenemaskflag 
        inseq2( (i+2):(i+shortlen-blocklen-1-2) ) = 'P' ;
    end
  end;
end;

if pseudogenemaskflag
    fprintf('MASKING PSEUDOGENES \n');
else
    fprintf('NOT MASKING PSEUDOGENES \n');
end

fprintf('Done searching\n');

% MERGE in case of additional sequence in the input FASTA file that needed
% to be filtered and appended to a previously filtered, shorter sequence
if strcmpi(merge_action,'append')
    inseq2 = [inseq2_old inseq2(1+shortlen-blocklen+1:end)];
    fprintf('Additional sequence was filtered and merged w/ stored version\n')
end

end
% END of local filtering steps for 'new' or 'append' merge_action.
% inseq2's value is now either freshly filtered, or reused.
%---------------------------------------------------------------

% Save inseq2 to a file for future reuse or reference
inseq2filename = [outfile '_inseq2.txt'];
fid = fopen(inseq2filename, 'w');
if fid ~= -1
    fprintf(fid,'%s',inseq2);
    fclose(fid);
else
    fprintf('Unable to write inseq2 to file: %s',inseq2filename);
end

% Now let's just mask with the repeat masking we got from before
inseq3 = inseq2;
inseq3(RMinseq=='n') = 'n';
inseq3(RMinseq=='x') = 'x';

inseq = inseq3;  % Now we have a properly masked input sequence

% We reset the goodlen value below after a potential merge of two sequences
goodlen = length(inseq)-shortlen;

%---------------------------------------
% BLAST filter via the NCBI BLAST servers
%---------------------------------------

[counts, seqB, RIDMap] = blastmask(inseq,outfile,dnasource,merge_action,shortlen,blocklen);

if length(inseq) ~= length(seqB)
    error('ERROR: Lengths of inseq and BLAST masked sequences don''t match');
else
    inseq = seqB;
end
    
%---------------------------------------
% Calculate properties of the sequences
%---------------------------------------

for i = 1:goodlen
  testseq = inseq(i:(i+shortlen-blocklen-1));
  testseq(ismember(testseq,'xXmMhHpPbB>')) = 'n';  % Just change this for the oligoprop
  seqprop = oligoprop(testseq,'Salt',0.33,'Temp',30);
  goodness(i) = seqprop.Tm(5);  % Use latest SantaLucia data
  goodness(i) = goodness(i) + 1000*sum(ismember(inseq(i:(i+shortlen-blocklen-1)),'nxXmMhHpPbB>'));
end;
percentgc = goodness;

% Target Tm of 67.5
goodness = (goodness-67.5).^2;

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

scores = [outstructure.score];

inds = find(scores < 100);
noligos = max(inds);

%noligos = length(pos);
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
  testoligo = curroligo{i}.rc;
  seqprop = oligoprop(testoligo,'Salt',0.33,'Temp',30);
  curroligo{i}.Tm = seqprop.Tm(5);
  curroligo{i}.seqprop = seqprop;
  curroligo{i}.gc = seqprop.GC;
end;

oligos = curroligo;

probeseq = blanks(length(inseq));
infoseq  = blanks(length(inseq));

for i = 1:noligos
  probeseq(oligopositions(i):oligopositions(i)+shortlen-blocklen-1) = curroligo{i}.comp;
  numstr = ['Prb# ' num2str(i) ',Tm ' num2str(curroligo{i}.Tm,'%2.3g') ',GC ' num2str(curroligo{i}.gc)];
  infoseq(oligopositions(i):oligopositions(i) + length(numstr)-1) = numstr;
end;

%probeseq
%infoseq

fid = fopen(outoligofile,'w');
for i = 1:noligos
  fprintf(fid,'%d,%d,%2.1f,%s,%s\n',i,round(curroligo{i}.gc),curroligo{i}.Tm,curroligo{i}.rc,[outfile '_' num2str(i)]);
end;
fclose(fid);



linelength = 75;

fid = fopen(outseqfile,'w');
fprintf(fid,'Processed from file %s, oligos in file %s\n\n',infile,outoligofile);
%fprintf(fid,'%d\n\n',length(oligos));

done = 0;
i = 1;
l = length(inseq);
while ~done
  idx = i:min(i+linelength-1,l);
  i;
  fprintf(fid,'%s\n%s\n%s\n\n',inseq(idx),probeseq(idx),infoseq(idx));
  if max(idx)==l
    done = 1;
  end
  i=i+linelength;
end

fclose(fid);


%matches = [];
%
%x = length(goodness);
%while x > 0
%  prevmatch = bmsf(x,1);
%  if ~prevmatch
%    break;
%  end;
%  matches = [matches prevmatch];
%  x = prevmatch-shortlen;
%end;




%%%%%%%%%%%%%%%%%%%%%%%%%%RAW VERSION OF ALGORITHM%%%%%%%%%%%%%%%%%
%
%
%instring = 'atcgggatgctgatctttctgatgctgatgtgcggcgagtttactgatcgatgcgcgagtctg';
%shortlen = 10;
%goodness = rand(length(instring)-shortlen,1);
%
%noligos = 7;
%
%clear bmsf;
%bmsf(1,1).position = 1;
%bmsf(1,1).score = goodness(1);
%for i = 2:noligos
%  bmsf(1,i).position = nan;
%  bmsf(1,i).score = 0;
%end;
%
%
%for x = 2:length(goodness)
%  bmsf = [bmsf ; bmsf(x-1,:)];
%  
%  for k = 1:noligos
%    potentialscore = 0;
%    if k == 1
%      potentialscore = goodness(x);
%    else
%      if x-shortlen > 0
%        if ~isnan(bmsf(x-shortlen,k-1).position)
%          potentialscore = goodness(x) + bmsf(x-shortlen,k-1).score;
%        end;
%      end;
%    end;
%    if potentialscore > bmsf(x,k).score
%      bmsf(x,k).position = x;
%      bmsf(x,k).score = potentialscore;
%    end;
%  end;
%end;
%


%%%%%IMPLEMENTATION OF TYLER'S OLD CODE BELOW%%%%%%%%%%%%%%

%
%
%
%instring = 'atcgggatgctgatctttctgatgctgatgtgcggcgagtttactgatcgatgcgcgagtctg';
%shortlen = 10;
%goodness = rand(length(instring)-shortlen,1);
%
%
%bmsf = [1,goodness(1)];
%
%for x = 2:length(goodness)
%  bmsf = [bmsf ; bmsf(x-1,:)];
%  potential_score = goodness(x);
%  if x-shortlen > 0
%    potential_score = potential_score + bmsf(x-shortlen,2);
%  end;
%  if potential_score > bmsf(x,2)
%    bmsf(x,:) = [x potential_score];
%  end;
%end;
%    
%  
%
%matches = [];
%
%x = length(goodness);
%while x > 0
%  prevmatch = bmsf(x,1);
%  if ~prevmatch
%    break;
%  end;
%  matches = [matches prevmatch];
%  x = prevmatch-shortlen;
%end;
