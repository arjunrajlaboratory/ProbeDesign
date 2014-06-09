function [outstructure,oligos] = findprobesHD(infile,varargin)
% findprobesHD  Create specific oligo probes against target RNA sequence
%
%   Input Args (required):
%   infile    - the template FASTA sequence you are designing probes for
%
%   Options (and defaults):
%
%   [os,ol] = findprobesHD('myinfile.fa',34);
%   Generates 34 probes (default is 48)
%
%   [os,ol] = findprobesHD('myinfile.fa',34,'outfilename','myout');
%   Puts output in 'myout_seq.txt' and 'myout_seq.txt'
%   Default is to use the input file name (in this case, 'myinfile')
%
%   [os,ol] = findprobesHD('myinfile.fa',34,'species','mouse');
%   Sets which species the sequence information comes from for repeat
%   masking, BLASTing and pseudogene masking.  Default is 'human', other
%   option is 'mouse'.
%
%   [os,ol] = findprobesHD('myinfile.fa',34,'miRNAmask',true);
%   Masks micro-RNA matching sequences. Species with miRNA data 
%   available are 'hsa','mmu','dme','cel'. Uses a file in the probe
%   design directory called 'mature.fa' (default is true).
%
%   [os,ol] = findprobesHD('myinfile.fa',34,'repeatmask',true);
%   [os,ol] = findprobesHD('myinfile.fa',34,'pseudogenemask',false);
%   [os,ol] = findprobesHD('myinfile.fa',34,'blastmask',true);
%   [os,ol] = findprobesHD('myinfile.fa',34,'humanfiltermask',false);
%   Tells the software whether to use repeat masking, pseudogene masking,
%   BLAST masking or human RNA filtering.  All are on by default.
%   
%   [os,ol] = findprobesHD('myinfile.fa',34,'targetTM',72);
%   [os,ol] = findprobesHD('myinfile.fa',34,'oligolength',25);
%   [os,ol] = findprobesHD('myinfile.fa',34,'spacerlength',3);
%   Target TM sets the target optimal TM for the oligos.
%   Oligolength sets the length of the desired oligos.
%   Spacerlength sets the length of the minimum spacer between oligos.
%
%
%   Output Args:
%   outstructure - stucture with fields: score, matches, gcs
%   oligos - structure with fields: seq, rc, comp, position, gc, Tm, seqprop
%
%----------------------------------------------------------------


%-----------------------------
% FIRST Parse input arguments 
%-----------------------------

p = inputParser;

p.addRequired('infile',@ischar);
p.addOptional('noligos',48,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));

p.addParamValue('targetTM',67.5,@(x)validateattributes(x,{'numeric'},{'positive'}));
p.addParamValue('spacerlength',2,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
p.addParamValue('oligolength',20,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));

[~,tmp,~] = fileparts(infile);
p.addParamValue('outfilename',tmp,@ischar);

otherSpecies = {'human','mouse','drosophila','elegans','off'};

p.addParamValue('species','human',@(x)any(strcmpi(x,otherSpecies)));
p.addParamValue('repeatmask',true,@islogical);
p.addParamValue('pseudogenemask',true,@islogical);
p.addParamValue('blastmask',true,@islogical);
p.addParamValue('miRNAmask','true',@islogical);
p.addParamValue('humanfiltermask',true,@islogical);

p.parse(infile,varargin{:});

outfile = p.Results.outfilename;
dnasource = p.Results.species;
noligos = p.Results.noligos;

%%%%%%%%%% DONE parsing input %%%%%%%%%%%%%%%%%%


%-------------------------------------------------------------------
% Read the FASTA input file containing the template strand sequences
%-------------------------------------------------------------------

[headers, seqs] = fastaread_blah(infile);

s = multiseq_singlestring(seqs);
    
inseq = s;

inseq = lower(inseq);  % convert the sequence to lowercase
inseq = inseq(ismember(inseq,'actgxn>'));  % keeps only characters that match a c t g x n or > 


%------------------------------------------------------------
% RepeatMask the input sequences using repeatmasker.org
%------------------------------------------------------------


if p.Results.repeatmask
    [headers, seqs] = repeatmask(infile,dnasource); % see repeatmask.m
    
    s = multiseq_singlestring(seqs);
    
    RMinseq = s;
    
    RMinseq = lower(RMinseq);  % convert the sequence to lowercase
    RMinseq = RMinseq(ismember(RMinseq,'actgxn>'));  % keeps only characters that match a c t g x n or >
else
    RMinseq = inseq;
end;


%%%%%%%%%%%%%%
% Some basic parameters...
shortlen = p.Results.oligolength + p.Results.spacerlength; %shortlen = 22;
blocklen = p.Results.spacerlength;                         %blocklen = 2;

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

% AR Edit out 8/27: goodlen = length(inseq2)-shortlen;
goodlen = length(inseq2)-(shortlen-blocklen);

% store an unmasked version of the sequence for creating testseq oligos.
% This could be the appending sequence, so doesn't have to be exact copy of
% original input sequence.
inseq2_unedited = inseq2;

if p.Results.humanfiltermask && strcmp(dnasource,'human')
    fprintf('Searching for matches to human RNAs...\n');
    [human_headers, human_seqs] = fastaread_blah('human_filter_reference.fasta');
    for i = 1:goodlen
        testseq = inseq2_unedited( (i+2):(i+shortlen-blocklen-1-2) );
        if ~isempty(cell2mat(strfind(human_seqs,testseq)))
            fprintf('hit on human RNAs! %d\n',i);
            inseq2( (i+2):(i+shortlen-blocklen-1-2) ) = 'H';
        end;
        
    end;
end;

if p.Results.miRNAmask
    speciesInd = strcmpi(dnasource,{'human','mouse','drosophila','elegans'});
    mirSpecies = {'hsa','mmu','dme','cel'};
    if any(speciesInd)
        mirSpeciesStr = cell2mat(mirSpecies(speciesInd));
    else
        error('Could not find a miRNA DB for the provided species: %s\n',dnasource);
    end
    fprintf(1,'Searching for matches to miRNAs for %s...\n',dnasource);
    mirseqs = load_mirdb('mature.fa',mirSpeciesStr);
    for i = 1:goodlen
        testseq = inseq2_unedited( (i+2):(i+shortlen-blocklen-1-2) );
        if ~isempty(cell2mat(strfind(mirseqs,testseq)))
            fprintf('\n hit (microRNA), %d!\n',i);
            inseq2( (i+2):(i+shortlen-blocklen-1-2) ) = 'M';
        end;
    end;
end;


if p.Results.pseudogenemask
    % Indexed values for allowed DNASource map to array of pseudogene filenames
    okdnasource = {'human','mouse','drosophila','elegans'};
    tf = strcmpi(dnasource,okdnasource);
    if any(tf)
        filename = ['pseudogeneDBs/' okdnasource{tf} '.fasta'];
    else
        error('Could not find a pseudogene DB for provided species: %s\n',dnasource);
    end
    fprintf('Searching for matches to %s pseudogenes...\n',filename);
    [pseudo_headers pseudo_seqs] = fastaread_blah(filename);
    pseudo_seqs = lower(pseudo_seqs);
    for i = 1:goodlen
        testseq = inseq2_unedited( (i+2):(i+shortlen-blocklen-1-2) );
        if ~isempty(cell2mat(strfind(pseudo_seqs,testseq)))
            fprintf('hit! (pseudogene) %d\n',i);
            inseq2( (i+2):(i+shortlen-blocklen-1-2) ) = 'P' ;
        end;
    end;
    
    fprintf('Done searching\n');
end;

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
goodlen = length(inseq)-(shortlen-blocklen);

%---------------------------------------
% BLAST filter via the NCBI BLAST servers
%---------------------------------------
if p.Results.blastmask
    okdnasource = {'human','mouse','drosophila','elegans'};
    if any(strcmpi(dnasource,okdnasource))
        [counts, seqB, RIDMap] = blastmask(inseq,outfile,dnasource,merge_action,shortlen,blocklen);
    else
        error('BLAST masking is not currently supported for species %s',dnasource);
    end
    
    if length(inseq) ~= length(seqB)
        error('ERROR: Lengths of inseq and BLAST masked sequences don''t match');
    else
        inseq = seqB;
    end
end;

%---------------------------------------
% Calculate properties of the sequences
%---------------------------------------

for i = 1:goodlen
  testseq = inseq(i:(i+shortlen-blocklen-1));
  testseq(ismember(testseq,'xXmMhHpPbB>')) = 'n';  % Just change this for the oligoprop
  seqprop = oligoprop(testseq,'Salt',0.33,'Temp',30);
  goodness(i) = seqprop.Tm(5);  % Use latest SantaLucia data
  goodness(i) = goodness(i) + 10000*sum(ismember(inseq(i:(i+shortlen-blocklen-1)),'nxXmMhHpPbB>'));
end;
percentgc = goodness;

% Target Tm of 67.5
goodness = (goodness-p.Results.targetTM).^2;
%goodness = (goodness-67.5).^2;

clear bmsf;
bmsf(1,1).position = 1;
bmsf(1,1).score = goodness(1);
for i = 2:noligos
  bmsf(1,i).position = nan;
  bmsf(1,i).score = inf;
end;


for x = 2:goodlen %length(goodness)
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
