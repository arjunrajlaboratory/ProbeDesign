function [output,oligos] = findprobesLocal(infile,varargin)
% findprobesLocal  Create specific oligo probes against target RNA sequence
% findprobesLocal is equivalent to findprobesHD except it uses bowtie
% locally for masking
%
%   Input Args (required):
%   infile    - the template FASTA sequence you are designing probes for
%
%   Options (and defaults):
%
%   output = findprobesLocal('myinfile.fa',34);
%   Generates 34 probes (default is 48)
%
%   output = findprobesLocal('myinfile.fa',34,'outfilename','myout');
%   Puts output in 'myout_seq.txt' and 'myout_seq.txt'
%   Default is to use the input file name (in this case, 'myinfile')
%
%   output = findprobesLocal('myinfile.fa',34,'species','mouse');
%   Sets which species the sequence information comes from for repeat
%   masking, BLASTing and pseudogene masking.  Default is 'human', other
%   options are 'mouse', 'elegans', 'drosophila'.
%
%   output = findprobesLocal('myinfile.fa',34,'repeatmask',true);
%   output = findprobesLocal('myinfile.fa',34,'pseudogenemask',false);
%   output = findprobesLocal('myinfile.fa',34,'genomemask',true);
%   output = findprobesLocal('myinfile.fa',34,'GCrunmask',true);
%   output = findprobesLocal('myinfile.fa',34,'GCmask',true);
%   Tells the software whether to use repeat masking, pseudogene masking,
%   BLAST masking, GCrun masking or GC masking.  All are on by default.
%   Note that GC mask doesn't actually "mask" the sequence, but rather
%   flags bad oligos.
%   
%   output = findprobesLocal('myinfile.fa',34,'targetTM',72);
%   output = findprobesLocal('myinfile.fa',34,'oligolength',25);
%   output = findprobesLocal('myinfile.fa',34,'spacerlength',3);
%   Target TM sets the target optimal TM for the oligos.
%   Oligolength sets the length of the desired oligos.
%   Spacerlength sets the length of the minimum spacer between oligos.
% 
%   output = findprobesLocal('myinfile.fa',32,'masksequences','maskseqs.fa');
%   This makes sure that no oligo in the output binds to any of the targets
%   in maskseqs.fa.  By binding, we mean that any one oligo doesn't have
%   any more than 14 bases of complementarity.
%   You can tweak this number:
%   output = findprobesLocal('myinfile.fa',32,'masksequences','maskseqs.fa','nummatchesformask',13);
%
%   Output Arg:
%   output - stucture with fields: score, matches
%
%----------------------------------------------------------------
if nargin == 0
    findprobesHD_GUI
    return
end

p = inputParser;

p.addRequired('infile',@ischar);
p.addOptional('noligos',48,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));

p.addParamValue('targetTM',67.5,@(x)validateattributes(x,{'numeric'},{'positive'}));
p.addParamValue('spacerlength',2,@(x)validateattributes(x,{'numeric'},{'nonnegative'}));
p.addParamValue('oligolength',20,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));

[~,tmp,~] = fileparts(infile);
p.addParamValue('outfilename',tmp,@ischar);

otherSpecies = {'human','mouse','drosophila','elegans','rat','cow','off'};
possibleThermoParams= {'RNA','DNA'};

p.addParamValue('species','human',@(x)any(strcmpi(x,otherSpecies)));
p.addParamValue('repeatmask',true,@islogical);
p.addParamValue('repeatmaskmanual',false,@islogical);
p.addParamValue('repeatmaskfile','',@ischar);
p.addParamValue('pseudogenemask',true,@islogical);
p.addParamValue('blastmask',true,@islogical);%unused
p.addParamValue('miRNAmask',true,@islogical); %unused
p.addParamValue('humanfiltermask',true,@islogical);%unused
p.addParamValue('genomemask',true,@islogical);
p.addParamValue('GCrunmask',false,@islogical);
p.addParamValue('GCmask',false,@islogical);
%BE 09/19/2022 - adding parameter below to filter probes with %GC outside
%bounds
p.addParameter('GCbounds', [], @(x)validateattributes(x,{'numeric'},{'>=', 0, '<=', 1, 'numel', 2}))

p.addParamValue('thermoparams','RNA',@(x)any(strcmpi(x,possibleThermoParams)));
p.addParamValue('targetGibbsFE',-23,@(x)validateattributes(x,{'numeric'},{'<', 0}));
p.addParamValue('allowableGibbsFE',[-26,-20],@(x)validateattributes(x,{'numeric'},{'size',[1 2]}));

p.addParamValue('masksequences','',@ischar);
p.addParamValue('nummatchesformask',14,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
p.addParamValue('writeoutput',true,@islogical);

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

%BE added to handle string inputs
if isempty(dir(infile))
    outfile = headers;
end
%%%%%%%%%%%%%%
% Some basic parameters...
oligolen = p.Results.oligolength; %shortlen = 22;
blocklen = p.Results.spacerlength;                         %blocklen = 2;

goodlen = length(inseq)-oligolen;

fullmask = zeros(length(inseq),1);
maskseqs = {};

%------------------------------------------------------------
% RepeatMask the input sequences using repeatmasker.org
%------------------------------------------------------------


if p.Results.repeatmask
    if p.Results.repeatmaskmanual
        [headers, seqs] = fastaread_blah(p.Results.repeatmaskfile); % feed in repeat mask file manually
    else
        [headers, seqs] = repeatmask(infile,dnasource); % see repeatmask.m
    end
    
    
    s = multiseq_singlestring(seqs);
    
    RMinseq = s;
    
    RMinseq = lower(RMinseq);  % convert the sequence to lowercase
    RMinseq = RMinseq(ismember(RMinseq,'actgxn>'));  % keeps only characters that match a c t g x n or >
    hitsRM = RMinseq == 'n';
    fullmask = fullmask + hitsRM';
    maskseqs = [maskseqs mask_string(inseq,hitsRM,'R')];
    
end

if ischar(infile)
    fileString = text_to_string_local(infile);
else
    fileString = ['>' infile.Header newline infile.Sequence];
end

if p.Results.pseudogenemask
    % Indexed values for allowed DNASource map to array of pseudogene filenames
    okdnasource = {'human','mouse','drosophila','elegans','rat'};
    tf = strcmpi(dnasource,okdnasource);
    if any(tf)
        db = [okdnasource{tf} 'Pseudo'];
    else
        error('Could not find a pseudogene DB for provided species: %s\n',dnasource);
    end
    
    if strcmpi(db,'elegansPseudo')  % This is a kludge--bowtie requires celegans, not elegans
        db = 'celegansPseudo';
    end
    
    hitsPseudo = hits_to_mask_local(fileString,16,db,0);
    % These parameters correspond to the exact same parameters we have used
    % before.  They seemed to work pretty well then (at least a few
    % documented saves), so we'll just go with it as is.
    
    % One thing is that you will end up with a lot of little 16-18 mer
    % masks.  I'm not convinced these are really important, but they do
    % sometimes block off so much sequence that it can be problematic.  So
    % we'll get rid of them this way:
    hitsPseudo = remove_short_runs(hitsPseudo,20,2);
    
    fullmask = fullmask + hitsPseudo;
    maskseqs = [maskseqs mask_string(inseq,hitsPseudo,'P')];

    fprintf('Done sequence searching, finding probe positions...\n');
end


if p.Results.genomemask
    % Indexed values for allowed DNASource map to array of pseudogene filenames
    okdnasource = {'human','mouse','drosophila','elegans','rat','cow'};
    tf = strcmpi(dnasource,okdnasource);
    if any(tf)
        db = [okdnasource{tf}];
    else
        error('Could not find a genome DB for provided species: %s\n',dnasource);
    end
    
    if strcmpi(db,'elegans')  % This is a kludge--bowtie requires celegans, not elegans
        db = 'celegans';
    end
    
    hitsGen1 = hits_to_mask_local(fileString,12,db,4000);
    hitsGen2 = hits_to_mask_local(fileString,14,db,500);
    hitsGen3 = hits_to_mask_local(fileString,16,db,20);
    
    hitsGenome = hitsGen1|hitsGen2|hitsGen3;
    
    % The idea here is to combine hits/thresholds for a number of different
    % seed lengths.  Of course, as you decrease the seed length, you will
    % get more hits, so the threshold will have to be higher.  In the end,
    % I arrived at these numbers by checking for "spikes" in HOTAIR and
    % LGR5 and also by looking to eliminate oligos in the shark tank.
    
    fullmask = fullmask + hitsGenome;
    
    maskseqs = [maskseqs mask_string(inseq,hitsGenome,'B')];
    
end

GCrunmask = zeros(length(inseq),1);
if p.Results.GCrunmask
    % Mask runs of Cs, Gs...
    %Cmask = mask_runs(inseq,'c',7,2);  % This is to mask out all potential runs of 'g' in the oligos
    %Gmask = mask_runs(inseq,'g',7,2);  % This is to mask out all potential runs of 'c' in the oligos
    %maskseqs = [maskseqs mask_string(inseq,Cmask|Gmask,'X')];
    %fullmask = fullmask + Cmask + Gmask;
    
    % NOTE: This is not really a mask. It gives a score per oligo.
                                        % In order to use this, it really
                                        % should be added to badness
                                        % directly without extending the
                                        % masking to the length of the
                                        % oligo.  The return argument
                                        % GCmask should be like this, but
                                        % just in case you're wondering why
                                        % X locations don't "make sense"...   
    Crunbadness = mask_oligos_with_runs(inseq,'c',7,2,oligolen);
    Grunbadness = mask_oligos_with_runs(inseq,'g',7,2,oligolen);
    maskseqs = [maskseqs mask_string(inseq,Crunbadness|Grunbadness,'X')];
    
    GCrunmask(Crunbadness>0) = inf;
    GCrunmask(Grunbadness>0) = inf;
end

if p.Results.GCmask
    % Mask "unbalanced oligos" with too many Gs, Cs.
    GCmask = GC_badness(inseq,oligolen); % NOTE: This is not really a mask. It gives a score per oligo.
                                        % In order to use this, it really
                                        % should be added to badness
                                        % directly without extending the
                                        % masking to the length of the
                                        % oligo.  The return argument
                                        % GCmask should be like this, but
                                        % just in case you're wondering why
                                        % there are only sporadic Zs...
    maskseqs = [maskseqs mask_string(inseq,GCmask,'Z')];
    %fullmask = fullmask + GCmask; %NOTE: This works differently, so you
    %want to put this directly into the badness, not as a mask.
else
    GCmask = zeros(length(inseq),1);
end
GCmask(GCmask>0) = inf;

if ~isempty(p.Results.masksequences)
    % Mask out the oligos that are too similar to the other sequences in
    % p.Results.masksequences
    
    maxmatches = mask_by_other_sequences(inseq,p.Results.masksequences,oligolen);
    other_seq_mask = zeros(size(maxmatches));
    other_seq_mask(maxmatches > p.Results.nummatchesformask) = Inf;
    maskseqs = [maskseqs mask_string(inseq,other_seq_mask==Inf,'M')];
else
    other_seq_mask = zeros(1,length(inseq));
end

%BE 09/19/2022 - Added new mask
if ~isempty(p.Results.GCbounds)
    % Mask oligos with %GC outside bounds
    % As with GCrunmask and GCmask
    % This is not really a mask. It gives a score per oligo.
    GCtotalmask = GCtotal_badness(inseq,oligolen,p.Results.GCbounds); 
    
    maskseqs = [maskseqs mask_string(inseq,GCtotalmask,'Z')];

else
    GCtotalmask = zeros(length(inseq),1);
end
GCtotalmask(GCtotalmask>0) = inf;

% Done with all masking
if strcmp(p.Results.thermoparams,'DNA')
    badness = Tm_badness(inseq,oligolen,p.Results.targetTM);
else
    badness = Tm_badness_RNA_DNA(inseq,oligolen,p.Results.targetGibbsFE,p.Results.allowableGibbsFE);
    maskseqs = [maskseqs mask_string(inseq,badness==inf,'F')];
end

badness = badness + mask_to_badness(fullmask,oligolen); % Want to fully block out all the truly masked sequences
badness = badness + GCmask + GCrunmask + other_seq_mask' + GCtotalmask; % In this case, we just block out these particular oligos.

output = find_best_matches(badness,goodlen,oligolen+blocklen,noligos);

noligos = length(output);

if noligos==0
    fprintf('\n************* COULD NOT FIND ANY OLIGOS for %s *************\n\n',outfile)
else
    currmatches = output(end).matches;
    
    
    outoligofile = [outfile '_oligos.txt'];
    outseqfile = [outfile '_seq.txt'];
    
    if strcmp(p.Results.thermoparams,'DNA')
        oligos = generate_oligo_file(inseq,oligolen,currmatches,outfile,outoligofile);
        align_oligo_file(inseq,maskseqs,oligolen,currmatches,outfile,outseqfile);
    else
        oligos = generate_oligo_file_RNA_DNA(inseq,oligolen,currmatches,outfile,outoligofile);
        align_oligo_file_RNA_DNA(inseq,maskseqs,oligolen,currmatches,outfile,outseqfile);
    end
    
    
    fprintf('Done. Fo real.\n');
end