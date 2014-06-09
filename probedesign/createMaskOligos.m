%% createMaskOligos
% Produces reverse complement sequences for a RNA FISH probe's oligos

%% Input
% *|probeOligosFilename|* - _oligos.txt file that contains FISH probe sequences 
% *|seqLength|* - length of the masking oligos (in base pairs)
% *|seqPosition|* - string equal to one of {'5prime','center','3prime'}
% *|oligoSubSet|* - provide an array of probe oligo numbers to limit what you want masked

%% Example Usage
% Create a set of 16bp masking oligos that hybridize to the 3-prime end of 
% SUZ12 mRNA probe's odd numbered oligos
%
%  >> createMaskOligos('SUZ12_mRNA_oligos.txt',16,'3prime','oligoSubSet',[1:2:15])
%  Creating 16-mer oligos positioned 3prime to mask the probe in SUZ12_mRNA_oligos.txt
%  ctccattttcggcttcttca SUZ12_mRNA_1
%      ||||||||||||||||
%      taaaagccgaagaagt
%  
%  gatctgtgttggcttctcaa SUZ12_mRNA_3
%      ||||||||||||||||
%      acacaaccgaagagtt
%  
%  tgtgcaaaaatattggtgct SUZ12_mRNA_5
%      ||||||||||||||||
%      gtttttataaccacga
%  
%  tgagattcttgctctccttt SUZ12_mRNA_7
%      ||||||||||||||||
%      taagaacgagaggaaa
%  
%  tggaagaaaccagtaaacgt SUZ12_mRNA_9
%      ||||||||||||||||
%      tctttggtcatttgca
%  
%  gcaataggagccgtagattt SUZ12_mRNA_11
%      ||||||||||||||||
%      atcctcggcatctaaa
%  
%  taactgaaccaggcttgttt SUZ12_mRNA_13
%      ||||||||||||||||
%      acttggtccgaacaaa
%  
%  tgttgccttgtattgttgtt SUZ12_mRNA_15
%      ||||||||||||||||
%      cggaacataacaacaa
%  
%  Wrote mask sequences to SUZ12_mRNA_16bp-3prime_mask.txt

%% Author
% Marshall J. Levesque 2011

function createMaskOligos(probeOligosFilename,seqLength,seqPosition,varargin)

%--------------------------------
% Make sure we have valid input
%--------------------------------

p = inputParser;
p.addRequired('probeOligosFilename', @(x)strfind(x,'_oligos.txt'));
p.addRequired('seqLength',@(x)isnumeric(x) && x>=1);
p.addOptional('seqPosition','center',@(x)any(strcmpi(x,{'5prime','center','3prime'})));
p.addOptional('oligoSubSet',[],@isnumeric);
p.parse(probeOligosFilename,seqLength,seqPosition,varargin{:});

seqLength = p.Results.seqLength;
seqPosition = p.Results.seqPosition; 
probeOligosFilename = p.Results.probeOligosFilename;
oligoSubSet = p.Results.oligoSubSet;
fprintf(1,'Creating %d-mer oligos positioned %s to mask the probe in %s\n',...
            p.Results.seqLength,p.Results.seqPosition,p.Results.probeOligosFilename);


%------------------------------
% Read in the probe oligos file
%-------------------------------
fastaStr = probeOligosToFASTA(probeOligosFilename);
seqs     = fastaread_blah(fastaStr);

if isempty(oligoSubSet)  % user provided no input, mask all sequences of this probe
    oligoSubSet = 1:length(seqs);
end


%----------------------------------------------------------------------
% For each oligo, generate reverse complement at the specified position
%----------------------------------------------------------------------

% check whether specified masking is longer than probe oligo length
if seqLength > length(seqs(1).Sequence)
    error('Length of masking oligos must be less than or equal to probe oligo length');
end


% Create output filename: SUZ12_mRNA_16bp-3prime_mask.txt
newSuffix = ['_' num2str(seqLength) 'bp-' seqPosition '_mask.txt'];
outputFileName = strrep(probeOligosFilename,'_oligos.txt',newSuffix);
fOut = fopen(outputFileName,'w');
fprintf(fOut,'MaskSequence\tSequenceID\n');

% iterate over the sequences in the probe oligo file
for i = 1:length(seqs)
    % check for any specified probe oligo sequence limits, skip if not specified
    if ~any(i==oligoSubSet); continue; end; % DEFAULT behavior is use all oligos

    seq = seqs(i);
    probeOligo = seq.Sequence;  % grab the sequence

    % cut the sequence to specified length and position
    whiteSpaceL = ''; whiteSpaceR = '';
    if strcmp(seqPosition,'5prime')
        maskOligo = probeOligo(1:seqLength);
        whiteSpaceR = repmat(' ',1,length(probeOligo) - seqLength);
    elseif strcmp(seqPosition,'center')
        % handle cases of odd or even number seqLength
        seqPadding = (length(probeOligo) - seqLength) / 2;
        startIndex = 1+ceil(seqPadding);
        endIndex   = length(probeOligo)-floor(seqPadding);
        maskOligo = probeOligo(startIndex:endIndex);
        whiteSpaceR = repmat(' ',1,ceil(seqPadding));
        whiteSpaceL = repmat(' ',1,floor(seqPadding));
    elseif strcmp(seqPosition,'3prime')
        startIndex = length(probeOligo) - seqLength +1;
        maskOligo = probeOligo(startIndex:end);
        whiteSpaceL = repmat(' ',1,length(probeOligo) - seqLength);
    end

    % format strings with white space padding
    probeSeqNHeader = [probeOligo ' ' seq.Header];
    basePairing     = [whiteSpaceL repmat('|',1,seqLength) whiteSpaceR];
    formatMaskOligo = [whiteSpaceL comp(maskOligo) whiteSpaceR];

    % Print to screen, example:
    % >> createMaskOligos('SUZ12_mRNA_oligos.txt',12,'5prime'):
    % ggtgggaatcaccaactttt  SUZ12_mRNA_1
    % ||||||||||||        
    % ccacccttagtg 
    % ...
    fprintf(1,'%s\n%s\n%s\n\n',probeSeqNHeader,basePairing,formatMaskOligo);

    % Print to output file, example file contents:
    % MaskSequence    SequenceID
    % gaaaatggag      SUZ12_mRNA_1_10bp-5prime_mask
    % gcttttcctc      SUZ12_mRNA_2_10bp-5prime_mask
    % ...
    fprintf(fOut,'%s\t%s\n',revcomp(maskOligo),[seq.Header strrep(newSuffix,'.txt','')]);
    
end
    
fprintf(1,'Wrote mask sequences to %s\n',outputFileName);
fclose(fOut);
    

