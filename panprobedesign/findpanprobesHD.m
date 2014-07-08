function [bestSolution,oligos] = findpanprobesHD(inputAlignmentFile,varargin)

p = inputParser;

p.addRequired('inputAlignmentFile',@ischar);
p.addOptional('maxNumOligos',64,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
p.addOptional('targetMinNumOligos',30,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));

%p.addParamValue('targetTM',67.5,@(x)validateattributes(x,{'numeric'},{'positive'}));
p.addParamValue('spacerLength',2,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
p.addParamValue('oligoLength',20,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));

[~,fileHeader,~] = fileparts(inputAlignmentFile);
p.addParamValue('outputFileHeader',fileHeader,@ischar);

otherSpecies = {'human','mouse','drosophila','elegans','rat','cow','off'};
%possibleThermoParams = {'RNA','DNA'};

p.addParamValue('species','human',@(x)any(strcmpi(x,otherSpecies)));
%p.addParamValue('repeatmask',true,@islogical);
%p.addParamValue('pseudogenemask',true,@islogical);
%p.addParamValue('blastmask',true,@islogical);
%p.addParamValue('miRNAmask',true,@islogical);
%p.addParamValue('humanfiltermask',true,@islogical);
%p.addParamValue('genomemask',true,@islogical);
%p.addParamValue('GCrunmask',false,@islogical);
%p.addParamValue('GCmask',false,@islogical);

%p.addParamValue('thermoparams','RNA',@(x)any(strcmpi(x,possibleThermoParams)));
p.addParamValue('targetGibbsFE',-23,@(x)validateattributes(x,{'numeric'},{'negative'}));
p.addParamValue('allowableGibbsFE',[-26,-20],@(x)validateattributes(x,{'numeric'},{'size',[1 2]}));

%p.addParamValue('masksequences','',@ischar);
%p.addParamValue('nummatchesformask',14,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));


p.parse(inputAlignmentFile,varargin{:});

outputFileHeader = p.Results.outputFileHeader;
dnasource = p.Results.species;
maxNumOligos = p.Results.maxNumOligos;
targetMinNumOligos = p.Results.targetMinNumOligos;

% Some basic parameters...
oligoLength = p.Results.oligoLength;
spacerLength = p.Results.spacerLength;


% done parsing input arguments, begin doing stuff.

alignedFasta = fastaread(inputAlignmentFile);
cleanedFasta = cleanAlignedFastas(alignedFasta); % changes all non ACTG to n and lowers case.
numberOfSeqs = numel(cleanedFasta);

targetableLength = length(cleanedFasta(1).Sequence) - oligoLength;

sequenceCrossMatches = getSequenceCrossMatches(cleanedFasta,oligoLength);

for i = 1:numberOfSeqs
    sequenceScores(i,:) = Tm_badness_RNA_DNA(cleanedFasta(i).Sequence,oligoLength,p.Results.targetGibbsFE,p.Results.allowableGibbsFE);
    whichOligosAreGood(i,:) = ~isinf(sequenceScores(i,:));
end


% Here comes the algorithm stuff:

bmsf(1,1).position = 1;
[bmsf(1,1).matchVector,bmsf(1,1).whichOligo] = ...
    selectBestPanOligoFromCandidates([],sequenceCrossMatches(1,:,:),whichOligosAreGood(1,:),targetMinNumOligos);
for i = 2:maxNumOligos  % For position 1, set score for all higher numbers of oligos to nan
    bmsf(1,i).position = nan;
    bmsf(1,i).matchVector = zeros(size(bmsf(1,1).matchVector));
end;



for x = 2:targetableLength   % length(allinds(:,1,1))
    fprintf('Finding optimal solution at position %d\n',x);
    bmsf(x,:) = bmsf(x-1,:); % By default, previous solution is best at next step
    
    for k = 1:maxNumOligos
        potentialMatches.matchVector = zeros(size(bmsf(1,1).matchVector));
        if k == 1
            [potentialMatches.matchVector,potentialMatches.whichOligo] = ...
                selectBestPanOligoFromCandidates([],sequenceCrossMatches(x,:,:),whichOligosAreGood(:,x),targetMinNumOligos);
        else
            if x-(oligoLength+spacerLength) > 0
                if ~isnan(bmsf(x-(oligoLength+spacerLength),k-1).position)
                    [potentialMatches.matchVector,potentialMatches.whichOligo] = ...
                        selectBestPanOligoFromCandidates(bmsf(x-(oligoLength+spacerLength),k-1).matchVector,sequenceCrossMatches(x,:,:),whichOligosAreGood(:,x),targetMinNumOligos);
                end;
            end;
        end;
        if computeOligoMatchScore(potentialMatches.matchVector,targetMinNumOligos) > computeOligoMatchScore(bmsf(x,k).matchVector,targetMinNumOligos)
            bmsf(x,k).position = x;
            bmsf(x,k).matchVector = potentialMatches.matchVector;
            bmsf(x,k).whichOligo = potentialMatches.whichOligo;
        end;
    end;
end;

% Now read out the solution:
for k = 1:maxNumOligos
    allMatches = [];
    x = targetableLength;
    currk = k;
    if ~isnan(bmsf(x,currk).position)
        bestSolution(k).finalScoreInfo = bmsf(x,currk);
        for counter = 1:k
            prevmatch.position = bmsf(x,currk).position;
            prevmatch.whichOligo = bmsf(x,currk).whichOligo;
            allMatches = [prevmatch allMatches];
            x = prevmatch.position-(oligoLength+spacerLength);
            currk = currk - 1;
        end;
        bestSolution(k).allMatches = allMatches;
    end
end



for i = 1:length(bestSolution)
    scores(i) = computeOligoMatchScore(bestSolution(i).finalScoreInfo.matchVector,targetMinNumOligos);
end;

inds = find(scores < 1000000);
finalNumOligos = max(inds);

fprintf('Found a total of %d oligos...\n',finalNumOligos);

% Now let's output the files...

outoligofile = [outputFileHeader '_oligos.txt'];
outmaskedseqfile = [outputFileHeader '_masked_seq.txt'];
outseqfile = [outputFileHeader '_seq.txt'];

oligos = generateOligoFilePanProbes(alignedFasta,oligoLength,bestSolution(finalNumOligos).allMatches,outputFileHeader,outoligofile);

for i = 1:numberOfSeqs
    maskedSequences(i).Sequence = mask_string(alignedFasta(i).Sequence,~whichOligosAreGood(i,:),'x');
end

generateAlignedFilePanProbe(maskedSequences,oligos,oligoLength,bestSolution(finalNumOligos).allMatches,outputFileHeader,outmaskedseqfile);
generateAlignedFilePanProbe(alignedFasta,oligos,oligoLength,bestSolution(finalNumOligos).allMatches,outputFileHeader,outseqfile);




