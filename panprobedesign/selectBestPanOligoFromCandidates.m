function [newMatchVector, whichOligo] = selectBestPanOligoFromCandidates(prevMatchVector,candidateMatchMatrix,whichOligosAreGoodVector,targetMinNumOligos)
% This function will take the prevMatchVector vector of matches and pick the best
% of the candidates to append.

% prevMatchVector is a vector with each entry telling you how many matches the
% previous oligo solution set has to the various target sequence variants.

% candidateMatchMatrix is a e.g. 27x27 matrix of matches.  It is symmetric.
whichOligosAreGoodVector = squeeze(whichOligosAreGoodVector);
candidateMatchMatrix = squeeze(candidateMatchMatrix);

numberOfSeqs = length(candidateMatchMatrix(:,1));

if isempty(prevMatchVector)
    prevMatchVector = zeros(numberOfSeqs,1);
end;

candidateMatchMatrix = squeeze(candidateMatchMatrix);
whichOligosAreGoodVector = squeeze(whichOligosAreGoodVector);

updatedMatchMatrix = zeros(size(candidateMatchMatrix));

for i = 1:numberOfSeqs
    updatedMatchMatrix(:,i) = prevMatchVector + candidateMatchMatrix(:,i);
    scores(i) = computeOligoMatchScore(updatedMatchMatrix(:,i),targetMinNumOligos);
    if ~whichOligosAreGoodVector(i)
        scores(i) = 0;
        updatedMatchMatrix(:,i) = 0;
    end;
end;

[~,idx] = max(scores);
newMatchVector = updatedMatchMatrix(:,idx);
whichOligo = idx;

