function sequenceCrossMatches = getSequenceCrossMatches(alignedFasta,oligoLength)

targetableLength = length(alignedFasta(1).Sequence) - oligoLength;
numberOfSeqs = numel(alignedFasta);

sequenceCrossMatches = zeros(targetableLength,numberOfSeqs,numberOfSeqs);

for i = 1:targetableLength
    for j = 1:numberOfSeqs
        testSeqs{j} = alignedFasta(j).Sequence(i:(i+oligoLength-1));
    end;
    for j = 1:numberOfSeqs
        sequenceCrossMatches(i,j,:) = strcmp(testSeqs{j},testSeqs);    
    end;
end;