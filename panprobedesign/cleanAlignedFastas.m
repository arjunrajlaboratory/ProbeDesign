function alignedFasta = cleanAlignedFastas(in_alignedFasta)

alignedFasta = in_alignedFasta; % Just copy over to output variable to copy structure.
for i = 1:numel(in_alignedFasta)
    tempSeq = in_alignedFasta(i).Sequence;
    tempSeq(~ismember(lower(tempSeq),'actg')) = 'n';
    alignedFasta(i).Sequence = lower(tempSeq);
end;
