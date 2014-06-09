function outseq = mask_runs(inseq,chr,runlength,mismatches)

inseq = lower(inseq);

%outseq = inseq;
outseq = zeros(length(inseq),1);

for i = 1:(length(inseq)-runlength+1)
    tmp = inseq(i:i+runlength-1) == chr;
    if sum(tmp) >= runlength-mismatches
        outseq(i:i+runlength-1) = 1;
    end;
end;
