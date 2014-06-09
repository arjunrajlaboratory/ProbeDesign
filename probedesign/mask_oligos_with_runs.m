function outseq = mask_oligos_with_runs(inseq,chr,runlength,mismatches,oligolen)

inseq = lower(inseq);
outseq = zeros(size(inseq));

for i = 1:(length(inseq)-oligolen)
    tmp = inseq(i:i+oligolen-1);
    rn = mask_runs(tmp,chr,runlength,mismatches);
    if sum(rn) > 0
        outseq(i) = 1;
    end;
end;

