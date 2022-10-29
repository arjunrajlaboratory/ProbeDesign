
function badness = GCtotal_badness(inseq,oligolen, GCbounds)

badness = zeros(length(inseq),1);

goodlen = length(inseq)-oligolen;

for i = 1:goodlen+1
  testseq = inseq(i:(i+oligolen-1));
  testseq(ismember(testseq,'xXmMhHpPbB>')) = 'n';  % Just change this for the oligoprop
  badness(i) = bad_gc(testseq, GCbounds);
end

