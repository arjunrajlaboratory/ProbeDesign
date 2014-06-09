
function badness = Tm_badness(inseq,oligolen,targetTM)

badness = zeros(length(inseq),1);

goodlen = length(inseq)-oligolen;

for i = 1:goodlen
  testseq = inseq(i:(i+oligolen-1));
  if sum(ismember(testseq,'xXmMhHpPbBnN>'))
      badness(i) = inf;
  else
    testseq(ismember(testseq,'xXmMhHpPbB>')) = 'n';  % Just change this for the oligoprop
    badness(i) = getTm(testseq);
  end;
end;

badness = (badness-targetTM).^2;

