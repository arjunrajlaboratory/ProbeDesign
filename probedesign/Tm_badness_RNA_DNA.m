
function badness = Tm_badness_RNA_DNA(inseq,oligolen,targetGibbsFE,allowableGibbsFE)

badness = zeros(length(inseq),1);

goodlen = length(inseq)-oligolen;

allowableRange = sort(allowableGibbsFE,'ascend');  % Let's just sort this to make sure the first number is smaller

for i = 1:goodlen
  testseq = inseq(i:(i+oligolen-1));
  if sum(ismember(testseq,'xXmMhHpPbBnN>'))
      badness(i) = inf;
  else
    badness(i) = getGibbs_RNA_DNA(testseq);
    if (badness(i) < allowableRange(1)) || (badness(i) > allowableRange(2))
        %fprintf('%d %g \n',i,badness(i));
        badness(i) = inf;
    end
  end
end

badness = (badness-targetGibbsFE).^2;

