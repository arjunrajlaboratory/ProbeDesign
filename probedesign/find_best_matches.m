function output = find_best_matches(badness, goodlen, probe_spacer_len, noligos)
% This function finds the matches (up to noligos) that minimize
% average "badness".  goodlen should be:
% the length of the sequence - oligo length
% Technically, goodlen should be this length, but we'll make the
% parameter explicit just in case.
% 
% If noligos is too large, it will find as many oligos as possible.
%
% Output:

bmsf(1,1).position = 1;
bmsf(1,1).score = badness(1);
for i = 2:noligos
  bmsf(1,i).position = nan;
  bmsf(1,i).score = inf;
end;


for x = 2:goodlen
  bmsf = [bmsf ; bmsf(x-1,:)];
  
  for k = 1:noligos
    potentialscore = inf;
    if k == 1
      potentialscore = badness(x);
    else
      if x-probe_spacer_len > 0
        if ~isnan(bmsf(x-probe_spacer_len,k-1).position)
          potentialscore = calcscore(bmsf(x-probe_spacer_len,k-1).score,k-1,badness(x));
        end;
      end;
    end;
    if potentialscore < bmsf(x,k).score
      bmsf(x,k).position = x;
      bmsf(x,k).score = potentialscore;
    end;
  end;
end;

for k = 1:noligos
  matches = [];
  x = goodlen;
  currk = k;
  if ~isnan(bmsf(x,currk).position)
    output(k).score = bmsf(x,currk).score;
    for counter = 1:k
      prevmatch = bmsf(x,currk).position;
      matches = [matches prevmatch];
      x = prevmatch-probe_spacer_len;
      currk = currk - 1;
    end;
    output(k).matches = matches;
  end;
end;

scores = [output.score];

inds = find(scores < 1000000);
output = output(1:max(inds));
for i = 1:length(output)
    output(i).matches = fliplr(output(i).matches);
end;

