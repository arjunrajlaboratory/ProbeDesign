function matchScore = computeOligoMatchScore(matchVector,targetMinNumOligos)
% This is the "objective function" that the algorithm tries to maximize. In
% this case, the goal is to get as many oligos onto as many different
% strains as possible, with a goal of a minimum of targetMinNumOligos
% oligos on each strain. We maximize this by maximizing the product, which
% will be biggest if you get a good balance of oligos between the strains,
% especially if the same oligo targets multiple strains.
% 
% Algorithm: for each strain, take the number of matches and multiply them
% all together (maxing out at targetMinNumOligos).  This will try to spread
% the oligos out and also really try to get a single oligo that hits
% multiple targets.

matchVector = squeeze(matchVector);
mx = targetMinNumOligos * ones(size(matchVector));
matchScore = prod(min(matchVector+1,mx+1));
