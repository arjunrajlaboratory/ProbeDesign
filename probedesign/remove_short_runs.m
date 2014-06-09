function out = remove_short_runs(inmask,n,tolerance)
% From an array of 1s and 0s, removes all runs of 1s <= n in length
% Try n=20 and tolerance = 2.  This removes all runs <18 in length, or all
% runs of 20 in length with 2 0s in the middle.

out = zeros(size(inmask));

for i = 1:length(inmask);
    if inmask(i) == 1;
        temp = inmask(i:min(length(inmask),(i+n-1)));
        if sum(temp) >= n-tolerance
            out(i:min(length(inmask),(i+n-1))) = inmask(i:min(length(inmask),(i+n-1)));
        end;
    end;
end;
        