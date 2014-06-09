function outbadness = mask_to_badness(mask,mer_length)

idx = find(mask > 0);

outbadness = zeros(size(mask));
for i = 1:length(idx)
    outbadness(max(idx(i)-(mer_length-1),1):idx(i)) = inf;
end;
