function outhts = hits_to_mask(inseq,mer_length,db,threshold)

hts = find_bowtie_hits(inseq,mer_length,db);
idx = find(hts > threshold);

outhts = zeros(size(hts));
for i = 1:length(idx)
    outhts(idx(i):(idx(i)+mer_length-1)) = 1;
end;
