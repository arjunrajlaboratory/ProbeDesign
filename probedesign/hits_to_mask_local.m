function outhts = hits_to_mask_local(inseq,mer_length,db,threshold)

hts = find_bowtie_hits_local(inseq,uint16(mer_length),db); %Note the mer_length needs to be integer for subsequent python scripts. Can change to uint32 for longer subsequences.  
idx = find(hts > threshold);

outhts = zeros(size(hts));
for i = 1:length(idx)
    outhts(idx(i):(idx(i)+mer_length-1)) = 1;
end