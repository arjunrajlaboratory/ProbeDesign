function outstr = mask_string(inseq,mask,chr)

outstr = inseq;
outstr(mask>0) = chr;