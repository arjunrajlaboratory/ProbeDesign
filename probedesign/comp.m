function outseq = comp(inseq)

%Everything in lower case:

outseq = lower(inseq);

outseq = strrep(outseq,'a','m');
outseq = strrep(outseq,'t','a');
outseq = strrep(outseq,'m','t');

outseq = strrep(outseq,'g','m');
outseq = strrep(outseq,'c','g');
outseq = strrep(outseq,'m','c');
