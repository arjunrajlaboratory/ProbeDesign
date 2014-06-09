function GC = getGC(seq)

seqprop = oligoprop(seq,'Salt',0.33,'Temp',30);
GC = seqprop.GC;  % Use latest SantaLucia data