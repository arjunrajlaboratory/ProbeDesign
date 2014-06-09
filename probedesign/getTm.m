function Tm = getTm(seq)

seqprop = oligoprop(seq,'Salt',0.33,'Temp',30);
Tm = seqprop.Tm(5);  % Use latest SantaLucia data