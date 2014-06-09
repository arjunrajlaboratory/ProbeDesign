function dG = getGibbs_RNA_DNA(seq)

[Tm,dG,dH,dS] = thermo_RNA_DNA(seq);
