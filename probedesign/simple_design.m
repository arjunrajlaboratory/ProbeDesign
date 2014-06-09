function simple_design(seq,prefix)

fprintf(1,'Oligo#\tGC\tGibbs RNA/DNA\tSeq\tOligoName\n');

optimum = -24;  % RNA/DNA hybridization Gibbs free energy (kcal/mol)
minlength = 10;
spacer_size = 2;
spacer = 0;
total_oligos = 0;
oligo = '';
for p = 1:length(seq)
    
    spacer = spacer + 1;
    
    if spacer <= spacer_size
        continue;
    end
    
    oligo = [oligo seq(p)];
  
    if length(oligo) < minlength
        continue
    end

    if any(strfind(lower(oligo),'n'))
        oligo = oligo(2:end);
        continue;
    end

    [Tm,dG,dH,dS] = thermo_RNA_DNA(oligo);
    ediff = abs(dG - optimum);
    
    if ediff < 2 || dG < -26

        if length(oligo) < minlength+12
            oligo = oligo(2:end);
            continue;
        else
            oligo_name = [prefix '_' num2str(total_oligos+1)];
            fprintf(1,'%d\t%.1f\t%.1f\t%s\t%s\n',...
                       total_oligos+1, getGC(oligo), dG, revcomp(oligo), oligo_name);
            oligo = [];  % clear the oligo
            spacer = 0;
            total_oligos = total_oligos + 1;
        end
    end
    
end

fprintf(1,'Found %d oligos\n',total_oligos);
