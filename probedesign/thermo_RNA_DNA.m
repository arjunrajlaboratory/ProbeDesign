%% thermo_RNA_DNA
% Calculate RNA/DNA base stack thermodynamic values (Sugimoto et al 1995)

%% Description
% Sugimoto 95 parameters for RNA/DNA Hybridization (Table 3)
% "Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid 
% Duplexes" in Biochemistry 1995
%
% And then SantaLucia PNAS 1998 Tm prediction equation with salt adjustment
% used by MATLAB (they get it from OligoCalc)

%% Input
% * *|inseq|* - RNA sequence of the RNA/DNA hybrid ( 5'->3' uracil->thymidine)

%% Output
% * *|Tm|* - melting temperature at 0.33M salt concentration 5e-5M probe conc
%
% * *|dG|* - Gibbs free energy of hybridization (kcal/mol)
%
% * *|dH|* - enthalpy (kcal/mol)
%
% * *|dS|* - entropy (cal/mol)

function [Tm,dG,dH,dS] = thermo_RNA_DNA(inseq)

    % uracil->thymidine
    % table is rows: acgt, columns: actg. (e.g. tbl(2,3) = 'cg')
    delH = [-7.8, -5.9, -9.1, -8.3;...   % aa ac ag at   ** ENTHALPY **
            -9.0, -9.3,-16.3, -7.0;...   % ca cc cg ct   ** kcal/mol **
            -5.5, -8.0,-12.8, -7.8;...   % ga gc gg gt   *** UNITS ***
            -7.8, -8.6,-10.4,-11.5];     % ta tc tg tt

    delS = [-21.9,-12.3,-23.5,-23.9;...  % aa ac ag at     **** ENTROPY ****
            -26.1,-23.2,-47.1,-19.7;...  % ca cc cg ct   ** cal/(mol*Kelvin) **
            -13.5,-17.1,-31.9,-21.6;...  % ga gc gg gt     ****  UNITS  ****
            -23.2,-22.9,-28.4,-36.4];    % ta tc tg tt  

    % sum enthalpy and entropy of RNA-DNA base stacks
    numseq = nt2int(inseq);
    seqlen = length(inseq);
    ind = sub2ind([4 4],numseq(1:seqlen-1),numseq(2:seqlen));
    dH = sum(delH(ind));  % kcal/mol
    dS = sum(delS(ind));  % cal/(mol*Kelvin)


    initH =  1.9; % kcal/mol
    initS = -3.9; % cal/(mol*Kelvin)
    dH = dH + initH;
    dS = dS + initS; 
    dG = dH*1000 - (37+273.15)*dS; % cal/mol
    dG = dG/1000; % kcal/mol
    
    salt = 0.33; % sodium molarity in 2XSSC used in Raj Lab FISH buffers
    ct = 50e-6;  % what MATLAB uses in oligoprop

    % SantaLucia PNAS 1998 eqn 3 for Tm calculation
    % and salt adjustment from OligoCalc website:
    %    http://www.basic.northwestern.edu/biotools/oligocalc.html
    Tm = ( dH*1000/(dS + (1.9872 * log(ct/4))) ) - 273.15 + 16.6*log10(salt);

