%% make_oligo_order
% Create a TSV file for importing into Biosearch 96-well plate format in M$OFT Excel

%% Description
% Lists of oligos in text files are the output of probe design (|findprobesHD|).
% These lists of oligos have a common filename format (e.g. GENESYMB_int_oligos.txt)
% so this script takes a list of GENESYMB, looks for the files in the current 
% working directory, reads their contents, and creates the order TSV.
%
% The oligo order format is:
%
%  ROW   COLUMN  Sequence    Sequence Name
%  (A-H) (1-12)  (5' to 3')  GENESYM_oligo#
%
% The probe oligo file format from probe design is:
%
%  oligo#   GC% TM  Sequence    OligoName
%   1       50  67.6    cagagtcaagagactcagga    RPL28_int_1
%
%
% This script is especially useful since it will pad probes with empty wells to
% fill a 96-well plate for easy pipetting.  For example for probes W, X, Y, Z
%     1  2  3  4  5  6  7  8  9 10 11 12
%    -----------------------------------
% A | W  W  X  X  Y  Y  Z  Z  Z         |
% B | W  W  X  X  Y  Y  Z  Z  Z         |
% C | W  W  X  X  Y  Y  Z  Z  Z         |
% D | W  W  X  X  Y  Y  Z  Z  Z         |
% E | W  W  X  X  Y  Y  Z  Z  Z         |
% F | W  W  X  X  Y  Y  Z  Z            |
% G | W  W  X     Y  Y  Z  Z            |
% H | W  W  X     Y  Y  Z  Z            |
%    -----------------------------------
 

%% Input
% * *|geneSymbols|* - cell array of strings with oligo file name prefixes 

%% Output
% * *|plates|* - cell array of cell arrays [96 x n] for n 96-well plates

%% Example Usage
%  >> ordergenes = {'MIDN','CIRBP','RPL36','P2RY11','SWSAP1',...
%                   'IER2','C19orf42','LOC148189','PAF1','SERTAD1',...
%                   'ZNF576','PVR','EMP3','PPP1R15A','RPL13A','IRF3',...
%                   'CLEC11A','LINC00085','IL11','UBE2S'};
%  >> plates = make_oligo_order(ordergenes);
%  >> fileread('order.tsv')

%% Author
% Marshall J. Levesque 2012

function plates = make_oligo_order(geneSymbols)

    seqs = {};  
    seqN = {};

    for k = 1:numel(geneSymbols)
        symb = geneSymbols{k}; 
        fname = [symb '_oligos.txt']; 
        A = exist(fname);
        if A ~= 2
            fprintf(1,'ERROR: could not find oligo file %s\n',fname);
        else
            fid = fopen(fname,'r');
            C = textscan(fid,'%d\t%f\t%f\t%s\t%s');
            fclose(fid);
            seqs(k) = C(4);
            seqNames(k) = C(5);
        end
    end

    rows = 'ABCDEFGH'; 
    cols = 1:12;          
    nRows = numel(rows); 
    nCols = numel(cols);     
    wells = createWells(rows,cols);
    plates = {};
    w = 1;
    for k = 1:numel(seqs)  % for each gene symbol with successful file read
        
        seqList = seqs{k};
        seqNms = seqNames{k};
        nSeqs = numel(seqList);
        colsNeeded = ceil(nSeqs/nRows);
        colsCompleted = floor(w/nRows);
        
        if colsNeeded > nCols - colsCompleted
            % fill out wells at end of plate with empties
            while w <= numel(wells)
                wells(w) = {[wells{w} [] []]};
                w = w + 1;
            end

            plates = [plates wells];  % add the plate
            wells = createWells(rows,cols);  % clear the wells
            w = 1;
        end

        for n = 1:nSeqs
            wells(w) = {[wells{w} seqList{n} seqNms{n}]};
            w = w + 1;
        end

        empties = colsNeeded*nRows - nSeqs;
        for n = 1:empties % fill out the rest of the last column
            wells(w) = {[wells{w} [] []]};
            w = w + 1;
        end
    
        if k == numel(seqs)  % finished the last input sequence
            plates = [plates wells];  % add the plate
        end
    
    end


    % Print the order to a TSV file
    fid = fopen('order.tsv','w');
    fprintf(fid,'Row\tColumn\tSequence\tOligoName\n');
    for p = 1:size(plates,2)  % for each plate
        for o = 1:size(plates,1)  % for each oligo / row
            r = plates{o,p};
            if numel(r) == 4
                fprintf(fid,'%s\t%d\t%s\t%s\n',r{1},r{2},r{3},r{4});
            else
                fprintf(fid,'%s\t%d\t\t\n',r{1},r{2});
            end
        end
    end
    fprintf(1,'* created * order.tsv *\n');
    fclose(fid);

            
    function wells = createWells(rows,cols)
        wells = cell(numel(rows)*numel(cols),1);
        k = 1;
        for c = cols
            for r = rows
                wells(k) = {{r,c}};
                k = k + 1;
            end
        end
