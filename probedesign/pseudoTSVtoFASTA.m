function pseudoTSVtoFASTA(inFilename,outFilename)
% pseudoTSVtoFASTA
% 
% The latest Pseudogene databases from Pseudogene.org have come in TSV file
% format instead of a FASTA list of sequences. This script converts 
% TSV->FASTA line-by-line and writes it to a new file
%
% Input Args:
%   inFilename  - the filename for the TSV ('.txt') file from pseudogene.org
%   outFilename - the filename for the output FASTA file version of the TSV
%
% Example Usage:
%-----------------------------------------------------------------
% >> pseudoTSVtoFASTA('Mouse56.txt','mouse.fasta')
% >>
% This will read Mouse56.txt and create mouse.fasta

% The columns of the TSV
%-----------------------
% ID        (Column 1)
% Chromosome
% Start Coordinate
% Stop Coordiante
% Strand
% Parent Protein
% Protein Start
% Protein Stop
% Parent Gene
% Fraction
% Num Insertions
% Num Deletions
% Num Shifts
% Num Stops
% E Value
% Identity
% PolyA
% Disablements
% Exons
% Introns
% Class
% Sequence  (Column 22)
% Link


fin  = fopen(inFilename,'r');
fout = fopen(outFilename,'w');

headerline = fgetl(fin);
while 1
    line = fgetl(fin);
    if ~ischar(line), break, end
    
    columns = regexp(line,'\t','split');
    ID  = columns{1};
    seq = columns{22};
    
    % Construct the output FASTA so a sequence line contains max 50 NTs
    format_seq = '';
    for i=1:length(seq)
        if mod(i,50)  % returns zero when at 51st character
            format_seq = [format_seq seq(i)];
        else
            format_seq = [format_seq seq(i) '\n'];
        end
    end

    fprintf(fout,['>' ID '\n' format_seq '\n']);

end

fclose(fin);
fclose(fout);
