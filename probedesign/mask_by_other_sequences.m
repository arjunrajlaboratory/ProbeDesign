
% This function finds, for each oligo substring in inseq, the maximum
% number of base matches it can find in the sequences stored in the file
% altsequence file.
%
% inseq: input sequence as a char array.  Assumed to be lower case and only
% to contain actgnx>, as defined in findprobesHD.
% altsequencefile: FASTA format file that contains all the sequences to
% mask against.
% oligolength: the length of the oligonucleotides to check.

function maxmatches = mask_by_other_sequences(inseq,altsequencefile, oligolength)

maxmatches = zeros(size(inseq));

seqs_to_avoid = fastaread(altsequencefile);

% Lowercase everything.
for i = 1:length(seqs_to_avoid)
    seqs_to_avoid(i).Sequence = lower(seqs_to_avoid(i).Sequence);
end;

for i = 1:(length(inseq)-oligolength + 1)
    
    % pick the oligo to search for in the rest.
    curroligo = inseq(i:(i+oligolength-1));
    
    maxmatch = 0;
    for j = 1:length(seqs_to_avoid)  % Look through the other seqs
        currseq = seqs_to_avoid(j).Sequence;
        
        for k = 1:length(currseq)-oligolength+1  % iterate through the seq itself
            testoligo = currseq(k:(k+oligolength-1));
            temp_max = sum(testoligo == curroligo);  % Number of matches
            maxmatch = max([temp_max maxmatch]);  % Keep max number of matches
        end;
    end;
    maxmatches(i) = maxmatch;  % Store data
end;

