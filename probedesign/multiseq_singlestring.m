function s = multiseq_singlestring(seqs)
% MULTISEQ_SINGLESTRING Concatenate a FASTA into a single string

% In order to handle the case where there is only one sequence in the FASTA file, in
% which case fastaread_blah returns a <1xseq_length> character array instead of a cell
% array, we need to check for cell array
s = [];
if iscell(seqs)
    % Concatenate the sequences from the multisequence FASTA file to the 's' variable,
    % while starting each sequences with the "greater-than" character
    for i = 1:length(seqs)
      s = [s '>' seqs{i}];
    end;
else
      s = ['>' seqs];
end;