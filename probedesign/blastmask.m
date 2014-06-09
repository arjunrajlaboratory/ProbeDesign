function [counts, seqB, RIDMap] = blastmask(inseq,outfile,database,merge_action,shortlen,blocklen)
% BLASTMASK
%
% Input Args:
%   outfile      - file prefix for the data file to read BLAST hit counts
%     and to write data to. Set to EMPTY string for no disk data read/write
%   merge_action - sent by calling script that used prior filtering, such
%     as in findprobes_introns_mjl.m  This function still runs if
%     merge_action is set to an empty string

%----------------------------
% Start filtering sequences
%----------------------------

counts = [];
countsfilename = [outfile '_counts.csv'];

if ~isempty(outfile)
    try 
        counts = csvread(countsfilename);
    catch
        fprintf(1,'Did not find previous BLAST hit count file: %s\n',countsfilename);
    end
end

% Check if we have the BLAST data. Start the BLAST jobs if there is no
% previously stored data. OR if instructed by a calling script via
% merge_action being set to 'append'
if isempty(counts)
    
    % Catch error of asking to append when there is no BLAST data
    if strmatch(merge_action,{'append','cut'})
        msg = sprintf(['ERROR: Cannot append/cut BLAST hit counts '...
                       'without data file %s\n'],countsfilename);
        error(msg);  % THIS SHOULD NOT HAPPEN, so hard quit the script
    end
    
    merge_action = 'new';
    
end
    
if strmatch(merge_action,{'new','append'})
    
    if strcmpi(merge_action,'append')
        fprintf(1,'Using previous BLAST hit count file: %s\n',countsfilename);
        old_counts = counts;
        inseqBLAST = inseq(length(counts)+1:end);
    else
        inseqBLAST = inseq;
    end
    
    %----------------------------------------------
    % Perform the BLAST batch jobs on NCBI servers
    %----------------------------------------------
    batch_size = 50;
    [counts, RIDMap] = blasthitsbatch(inseqBLAST,database,batch_size,shortlen,blocklen);
   
    if strcmpi(merge_action,'append')  % Append the new set of hit counts
        counts = [old_counts; counts];
    end
else
    
    % This is the merge_action case of 'reuse' or 'cut'
    fprintf(1,'Using previous BLAST hit count file: %s\n',countsfilename);
    if strcmpi(merge_action,'cut')  % Shorten the set of hit counts
        counts = counts(1:length(inseq));
    end
    
    RIDMap = 0;
end
    
% Store the results of BLASTHITCOUNTS() for later use and/or review
if (~isempty(counts) && ~isempty(RIDMap) && ~isempty(outfile))
    try
        csvwrite(countsfilename,counts);
    catch
        fprintf('ERROR: could not write to %s, continuing...',countsfilename);
    end
end


% seqB = the masked sequence with 'B' in place of nucleotides making up
%        a sequence that has a hit count exceeding THRESHOLD value
THRESHOLD = 10000; 
seqB = lower(inseq);  % create a 'working version' of the input sequence

% column vector of start indicies for sequences w/ hits > THRESHOLD
% *Note: indicies are 1-to-1 with inseq nucleotide indicies
hit_seq_inds = find(counts>THRESHOLD); 

for i = hit_seq_inds'  % iterate over rows w/ hit counts above THRESHOLD
  hit_value = counts(i);
  displayseq = inseq( (i+2):(i+shortlen-blocklen-1-2) );
  fprintf('Seq: %s had %.0f hits!\n',displayseq,hit_value);
  seqB( (i+2):(i+shortlen-blocklen-1-2) ) = 'B';
end