%% BLASTHITSBATCH 
% BLAST the input sequence for hits as a batch job

%% Description
% Constructs batch NCBI BLAST job requests for all potential oligo 
% sequences and then retreives BLAST report overview HTML to extract the
% total number of hits for each oligo.  
%
% Keeping track of submitted jobs via the RID (aka BLAST jobID) is
% through a struct() RIDMap. This key/value model chosen over sequential
% job submission and retrieval because out-of-order-sequence handling of
% pre-masked regions or errors with BLAST jobs.
%
% The BLAST report retrieval process is error prone due to failed jobs, 
% so there is a built-in resubmission process of individual oligos and 
% faulty oligos are given a hit count of ZERO.

%% Input
% * |inseq|      - single or multiple sequence FASTA formated single string
% * |batch_size| - the number of sequences submitted per BLAST request
% * *|shortlen|* & *|blocklen|* - basic parameters used to determine oligo length

%% Output
% * counts = 2-column matrix. [start_seq_index & hit_count]  Column 1 has
%     integer values which specify the starting index of a tested sequence.
%     Column 2 is the resulting hit count from the BLAST report.
% * RIDMap = struct()
%     RIDMap.(RID_RIDValue) : Stores a array of integers which are the 
%     starting indicies of each XX-mer oligo submitted in the Nucleotide 
%     BLAST job batch'd with ~50 sequences. We use the 'RID\underscore\11CharRIDValue' 
%     just in case the RID starts with a numeric and so Matlab doesn't spit
%     out an error about invalid field name. This adds some extra code but
%     should prevent errors and future-proof any changes to RID
%     construction by NCBI
%

%% Example Usage:
% |>> [counts,RIDMap] = blasthitsbatch(inseq,'human',50,22,2)|


function [counts, RIDMap] = blasthitsbatch(inseq,database,batch_size,shortlen,blocklen)

% goodlen is the length of the input sequence minus the query seq length
goodlen = length(inseq)-shortlen;

%----------------------------
% Start filtering sequences
%----------------------------

sendseqs = '';
sendseqinds = [];
MAX_WAIT = 0;
fprintf('***\nSearching for hits to %s genome and mRNA via BLAST...\n',database);
for i = 1:goodlen
  testseq = inseq( (i+2):(i+shortlen-blocklen-1-2) );
  
  % for all masked (Psuedo, Repeats, or FASTA header char '>') characters,
  % then skip to next character in the sequence (eg. next in FOR loop). 
  % By skipping them, we avoid the BLAST workload of pre-masked sequences
  if regexp(testseq,'[^gGtTcCaA]')  % if contains chars that aren't G C A or T
      continue
  end
  
  % Construct the batched multi-sequence FASTA to send to BLAST
  sendseqs = [sendseqs sprintf('>seq%.0f\n%s\n',i,testseq)];
  sendseqinds = [sendseqinds i];
  
  if length(sendseqinds) >= batch_size || i==goodlen
      % Submit the BLAST job to NCBI
      [RID, RTOE] = blastncbishortseqs(sendseqs,database);
      if ~isempty(RID)
          % Store the RID and submitted sequence indicies
          RIDMap.(['RID_' RID]) = sendseqinds; 
          if RTOE > MAX_WAIT
            MAX_WAIT = RTOE;
          end
          fprintf(1,'BLASTing %.0f seqs, RID: %s, RTOE %.0f seconds\n',length(sendseqinds),RID,RTOE)
      end
      
      sendseqs = '';
      sendseqinds = [];
  end
  
end

fprintf(1,'Waiting %.0f seconds...\n',MAX_WAIT*3)
pause(MAX_WAIT*3);


%-----------------------------------------------------------------
% Get the BLAST reports for each RID and store the number of HITS
%-----------------------------------------------------------------

counts = zeros(0,2);  % This is to handle the case w/ no BLAST jobs

% Check if there were any BLAST jobs. Its possible that all testseqs had
% masked regions so no BLAST jobs were started. This usually happens after
% appending more sequence to a previous probe design and reruning the job
if exist('RIDMap') 
    fprintf(1,'Now retrieving BLAST reports from each submitted job...\n')
    RIDS = fieldnames(RIDMap);  % TODO: Make sure RIDMap exists, could be empty
else
    RIDMap = [];  % Assign an empty variable so there is something to return
    RIDS = [];  % This causes skipping of the following FOR loop
end

for i = 1:numel(RIDS)
    
    % Store RID value, exclude the prefix 'RID_' created for the field name
    RID = RIDS{i}(5:end);
    
    % A row vector of start indicies for the oligos submitted for this RID
    seqinds = RIDMap.(RIDS{i}); 
    
    hits = [];
    % go get the report and count the number of hits
    hits = getblasthitcount(RID,MAX_WAIT*2);
    
    % handle error in BLAST report retrieval/job, resubmit as individuals
    if isempty(hits)
       
        fprintf(1,'ERROR with RID: %s, resubmitting as individual sequneces...\n',RID)
        single_RIDS = {}; % Used only in this single submission section
        single_seqs = {}; % Used only in this single submission section
        MAX_WAIT = 0;
        for k=1:numel(seqinds)  % Loop over individual sequences, submit one-by-one
        	
            ind = seqinds(k);
            sendseq = inseq( (ind+2):(ind+shortlen-blocklen-1-2) );
            
            [RID,RTOE] = blastncbishortseqs(sendseq,database);
            
            if ~isempty(RID)
                RIDMap.(['RID_' RID]) = sendseq;
                single_RIDS = [single_RIDS RID];
                single_seqs = [single_seqs sendseq];
                fprintf(1,'\tBLASTing %s, RID: %s, RTOE: %.0f seconds\n',sendseq,RID,RTOE)
                if RTOE > MAX_WAIT
                    MAX_WAIT = RTOE;
                end
            end
        end
        
        fprintf(1,'\tWaiting %.0f seconds for jobs to finish...\n',MAX_WAIT*3)
        pause(MAX_WAIT*3);

        for j = 1:numel(single_RIDS)
            
            % go get the report and count the number of hits
            hit = getblasthitcount(single_RIDS{j},MAX_WAIT*2); % hit variable could be an array
            
            % don't handle error, just store a ZERO for the hits
            if isempty(hit)
                hit = 0;
                fprintf(1,'\tERROR: BLAST unsucessful for RID:%s and Seq:%s\n',single_RIDS{j},single_seqs{j})
            end
            
            hits = [hits hit];
        end
        
    end
    
    % Construct 2-column matrix. [start_seq_index & hit_count]
    counts = [counts; seqinds' hits'];
end
  
% TRANSFORM COUNTS INTO A ONE-column vertical vector
% Fill-in gaps in the counts matrix with hit_value=0 for places that were 
% masked locally for mRNA, uRNA, Pseudogenes, and repetitive elements
fixed_counts = [];
for i=1:length(inseq)
    % check if there is already a value at this sequence index
    hit_row = find(counts(:,1)==i); 
    if isempty(hit_row)
        fixed_counts = [fixed_counts; 0];
    else
        fixed_counts = [fixed_counts; counts(hit_row,2)];
    end
end
counts = fixed_counts;
clear fixed_counts;


end
