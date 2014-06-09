function [inseq2,inseq2_old,merge_action] = readinseq2(inseq,outfile,shortlen,blocklen)
% READINSEQ2 - Check if there is a previously filtered inseq2 stored file 
%   Pre-filtered sequences should be used when available to avoid re-doing
%   the process of checking for mRNA, uRNA, and pseudogene sequences
%
% Input Args:
%   inseq   = single/multi-sequence FASTA as a single string. The latest
%     version of the sequence to design probes against.
%   outfile = a file prefix used for all log/results files in PROBE DESIGN
%   shortlen & blocklen - basic parameters set in PROBE DESIGN
%
% Output Args:
%   inseq2 = the sequence to be used for probe design, may or may not
%     need filtering which is specified by merge_action
%   inseq2_old = in the case of adding more length to the input sequence
%     pre-filtered upstream sequence can still be used. This will be
%     upstream of inseq2 (consisting only of additional non-filtered
%     downstream sequence) after it is filtered and appended to create the
%     full-length, filtered version of inseq.
%   merge_action = a status variable to indicate what further steps should
%     be taken to create a full-length filtered inseq. Can be one of:
%     {'reuse', 'append', 'cut', 'new', 'quit'}
%
%


inseq2filename = [outfile '_inseq2.txt'];
merge_action = '';  % reuse, append, cut, new, quit
inseq2_old   = '';  % set only on 'append' action
try 
    inseq2 = fileread(inseq2filename);
    fprintf('Found %s to use as filtered ''H-M-P'' sequence...\n',inseq2filename);
    
    % Ask for user input about whether or not to use this file
    reply = input('Use this file (y/n)?: ','s');
    if regexp(reply,'[^yY]')  % If reply doesn't contain a 'y' or 'Y'
        error('msg') % prints nothing, goes to catch expression
        % User chose not to use inseq2filename, it will be overwritten
    else
        merge_action = 'reuse';
    end
    
    % Check length of file. If empty, quit. decide whether to append, cut, or 
    seq2len = length(inseq2);
    if seq2len == 0
        fprintf('*** NOTE *** %s is EMPTY...\n',inseq2filename);
        merge_action = 'new';
        reply = input('Continue (y/n)??? ','s');
        if regexp(reply,'[^yY]')  % If reply doesn't contain a 'y' or 'Y'
            merge_action = 'quit';
        end
        error('msg') % prints nothing, goes to catch expression
    
    elseif seq2len < length(inseq) % APPEND
        msg =['-------------------------------------------\n*NOTICE:\n'...
              'Input sequence is LONGER than stored version which ' ...
              'was previously filtered.\nWill filter and APPEND ' ...
              'the additional sequence length to %s...\n'];
        fprintf(1,msg,inseq2filename);
        reply = input('Continue (y/n)??? ','s');
        if regexp(reply,'[^yY]')  % If reply doesn't contain a 'y' or 'Y'
            merge_action = 'quit';
            error('msg')  % prints nothing, goes to catch expression
        else
            merge_action = 'append';
        end
        
        % In order to make sure masking happens properly at the border of
        % the new and old sequence, we cut back the shorter stored
        % sequence and extend one testseq length into the new inseq overlap
        shift = shortlen-blocklen;
        inseq2_old = inseq2; % save filtered version for later append
        inseq2 = inseq(seq2len-shift:end);% use overlap & additional seq for filtering
           
    elseif seq2len > length(inseq)  % CUT
        msg =['Input sequence is SHORTER than stored version which ' ...
              'was previously filtered.\nThis excess sequence will ' ...
              'be DELETED in %s\n'];
        fprintf(1,msg,inseq2filename);
        merge_action = 'cut';
        reply = input('Continue (y/n)??? ','s');
        if regexp(reply,'[nNoO]')  % If reply doesn't contain a 'y' or 'Y'
            merge_action = 'quit';
            error('msg')  % prints nothing, goes to catch expression
        else
            inseq2 = inseq2(1:length(inseq));
        end
    end
     
catch
    if strcmpi(merge_action,'quit')
        error('User quit...'); 
    end
    fprintf('\nPerforming NEW filtering...\n\n');
    merge_action = 'new';
    inseq2 = inseq;  % create a 'working version' of the input sequence
    
end
