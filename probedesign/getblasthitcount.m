function hits = getblasthitcount(RID,RTOE)
% getblasthitcount  Matlab interface to short (<30nt) nucleotide BLAST 
%   Uses the BLAST URL API documentation to construct a URL that requests a
%   BLAST report from NCBI's servers. 
%        http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/index.html
%   
%  Input Args:
%   RID      - BLAST report ID 
%   RTOE     - Request Time Of Execution, use this to pause() before
%              retrieving the report
%
%  Output Args:  (these are same as calling [h,s] = FASTAREAD() directly)
%     hits - row vector of BLAST hit counts for the input sequence(s)
%
% Example Usage:
%-------------------------------------------------------------
% >> hits = getblasthitcount('0AYH5WBV014',12)
%--------------------------------------------------------------------------

URL = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi';


params = {'CMD' 'Get' 'RID' RID ...
          'DESCRIPTIONS' '0' ...
          'ALIGNMENTS' '0' ...
          'SHOW_OVERVIEW' 'yes' ...
          'FORMAT_TYPE' 'HTML'
         };

status = 0;
retry=0;
while (status ==0 && retry <5)
    [contents,status] = urlread(URL,'post',params);
    if status==0, pause(2); end
    retry = retry+1;
end

if status == 0
    fprintf('Could not connect with BLAST service to receive RID: %s', RID);
    hits = [];
    return
end

status = regexp(contents,'QBlastInfoBegin.*Status=(\w+).*QBlastInfoEnd','tokens');
if isempty(status)
    fprintf('Unable to parse a status string from BLAST report, retrying...')
    status = 'waiting';
else
    status = cell2mat(status{1});  % may not have a string here, check for error options
    fprintf(1,'BLAST report status: %s\n', status)
end

% Sometimes there is an unresolvable error on NCBIs servers, the calling
% script can interpret this as a reason to resubmit with individual
% sequences, or just mark it as an error in the process
unresolve_errormsg = ['An error has occurred on the server. Contact ' ...
    'Blast-help@ncbi.nlm.nih.gov and include your RID: (\w+)'];

retry_wait=0;
while (strcmpi(status,'waiting') && retry_wait < 3)
    fprintf('Will check with BLAST server again in %.0f seconds...\n',RTOE)
    pause(RTOE);

    status = 0;
    retry=0;
    while (status ==0 && retry <5)
        % Try the retrieval again, retrying for URLREAD errors, not BLAST
        % errors or if the job is still running
        [contents,status] = urlread(URL,'post',params);
        if status==0, pause(2); end
        retry = retry+1;
    end
    if status == 0
        fprintf('Could not connect with BLAST service to receive RID: %s\n', RID);
        hits = [];
        return
    end
    
    errorRID = regexp(contents,unresolve_errormsg,'tokens');
    if ~isempty(errorRID)  % we found the error message
        fprintf(1,'There was an NCBI error with RID: %s\n',RID);
        hits = [];
        return
    end
        
    status = regexp(contents,'QBlastInfoBegin.*Status=(\w+).*QBlastInfoEnd','tokens');
    if isempty(status)
        fprintf('Unable to parse a status string from BLAST report, retrying...')
        status = 'waiting';
    else
        status = cell2mat(status{1});  % may not have a string here, check for error options
        fprintf(1,'BLAST report status: %s\n', status)
    end
    retry_wait = retry_wait+1;
    
end
    

hits = [];
if strcmpi(status,'waiting')
    fprintf(1,'Giving up on waiting for RID: %s\n', RID)
    return
end

hit = regexp(contents,'Distribution of (\d+) Blast Hits','tokens');
for i=1:numel(hit)
    hits = [hits str2num(cell2mat(hit{i}))];
end

% % PLOT
% bar([3 4 5]) 
% set(gca,'XTickLabel',{'me','myself','&I'})

end
