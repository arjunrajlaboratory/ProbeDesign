function [RID,RTOE] = blastncbishortseqs(inseqs,database)
% blastncbishortseqs  Matlab interface to short (<30nt) nucleotide BLAST 
%   Uses the BLAST URL API documentation to construct a URL that requests a
%   BLAST job be submitted to NCBI's servers. 
%        http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/index.html
%   
%   Input Args:
%   inseqs    - SHOULD BE A VALID FASTA FILE w/ headers, a Matlab STRING
%   database  - string value, either 'human' or 'mouse'
%
%   Output Args:  (these are same as calling [h,s] = FASTAREAD() directly)
%   RID   - BLAST report ID 
%   RTOE  - Request Time Of Execution, use this to pause() before
%           retrieving the report
%
% Example Usage:
%-------------------------------------------------------------
% >> [RID,RTOE] = blastncbishortseqs('atgcattgtatctcgtagctag','human')
%--------------------------------------------------------------------------

% BLAST URL API to submit a BLAST job
baseURL = 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';

% The special DB value found in the <form> element on the BLAST site that 
% will screen against the human genome and RNA transcript libraries
genometranscriptDB.human = ...
  'dbindex/9606/ref_contig dbindex/9606/alt_contig_HuRef dbindex/9606/rna';
genometranscriptDB.mouse = ...
  'dbindex/10090/alt_contig dbindex/10090/ref_contig dbindex/10090/rna';
genometranscriptDB.drosophila = ...
  'gp/7227.9554/dm_refc gp/7227.9554/dm_refm';  % drosophila genome refseq contig and mRNA
genometranscriptDB.elegans = ...
  'genomes/c_elegans genomes/c_elegans_chr';  % elegans mRNA and chromsomes

if nargin == 1
    database = 'human';
end

% These parameters are the same ones used by NCBI BLAST for short
% nucleotide sequences (<30nt) that are automatically set by NCBI web
% interface. There is a URL param for this autosetting, but it didn't work
% in my experience:
%                      &SHORT_QUERY_ADJUST=true
%      http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node66.html
%
params = ['?CMD=Put' ...
          '&QUERY=' urlencode(inseqs) ...
          '&PROGRAM=blastn' ...
          '&DATABASE=' urlencode(genometranscriptDB.(database)),...
          '&EXPECT=1000' ...
          '&WORD_SIZE=7' ...
          '&FILTER=F' ...
          '&GAPCOSTS=5%202' ... % [5 2]
          '&MATCH_SCORES=1%2C-3' ... % [1 -3]
         ];
      
% request to the CGI while providing the input params
status = 0;
retry=0;
while (status ==0 && retry <5)
    [page,status] = urlread([baseURL params]);
    if status==0, pause(2); end
    retry = retry+1;
end
if status == 0
    error('Could not connect with BLAST service');
end

% Search the HTML output to find the RID and expected time to execute
info = regexp(page,'QBlastInfoBegin.*RID = (\w+).*RTOE = (\d+)','tokens');
if numel(info{1}) ~= 2
    error('Unable to parse out the RID and RTOE from the BLAST response');
end

% regexp() using 'tokens' parameter returns a crazy cell of a cell array
RID = info{1}{1};
RTOE = str2num(info{1}{2});

end
