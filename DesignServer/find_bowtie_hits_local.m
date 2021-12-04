% Runs bowtie on local server to see if subsequences of input sequences (inseq)
% align to specified databases. Uses bowtie version () available at *** and 
% on Dropbox at . Pre-indexed genomes are available on Dropbox. 

% inseq is the fasta formatted sequence (not the stripped seq)
% merlength is the length of substring used to search
% database is the database to search (e.g., 'human','humanMito')
% Output is hts, which is number of hits for all these things.
function hts = find_bowtie_hits_local(inseq,merlength,database)
    
    %Add DesignServer directory path to pythonpath
    pythonScript = which('bowtie_local.py');
    [pythonScriptPath,~,~] = fileparts(pythonScript);
    insert(py.sys.path, int32(0), pythonScriptPath)
    
%     %Import python module bowtie_local
%     bowtie_local = py.importlib.import_module('bowtie_local');
    tmpHts = py.bowtie_local.screen_seqence(inseq, merlength, database);
    tmpHts = cell(tmpHts);
    hts = cellfun(@double, tmpHts'); %transpose and convert to array
end