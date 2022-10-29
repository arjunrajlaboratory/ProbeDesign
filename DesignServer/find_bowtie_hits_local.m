% Runs bowtie on local server to see if subsequences of input sequences (inseq)
% align to specified databases. Uses bowtie version (0.12.7) available at  
% https://sourceforge.net/projects/bowtie-bio/files/bowtie/ and 
% on Dropbox at https://www.dropbox.com/sh/6q1bnpl60r24byz/AAAGFIxAGV4PKK6qmk--8WDZa?dl=0.
% Pre-indexed genomes are available on Dropbox. 

% inseq is the fasta formatted sequence (not the stripped seq)
% merlength is the length of substring used to search
% database is the database to search (e.g., 'human','humanMito')
% Output is hts, which is number of hits for all these things.
function hts = find_bowtie_hits_local(inseq,merlength,database)
    
    screen_seq_path = which('screen_sequence.py');
    %Get Python executable 
    pe = pyenv;
    pyexe = pe.Executable;
%     syspath = getenv('PATH');
%     if count(syspath, pyexe) == 0
%         newpath = sprintf('%s:%s', syspath, pyexe);
%         setenv('PATH', newpath); 
%     end
    %Remove '' and "" from inseq since these will mess up the shell command
    inseq = strrep(inseq, '''', '');
    inseq = strrep(inseq, '"', '');
    cmd = sprintf("%s %s '%s' '%d' '%s'", pyexe, screen_seq_path, inseq, merlength, database);
    
    fprintf('Screen sequence request %s\n', datetime) %replacing print statements in bowtie_local.py
    [~,stdout] = system(cmd);
    fprintf('Successful screen\n')
    
    tmpHts = regexp(stdout, '[0-9]+', 'match');
    hts = cellfun(@str2double, tmpHts');
end
