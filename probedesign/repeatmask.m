function [headers, masked_sequences] = repeatmask(infile, dnasource)
% REPEATMASK  Matlab interface to RepeatMasker.org 
%   Sequences contained in the input FASTA file will be screened against
%   libraries of repetitive elements known for the specified source
%   species. Uses Matlab's URLREAD() function with POST method request.
%   Parameters in the request mirror the HTML form found on:
%         http://www.repeatmasker.org/cgi-bin/WEBRepeatMasker
%
%   Input Args:
%   infile    - the FASTA sequence(s) to be masked
%   dnasource - a string value of species to use as a repeat source
%
%   Output Args:  (these are same as calling [h,s] = FASTAREAD() directly)
%   headers   - the header strings separating multiple FASTA sequences
%   masked_sequences - the masked sequences
%
% Example Usage:
%-------------------------------------------------------------
% % >> [h,seqs] = repeatmask('UBL5.fasta','mouse');
% DNA source to mask for: mouse
% RepeatMasker job is running for UBL5.fasta
% Results will be available at URL: http://www.repeatmasker.org/tmp/cda7a6b0f11075f39fe59b95b6d95cad.html
% Results not available yet, retrying in 5 sec
% *** SUMMARY OF REPEATMASKER JOB ***
% 
% ==================================================
% file name: RM2sequpload_1274302088  
% sequences:             9
% total length:       2230 bp  (2230 bp excl N/X-runs)  
% GC level:         52.33 %
% bases masked:        623 bp ( 27.94 %)
% ==================================================
%                number of      length   percentage
%                elements*    occupied  of sequence
% --------------------------------------------------
% SINEs:                 5          623 bp   27.94 %
%       Alu/B1           4          543 bp   24.35 %
%       B2-B4            0            0 bp    0.00 %
%       IDs              0            0 bp    0.00 %
%       MIRs             1           80 bp    3.59 %
% 
% LINEs:                 0            0 bp    0.00 %
%       LINE1            0            0 bp    0.00 %
%       LINE2            0            0 bp    0.00 %
%       L3/CR1           0            0 bp    0.00 %
% 
% LTR elements:          0            0 bp    0.00 %
%       ERVL             0            0 bp    0.00 %
%       ERVL-MaLRs       0            0 bp    0.00 %
%       ERV_classI       0            0 bp    0.00 %
%       ERV_classII      0            0 bp    0.00 %
% 
% DNA elements:          0             0bp    0.00 %
%       hAT-Charlie      0            0 bp    0.00 %
%       TcMar-Tigger     0            0 bp    0.00 %
% 
% Unclassified:          0            0 bp    0.00 %
% 
% Total interspersed repeats:       623 bp   27.94 %
% 
% 
% Small RNA:             0            0 bp    0.00 %
% 
% Satellites:            0            0 bp    0.00 %
% Simple repeats:        0            0 bp    0.00 %
% Low complexity:        0            0 bp    0.00 %
% ==================================================
% 
% * most repeats fragmented by insertions or deletions
%   have been counted as one element
%                                                       
% 
% The query species was assumed to be mus musculus  
% RepeatMasker version development-$Id: RepeatMasker,v 1.12 2009/07/16 20:25:11 rhub
% run with blastp version 3.0SE-AB [2009-10-30] [linux26-x64-I32LPF64 2009-10-30T17:06:09]
% RepBase Update 20090604, RM database version 20090604
% Retrieving file http://www.repeatmasker.org/tmp/RM2sequpload_1274302088.masked.txt
% >>
%--------------------------------------------------------------------------

% GN changed this so that it works if you put in a sequence string instead
% of a filename too.
if ~isempty(dir(infile(1:min(256,length(infile)))));
    contents = fileread(infile);
else
    contents = infile;
end
    
% acceptable options for DNA source (from the HTML of RepeatMasker form)
oksource = {'vertebrate','mammal','human','rodent','mouse','rat',...
            'artiodactyl','cow','pig','carnivore','cat','dog','chicken',...
            'xenopus','fugu','danio','cionasav','drosophila','anopheles',...
            'elegans','diatom','arabidopsis','wheat','maize','oryza'};
if (ischar(dnasource) && any(strcmpi(dnasource,oksource)))
    disp(['DNA source to mask for: ' dnasource])
else
    error(['DNA source argument: ' dnasource ' is not supported']);
end

% Construct the param/value pairs in an cell array of strings
params = {'sequence',contents,...
          'engine','wublast',...
          'speed','default',...
          'dnasource',dnasource,...
          'ReturnFormat','links',...
          'ReturnMethod','html'};

% URLs to point to the CGI RepeatMasker web service page
baseURL = 'http://www.repeatmasker.org';
actionURL = '/cgi-bin/WEBRepeatMasker';

% POST request to the CGI while providing the input params & sequences
[page,status] = urlread([baseURL actionURL],'post',params);
if status == 0
    error('Could not connect with RepeatMasker service');
end

% Search the HTML output to find the name of the masked FASTA file
% square brackets used for "match one of these characters", this will
% handle the ambiguity of single or double quote usage in HTML attributes
result = regexp(page,'href\s*=.*(.*/tmp/\w+\.html).*','tokens');
if length(result) ~= 1
    error('Found zero or >1 match when searching for the result page URL');
end
resultURL = cell2mat(result{1});
if isempty(strfind(resultURL,'repeatmasker.org'))
    resultURL = [baseURL resultURL];
end

% GN added this line to change depending on whether a sequence or a
% filename is specified.
if ~isempty(dir(infile(1:min(256,length(infile)))));
    disp(['RepeatMasker job is running for ' infile])
else
    disp(['RepeatMasker job is running for entered sequence'])
end



disp(['Results will be available at URL: ' resultURL])
pause(5);

% Get the results page HTML
[page,status] = urlread(resultURL);
if status == 0
    error('Unable to retrieve results from RepeatMasker.org');
else
    summary_text = regexp(page,'</h2><PRE>(.+ version \d+).*</PRE>','tokens');
    % If the job is still running, the URL should have pointed to an 
    % empty endpoint, so we can just pause for 5sec and retry.
    % What can also happen is that the page can show text
    % that is not yet the results (Ie. summary_text is still empty)
    % so in these cases we should wait 5sec and retry
    
    
    while or(isempty(page),isempty(summary_text))
        
        if regexp(page,'No repetitive sequences were detected')
            disp('No repetitive sequences were detected for the input file')
            headers = '';
            masked_sequences = '';
            return
        end
        
        disp('Results not available yet, retrying in 5 sec');
        pause(5);
        [page,status] = urlread(resultURL);
        if status == 0
            error('Unable to retrieve results from RepeatMasker.org');
        end
        summary_text = regexp(page,'</h2><PRE>(.+ version \d+).*</PRE>','tokens');
    end
end

%----------------------------------------------------------------------
% Search the results page URL to find the URL of the masked FASTA file
%----------------------------------------------------------------------

% Let's also print out the summary information
summary_text = regexp(page,'</h2><PRE>(.+ version \d+).*</PRE>','tokens');
disp('*** SUMMARY OF REPEATMASKER JOB ***')
disp(cell2mat(summary_text{1}))

% square brackets used for "match one of these characters", this will
% handle the ambiguity of single or double quote usage in HTML attributes
result = regexp(page,'href\s*=["''](\S+\.masked\.txt)["'']','tokens');
if length(result) ~= 1
    error('Found zero or >1 match when searching for the masked infile');
end

% strcat returns a cell array of str since result is cell array of str so
% need to convert to regular character array
fileURL = cell2mat(strcat(baseURL, result{1})); 
fprintf('Retrieving file %s...\n', fileURL)
[headers, masked_sequences] = fastaread_blah(fileURL);

fprintf(['*** Repeat masking finished ***\n'...
        '-------------------------------\n\n']);

end
