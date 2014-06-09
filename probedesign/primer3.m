function [left,right,page] = primer3(infile, varargin)
% PRIMER3  Matlab interface to primer3 web service.org 
%   Uses Matlab's URLREAD() function with POST method request.
%   Parameters in the request mirror the HTML form found on:
%         http://frodo.wi.mit.edu/primer3/
%
%   Input Args:
%   infile    - the FASTA sequences filename to design primers for (str)
%   minlength - minimum length of the desired amplicon (int)
%   maxlength - maximum length of the desired amplicon (int)
%   length_to_use - length of template sequence (starting index 1) to use
%   primer_left - specific primer sequence to use as left primer
%   primer_right - specific primer sequence to use as right primer
%   primer_internal - specific primer sequence to use as internal primer
%   misprime_lib - indicates what mispriming library (if any) Primer3 
%                  should use to screen for interspersed repeats or for 
%                  other sequence to avoid as a location for primers
%     'HUMAN','RODENT_AND_SIMPLE','RODENT','DROSOPHILA','NONE'
%
%   Output Args:  
%   left - forward (aka sense, left, 5' ) primer (struct)
%   right - reverse (aka anti-sense, right, 3' ) primer (struct)
%   eg: right = 
%           start: 89
%             len: 19
%              tm: 59.79
%              gc: 57.89
%             seq: 'gaagcctcaacccgttctc'
%
% Example Usage:
%-------------------------------------------------------------
% % >> [left,right] = primer3('COX6B1_int.fasta');
%     or
% % >> [left,right] = primer3('COX6B1_int.fasta',
%                             'minlength',80,
%                             'maxlength',150,
%                             'misprime_lib','HUMAN',
%                             'length_to_use',300,
%                             'primer_left','tgtggagggttctgtaacctg');
% 

%-------------------------------------------------------------------
% Parse the inputs 
%-------------------------------------------------------------------

p = inputParser;

% our input sequence should be in the form of a FASTA file, only the first
% sequence in a multi-sequence FASTA will be used for primer design
p.addRequired('infile',@ischar);
% Use these input variables to set the range of acceptable lengths
p.addOptional('minlength',180,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
p.addOptional('maxlength',220,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
% length of our template sequence (starting from base position 1) to use
% for primer design
p.addOptional('length_to_use',250,@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
% specify left, right, hybridization primer sequence for Primer3 to use.
% Make sure the user is providing a sequence with valid nucleotide chars
p.addParamValue('primer_left','',@(x)~any(regexp(x,'[^AaCcGgTtNn]')));
p.addParamValue('primer_right','',@(x)~any(regexp(x,'[^AaCcGgTtNn]')));
p.addParamValue('primer_internal','',@(x)~any(regexp(x,'[^AaCcGgTtNn]')));
% acceptable options for mispriming library (from the HTML of Primer3 form)
mispriming_libraries = {'HUMAN','RODENT_AND_SIMPLE','RODENT','DROSOPHILA','NONE'};
p.addParamValue('misprime_lib','NONE',@(x)any(strcmpi(x,mispriming_libraries)));

p.parse(infile,varargin{:});


%%%%  DONE PARSING INPUT %%%%%
%-------------------------------------------------------------------

%-------------------------------------------------------------------
% Read the FASTA input file containing the template strand sequences
%-------------------------------------------------------------------

[headers, seqs] = fastaread_blah(infile);
if iscell(seqs) 
    sequence = seqs{1};  % Use only the first sequence in a multi-seq FASTA file
else
    sequence = seqs;  % a one sequence FASTA returns a char array, not cell
end
trim_length = p.Results.length_to_use;
if length(sequence) < trim_length 
    disp('Provided sequence is shorter than specified length to use as template');
    return;
elseif p.Results.maxlength > trim_length
    disp('maxlength of product is longer than length_to_use of input sequence');
    return;
else
    disp(['Trimming input sequence to ',int2str(trim_length),'bp long']);
    sequence = sequence(1:trim_length);
end

% Construct the param/value pairs in an cell array of strings
% Checkbox <input> elements consider anything but an empty string as TRUE
params = {'SEQUENCE',sequence,...
          'PRIMER_MISPRIMING_LIBRARY', p.Results.misprime_lib,...
          'MUST_XLATE_PICK_LEFT','true',...
          'PRIMER_LEFT_INPUT', p.Results.primer_left,...
          'MUST_XLATE_PICK_HYB_PROBE','',...
          'PRIMER_INTERNAL_OLIGO_INPUT', p.Results.primer_internal,...
          'MUST_XLATE_PICK_RIGHT','true',...
          'PRIMER_RIGHT_INPUT', p.Results.primer_right,...
          'PRIMER_PRODUCT_SIZE_RANGE', [int2str(p.Results.minlength),'-', int2str(p.Results.maxlength)],...
          'PRIMER_NUM_RETURN','1',...
          'PRIMER_MAX_END_STABILITY','9.0',...
          'PRIMER_MAX_MISPRIMING','12.00',...
          'PRIMER_PAIR_MAX_MISPRIMING','24.00',...
          'PRIMER_MAX_TEMPLATE_MISPRIMING','12.00',...
          'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING','24.00',...
            ...% General Primer Picking Conditions %
          'PRIMER_MIN_SIZE','18',...
          'PRIMER_OPT_SIZE','20',...
          'PRIMER_MAX_SIZE','27',...
          'PRIMER_MIN_TM','57.0',...
          'PRIMER_OPT_TM','60.0',...
          'PRIMER_MAX_TM','63.0',...
          'PRIMER_MAX_DIFF_TM','100.0',...
          'PRIMER_TM_SANTALUCIA','0',... % Thermo params, see Primer3 page for info
          'PRIMER_MIN_GC','20.0',...
          'PRIMER_MAX_GC','80.0',...
          'PRIMER_SELF_ANY','8.00',...
          'PRIMER_SELF_END','3.00',...
          'PRIMER_NUM_NS_ACCEPTED','0',...
          'PRIMER_MAX_POLY_X','5',...
          'PRIMER_INSIDE_PENALTY','',...
          'PRIMER_OUTSIDE_PENALTY','0',...
          'PRIMER_FIRST_BASE_INDEX','1',...
          'PRIMER_GC_CLAMP','0',...
          'PRIMER_SALT_CONC','50.0',...
          'PRIMER_SALT_CORRECTIONS','0',... % Schildkraut and Lifson 1965
          'PRIMER_DIVALENT_CONC','0.0',...
          'PRIMER_DNTP_CONC','0.0',...
          'PRIMER_LIBERAL_BASE','1',...
          'MUST_XLATE_PRINT_INPUT','1',...
          'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS','0',...
          'PRIMER_LOWERCASE_MASKING','0',...
            ...% Other Per-Sequence Inputs
          'INCLUDED_REGION', '',...
          'PRIMER_START_CODON_POSITION','',...
          'PRIMER_SEQUENCE_QUALITY','',...
          'PRIMER_MIN_QUALITY','0',...
          'PRIMER_MIN_END_QUALITY','0',...
          'PRIMER_QUALITY_RANGE_MIN','0',...
          'PRIMER_QUALITY_RANGE_MAX','100',...
            ...% Objective Function Penalty Weights for Primers %
          'PRIMER_WT_TM_LT','1.0',...
          'PRIMER_WT_TM_GT','1.0',...
          'PRIMER_WT_SIZE_LT','1.0',...
          'PRIMER_WT_SIZE_GT','1.0',...
          'PRIMER_WT_GC_PERCENT_LT','0.0',...
          'PRIMER_WT_GC_PERCENT_GT','0.0',...
          'Pick Primers','Pick Primers'    };

% URLs to point to the CGI RepeatMasker web service page
baseURL = 'http://frodo.wi.mit.edu/';
actionURL = 'cgi-bin/primer3-web-cgi-bin-0.4.0/primer3_results.cgi';

disp(['Submitting the input sequence from: ', infile]);
disp(['Using mispriming library: ',p.Results.misprime_lib])

% POST request to the CGI while providing the input params & sequence
[page,status] = urlread([baseURL actionURL],'post',params);
if status == 0
    error('Could not connect with Primer3 service');
else
    % If the job is still running, the URL should have pointed to an 
    % empty endpoint, so we can just pause for 5sec and retry
    while isempty(page)
        disp('Results not available yet, retrying in 5 sec');
        pause(5);
        [page,status] = urlread(resultURL);
        if status == 0
            error('Unable to retrieve results from Primer3');
        end
    end
end

%-----------------------------------------------------------------------------
% The output from Primer3 that we are interested in:
%-----------------------------------------------------------------------------
% OLIGO            start  len      tm     gc%   any    3'   rep seq 
% LEFT PRIMER        648   20   59.99   55.00  4.00  0.00 10.00 acatgggcttagggtgtgag
% RIGHT PRIMER       888   20   60.01   50.00  4.00  0.00 10.00 ccccaaattatccccactct
%-----------------------------------------------------------------------------

% these regular expression statements return structure arrays with primer
% sequence and properties:
left = regexp(page,'LEFT PRIMER\s+(?<start>\d+)\s+(?<len>\d+)\s+(?<tm>\S+)\s+(?<gc>\S+)\s+\S+\s+\S+\s+\S+\s+(?<seq>\w+)','names');
right = regexp(page,'RIGHT PRIMER\s+(?<start>\d+)\s+(?<len>\d+)\s+(?<tm>\S+)\s+(?<gc>\S+)\s+\S+\s+\S+\s+\S+\s+(?<seq>\w+)','names');

if isempty(left) || isempty(right)
    fprintf(' #! ERROR: Unable to find any primers !#\n');
else
    fprintf(['-----------------------------\n'...
             '******* Primer3 found *******\n'...
             'Left primer:\t%s\n'...
             'Right primer:\t%s\n'...
             '--------------END------------\n']...
             ,left.seq,right.seq);
    % convert strings to useful numbers
    left.start = str2num(left.start);
    left.len   = str2num(left.len); 
    left.gc    = str2num(left.gc);
    left.tm    = str2num(left.tm);
    right.start = str2num(right.start);
    right.len   = str2num(right.len); 
    right.gc    = str2num(right.gc);
    right.tm    = str2num(right.tm);
end



end
