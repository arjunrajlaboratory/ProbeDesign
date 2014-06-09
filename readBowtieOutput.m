%% readBowtieOutput
% parses the _default_ output from a |bowtie| alignment

%% Description
% The |bowtie| manual (http://bowtie-bio.sourceforge.net/manual.shtml) as of Aug23 2011
% describes the default output as each line being 8 fields separated by tabs. The fields:
%
% # Name of the read that aligned (in our case an index for sub-sequence position)
% # Reference strand aligned to, *+* for forward strand, *-* for reverse
% # Name of reference sequence where alignment occurs, or numeric ID if no name was provided
% # 0-based offset into the forward strand where leftmost character of the alignment occurs
% # Read sequence (reverse-complemented if orientation is *-*)
% # ASCII-encoded read qualities (reversed if orientation is *-*). The encoded quality values
% are on the Phred scale and the encoding is ASCII-offset by 33 (ASCII char *!*).
% # Contains the number of valid alignments for this read besides the one reported by the
% current line of output itself (i.e. a value of 0 means 1 hit, a value of n means n+1 hits)
% # comma separated list of mismatch descriptors in the format |offset:reference-base>read-base|

%% Input
% * *|filename|* - path to the |bowtie| output text file 

%% Output
% * *|outfile|* - MATLAB |struct| with fields corresponding to columns of |bowtie| output 
% 
% Fields as |outfile.*| (see input description for details):
%
% # readNames
% # refStrands
% # refNames
% # offsets
% # readQual
% # numAlign
% # misMatches

%% Example Usage
%  >> out = readBowtieOutput('/path/to/MKNK2_intron_vs_hg19.out');
%  >> readIndex = out.readNames;  % |bar()| works if reads were unamed and numbered by |bowtie|
%  >> hitCounts = out.numAlign;
%  >> bar(readIndex,hitCounts)  % bar graph to visualze number of hits per read
%  >> xlim( [ min(outC{1}) max(outC{1}) ] )

%% Author
% Marshall J. Levesque 2011

function outfile = readBowtieOutput(filename)

fid = fopen(filename,'r');

% each cell in the cell array |outfileC| corresponds to a column of |bowtie| output
outfileC = textscan(fid,'%d%c%s%d%s%s%d%s','Delimiter','\t');
fclose(fid);

% put the input in a useful/descriptive |struct|
outfile.readNames  = outfileC{1};
outfile.refStrands = outfileC{2};
outfile.refNames   = outfileC{3};
outfile.offsets    = outfileC{4};
outfile.readSeqs   = outfileC{5};
outfile.readQual   = outfileC{6};
outfile.numAlign   = outfileC{7};
outfile.misMatches = outfileC{8};
clear outfileC

