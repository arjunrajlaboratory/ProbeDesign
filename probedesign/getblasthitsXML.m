function hits = getblasthitsXML(RID)
% getblasthits  Matlab interface to retrieve XML BLAST reports and stores hit E-Values
%               ** RECOMMEDED FOR SHORT SETS OF RESULTS ONLY, XML VERY VERBOSE!!! **
%   Uses the BLAST URL API documentation to construct a URL that requests a
%   BLAST report from NCBI's servers. 
%        http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/index.html
%   BLAST report is requested as XML and parsed using Matlab's DOM parsing methods
%   that are buried in it's JAVA-based underbelly.
%
%   TODO: Extend this to parse out even more details from BLAST reports about the hits
%         we find. Interesting info may include: accession numbers, start/end sequence
%         positions, definitions, and the alignment text. The XML parsing is robust and 
%         should be reliable as long as NCBI BLAST doesn't change report formatting
%   
%  Input Args:
%   RID      - BLAST report ID 
%
%  Output Args: 
%     hits - a struct that contains fields BLAST report E-Value results for each 
%            oligo of the query input sequence for this BLAST job
%
% Example Usage:
%-------------------------------------------------------------
% >> RPL18A = getblasthits('0AYH5WBV014',12)
% >> numHitsOligo2 = length(RPL18A.seq2);  
% >> eValuesOligo2 = RPL18A.seq12;  % an array of doubles
%--------------------------------------------------------------------------

URL = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi';

params = ['?CMD='           'GET'...
          '&RID='           RID ...
          '&FORMAT_TYPE='   'XML'
         ];

URL = [URL params];
fprintf(1,'Retreiving the XML formatted BLAST report using request URL:\n%s\n',URL);

DOM = '';
retry=0;
while (isempty(DOM) && retry <5)
    % The URLREAD() method inside of XMLREAD() doesn't seem to check for any failed
    % HTTP responses and will throw and empty result into the XML parser method, thus
    % causing a fatal error and a quitting of the script. We aim to prevent it here
    try   
        DOM = xmlread(URL);
    catch  % no alternative command needed, just check if we have a result
        fprintf(1,'xmlread failed, retrying with BLAST server\n');
    end
    if isempty(DOM), pause(2); end
    retry = retry+1;
end

if isempty(DOM)
    fprintf(1,'Could not connect with BLAST service to receive RID: %s', RID);
    hits = [];
    return
end

% Get the root node out of the full XML DOM (Document Object Model)
BlastOutput = DOM.getDocumentElement;

% Get the gene symbol from the prefix of the oligo name (eg 'EEF2_int_1')
geneSym = BlastOutput.getElementsByTagName('BlastOutput_query-def').item(0).getTextContent;
geneSym = char(geneSym);  % convert from java.lang.String to char array
geneSym = strrep(geneSym,'_1','');  % remove the suffix and keep the Gene Symbol

% Get the set of <Iteration> nodes that contain the info & results for each oligo
Iterations = BlastOutput.getElementsByTagName('Iteration');
numOligos = Iterations.getLength;
fprintf(1,'This BLAST job had %d query sequences\n\n',numOligos);
fprintf(1,'OligoID\tHitCount\tMedianEVal\n');
fprintf(1,'------------------------------\n');
for i = 0:numOligos-1  % using Java array indexing that starts with ZERO
    % this first step returns a java.lang.String so we need to convert
    oligoID = Iterations.item(i).getElementsByTagName('Iteration_query-def').item(0).getTextContent;
    oligoID = char(oligoID);
    hits = Iterations.item(i).getElementsByTagName('Hit');
    numHits = hits.getLength;
    eVals = zeros(0,numHits);
    for j = 0:numHits-1  % again using Java array indexing
        eVal = hits.item(j).getElementsByTagName('Hsp_evalue').item(0).getTextContent;
        eVals(j+1) = str2num(char(eVal));
    end
    fprintf('%s\t%d\t%.4f\n',oligoID,numHits,median(eVals));
    
    % We store the results in a struct named after the gene symbol and a field for each oligo
    % tested in BLAST. Each of these '.int#' fields contain an array with eValues for the hits
    % and the length of this array provides us with the number of hits
    eval([geneSym '.seq' num2str(i+1) ' = eVals;']);  % make struct (eg EEF2.EEF2_int_4 = [0.0145 0.2])
end

hits = eval(geneSym);  % return the struct named after the gene symbol

end  % END OF THE FUNCTION
