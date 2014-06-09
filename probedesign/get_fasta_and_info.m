function [genesymbol, FASTAstring] = get_fasta_and_info(hgsid,NCBI_transcriptID)
% GET_FASTA_AND_INFO - download a FASTA file from UCSC genome browser
%   Currently the formatting of the FASTA, including upper/lowercase
%   intron/exons and splitting of coding and non-coding regions, is hard
%   coded into this script.
%
% ** IMPORTANT **  UCSC Genome browser seems to use a browser cookie that I
% assume expires over some period. This is in the form of the 'hgsid'
% parameter. To obtain one, open a web browser window and perform any kind
% of search on http://genome.ucsc.edu/cgi-bin/hgGateway , then check the
% URL field and look for the 'hgsid' value.
%
% Example Usage:
%-------------------------------------------------------------
% >> [sym, chr, fasta] = get_fasta_and_info('160980111','NM_017572');
% Found info for MKNK2 on Chr 19:  Homo sapiens MAP kinase interacting serine/threonine kinase 2 (MKNK2), transcript variant 1, mRNA.
% See gene tracks on UCSC Genome Browser: http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=160980111&db=hg19&position=chr19%3A2037470-2051243
% Found FASTA for NM_017572
% >>
%

FASTAstring = '';
chromosome  = '';
genesymbol  = '';

% http://genome.ucsc.edu/cgi-bin/hgc
% ?g=refGene
% &i=NM_000146

baseURL = 'http://genome.ucsc.edu';
actionURL = '/cgi-bin/hgc';
params  = {'hgsid', hgsid,...
           'g', 'refGene',...
           'i', NCBI_transcriptID};

% GET request to the CGI while providing the input params
[page,status] = urlread([baseURL actionURL],'get',params);
if status == 0
    error('Could not connect with UCSC Genome Browser service');
end
    

% Parse out the required info to construct the FASTA and return info
% <B>Position:</B> <A HREF="../cgi-bin/hgTracks?hgsid=160980111&db=hg19&position=chr19%3A49468566-49470135">chr19:49468566-49470135</A><BR>

% Capture 1)Gene Symbol 2)Description 3)URL to see gene tracks on UCSC 
%         4)Chromosome 5)seq_from 6)seq_to
pattern = ['<H2>RefSeq Gene (\w+)</H2>.*<B>RefSeq:</B>.*' ...
           '<B>Description:</B>(.+)<BR>\s<B>CDS:'];
info = regexp(page,pattern,'tokens');
if isempty(info)
    fprintf('ERROR: Unable to parse info from UCSC for %s\n',NCBI_transcriptID);
else
    genesymbol = info{1}{1}
    gene_desc  = info{1}{2}
%     tracksURL  = [baseURL info{1}{3}];
%     chromosome = info{1}{4};
%     seq_from   = info{1}{5};
% %     seq_to     = info{1}{6};
%     fprintf('Found info for %s on Chr %s: %s\n',...
%                 genesymbol,chromosome,gene_desc)
%     fprintf('See gene tracks on UCSC Genome Browser: %s\n', tracksURL)
end


% Get the FASTA file
       
getFASTAURL = ['http://genome.ucsc.edu/cgi-bin/hgc'...
               '?hgsid=' hgsid...
               '&g=htcDnaNearGene&i=' NCBI_transcriptID ...
               '&o=refGene&boolshad.hgSeq.promoter=0'...
               '&hgSeq.promoterSize=1000&hgSeq.utrExon5=on&'...
               'boolshad.hgSeq.utrExon5=0&hgSeq.cdsExon=on&'...
               'boolshad.hgSeq.cdsExon=0&hgSeq.utrExon3=on&'...
               'boolshad.hgSeq.utrExon3=0&hgSeq.intron=on&'...
               'boolshad.hgSeq.intron=0&boolshad.hgSeq.downstream=0&'...
               'hgSeq.downstreamSize=1000&hgSeq.granularity=feature&'...
               'hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.splitCDSUTR=0&'...
               'boolshad.hgSeq.splitCDSUTR=0&hgSeq.casing=exon&'...
               'boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&'...
               'submit=submit'] ;      
% GET request to the CGI while providing the input params
% [page,status] = urlread([baseURL actionURL],'get',params);
[page,status] = urlread(getFASTAURL);
if status == 0
    error('Could not connect with UCSC Genome Browser service');
end

FASTAstring = regexp(page,'<PRE>\n(.*)</PRE>','tokens');
if isempty(FASTAstring)
    fprintf('ERROR: Unable to parse FASTA from UCSC for %s\n',NCBI_transcriptID);
else
    FASTAstring = FASTAstring{1}{1};
    fprintf('Found FASTA for %s\n',NCBI_transcriptID)
end
       
% URL FOR THE FASTA FILE,
% http://genome.ucsc.edu/cgi-bin/hgc?
% hgsid=160980111
% &g=htcDnaNearGene
% &i=NM_000146
% &c=chr19
% &l=49468565
% &r=49470135
% &o=refGene
% &boolshad.hgSeq.promoter=0
% &hgSeq.promoterSize=1000
% &hgSeq.utrExon5=on
% &boolshad.hgSeq.utrExon5=0
% &hgSeq.cdsExon=on
% &boolshad.hgSeq.cdsExon=0
% &hgSeq.utrExon3=on&boolshad.hgSeq.utrExon3=0
% &hgSeq.intron=on
% &boolshad.hgSeq.intron=0
% &boolshad.hgSeq.downstream=0
% &hgSeq.downstreamSize=1000
% &hgSeq.granularity=feature
% &hgSeq.padding5=0
% &hgSeq.padding3=0
% &hgSeq.splitCDSUTR=on          %SPLIT the Coding/UTR 
% &boolshad.hgSeq.splitCDSUTR=0
% &hgSeq.casing=exon
% &boolshad.hgSeq.maskRepeats=0
% &hgSeq.repMasking=lower
% &submit=submit


end