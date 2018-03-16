% Runs bowtie on remote server to see if any of the sequences align
% to particular hits.
% inseq is the fasta formatted sequence (not the stripped seq)
% merlength is the length of substring used to search
% database is the database to search (e.g., 'human','humanMito')
% Output is hts, which is number of hits for all these things.

function hts = find_bowtie_hits(inseq,merlength,database)


dom = xmlread('sample_bowtie_search.txt');

%endpoint = 'http://192.168.1.74:8000/service';
%endpoint = 'http://158.130.146.251:8000/service';
%endpoint = 'http://158.130.14.49:8000/service';
endpoint = 'http://158.130.14.239:8000/service';


ch = dom.getElementsByTagName('bow:mer_length');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(num2str(merlength));


ch = dom.getElementsByTagName('bow:seq_database');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(database);


ch = dom.getElementsByTagName('bow:inseq');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(inseq);

resp = callSoapService(endpoint,'screen_sequence',dom);

t = parseSoapResponse(resp);
hts = t.hits.integer;



