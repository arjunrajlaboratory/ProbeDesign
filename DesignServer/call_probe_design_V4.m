
function [oligos,alignment] = call_probe_design_V4(filename,noligos,oligo_length,maskingflag,probeprefix,spacer_length,species)
% species is either human, mouse, rat, drosophila, elegans

fid = fopen(filename);  % Read the file.
seq = fread(fid,'*char');
fclose(fid);
seq = seq';

dom = xmlread('sample_findprobesV4_call.txt');

endpoint = 'http://158.130.146.251:8083/service';


ch = dom.getElementsByTagName('fin:inseq');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(seq);


ch = dom.getElementsByTagName('fin:oligo_length');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(num2str(oligo_length));

ch = dom.getElementsByTagName('fin:noligos');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(num2str(noligos));

ch = dom.getElementsByTagName('fin:maskingflag');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(num2str(maskingflag));


ch = dom.getElementsByTagName('fin:probeprefix');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(probeprefix);


ch = dom.getElementsByTagName('fin:spacer_length');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(num2str(spacer_length));

ch = dom.getElementsByTagName('fin:species');
it = ch.item(0);
ch2 = it.getFirstChild;
ch2.setTextContent(species);

resp = callSoapService(endpoint,'screen_sequence',dom);

t = parseSoapResponse(resp);
oligos = t.oligos.ProbeOligo;
alignment = [t.alignment.raw_sequence;[t.alignment.probe_oligos ' '];[t.alignment.labels ' ']];



