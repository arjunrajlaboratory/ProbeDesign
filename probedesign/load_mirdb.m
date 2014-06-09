
function seqs = load_mirdb(thefile,species)

% E.g., seqs = load_mirdb('mature.fa','hsa');


%species = 'hsa';  % Could be mmu for mouse, cel for elegans, etc.

%s = fastaread('mature.fa');
s = fastaread(thefile);

allhead = {s.Header};
seqs = {s.Sequence};

tf = strncmp(species,allhead,3);

seqs = {seqs{tf}}; % select all sequences matching the species

for i = 1:length(seqs)
  seq = seqs{i};
  seq = lower(seq);
  seq(find(seq == 'u')) = 't';
  seqs{i} = seq;
end;


