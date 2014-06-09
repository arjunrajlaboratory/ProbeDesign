function bd = bad_acgt(inseqs)

tmp = get_acgt(inseqs);

%bd = tmp(:,2:3) >= 0.40;  % Get rid of seqs with >40% C or G
bd = (tmp(:,2) >= 0.40) | (tmp(:,3) >= 0.40);
bd2 = (tmp(:,2) <= 0.10) | (tmp(:,3) <= 0.10);
bd = bd|bd2;

%bd = tmp(:,2) >= 0.40;  % Get rid of seqs with >40% C

%bd = tmp(:,3) <= 0.10;

%bd = max( [max(tmp(:,2:3) >= 0.40,[],2), max(tmp(:,2:3) <= 0.10,[],2)], [], 2 );