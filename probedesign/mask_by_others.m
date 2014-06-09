% masks oligos in one sequence by a set of other sequences

S = fastaread('HPV16E1.txt');
seqs{1} = S.Sequence(ismember(S.Sequence,'actg'));

S = fastaread('HPV33E1.txt');
seqs{2} = S.Sequence(ismember(S.Sequence,'actg'));

S = fastaread('HPV45E1.txt');
seqs{3} = S.Sequence(ismember(S.Sequence,'actg'));

S = fastaread('HPV59E1.txt');
seqs{4} = S.Sequence(ismember(S.Sequence,'actg'));

%%%%%%%%%%%% HPV16

currseq = seqs{1};
maskseq = [seqs{2} seqs{3} seqs{4}];

for i = 1:(length(currseq)-20)
    fprintf('i = %d\n',i);
    testseq = currseq(i:i+20-1);
    for j = 1:length(maskseq)-20
        dist(i,j) = sum(testseq ~= maskseq(j:j+20-1));
    end;
end;

mi = min(dist,[],2);
emask = mi > 5;

[os,ol] = findprobesHD_mask('HPV16E1.txt',40,...
    'blastmask',false,'pseudogenemask',false,'externalmask',emask);

%%%%%%%%%%%% HPV33

currseq = seqs{2};
maskseq = [seqs{1} seqs{3} seqs{4}];

for i = 1:(length(currseq)-20)
    fprintf('i = %d\n',i);
    testseq = currseq(i:i+20-1);
    for j = 1:length(maskseq)-20
        dist(i,j) = sum(testseq ~= maskseq(j:j+20-1));
    end;
end;

mi = min(dist,[],2);
emask = mi > 5;

[os,ol] = findprobesHD_mask('HPV33E1.txt',40,...
    'blastmask',false,'pseudogenemask',false,'externalmask',emask);


%%%%%%%%%%%% HPV45

currseq = seqs{3};
maskseq = [seqs{1} seqs{2} seqs{4}];

for i = 1:(length(currseq)-20)
    fprintf('i = %d\n',i);
    testseq = currseq(i:i+20-1);
    for j = 1:length(maskseq)-20
        dist(i,j) = sum(testseq ~= maskseq(j:j+20-1));
    end;
end;

mi = min(dist,[],2);
emask = mi > 5;

[os,ol] = findprobesHD_mask('HPV45E1.txt',40,...
    'blastmask',false,'pseudogenemask',false,'externalmask',emask);

%%%%%%%%%%%% HPV59

currseq = seqs{4};
maskseq = [seqs{1} seqs{2} seqs{3}];

for i = 1:(length(currseq)-20)
    fprintf('i = %d\n',i);
    testseq = currseq(i:i+20-1);
    for j = 1:length(maskseq)-20
        dist(i,j) = sum(testseq ~= maskseq(j:j+20-1));
    end;
end;

mi = min(dist,[],2);
emask = mi > 5;

[os,ol] = findprobesHD_mask('HPV59E1.txt',40,...
    'blastmask',false,'pseudogenemask',false,'externalmask',emask);




