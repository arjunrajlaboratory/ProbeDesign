function out = get_acgt(inseq)

if isstr(inseq)
    out = getacgt(inseq);
elseif iscellstr(inseq)
    for i = 1:length(inseq)
        out(i,:) = getacgt(inseq{i});
    end;
end;

function ot = getacgt(inseq)

inseq = lower(inseq);
a = length(find(inseq == 'a'));
c = length(find(inseq == 'c'));
g = length(find(inseq == 'g'));
t = length(find(inseq == 't'));

ot = [a c g t]/length(inseq);