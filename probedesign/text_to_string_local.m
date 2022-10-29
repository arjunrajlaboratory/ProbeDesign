function outstr = text_to_string_local(infile)

% BE changed this so that it works if you put in a sequence string instead
% of a filename too.
if ~isempty(dir(infile(1:min(256,length(infile)))))
    outstr = fileread(infile);
else
    outstr = infile;
end

