function [outstring, mismatches] = readfasta(infile)

fid = fopen(infile);

outstring = [];

while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  if ~isempty(tline)
    if tline(1) == '>'
      if ~isempty(outstring)
        outstring = [outstring 'M'];
      end;
    else
      outstring = [outstring tline];
    end;
  end;
end;

mismatches = outstring=='M';
outstring(mismatches) = 'T';

fclose(fid);
    


