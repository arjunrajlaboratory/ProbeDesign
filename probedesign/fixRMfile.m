function outfile = fixRMfile(infile)

outfile = ['fixed_' infile];

fid = fopen(infile,'r');
fout = fopen(outfile,'w');

while 1
  tline = fgetl(fid);
  if ~isempty(tline)
    fprintf(fout,'%s\n',tline);
  end;
    
  if ~ischar(tline), break, end
end;

fclose(fid);
fclose(fout);
