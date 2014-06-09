function outstr = text_to_string(infile)

fid = fopen(infile,'r');
outstr = fread(fid,'uint8=>char')';
fclose(fid);