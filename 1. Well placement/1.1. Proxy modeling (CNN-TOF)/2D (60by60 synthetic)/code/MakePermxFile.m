function MakePermxFile(permx, filename)
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', 'PERMX');
fprintf(fid, '%d\n', permx);
fprintf(fid, '/');
fclose(fid);
end