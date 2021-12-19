function [Permx]=GetPermx(filename)
 fid = fopen(filename, 'r');
 if fid == -1
     disp('File open not successful');
 else
     fgetl(fid);
     C = textscan(fid,'%f');
     Permx = C(1:end,1);
     Permx = cell2mat(Permx);
 end
 fclose(fid);
end