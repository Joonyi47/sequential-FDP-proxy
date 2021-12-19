function [Active]=GetActive(filename)
 fid = fopen(filename, 'r');
 if fid == -1
     disp('File open not successful');
 else
     fgetl(fid);
     C = textscan(fid,'%f');
     Active = C(1:end,1);
     Active = cell2mat(Active);
 end
 fclose(fid);
end