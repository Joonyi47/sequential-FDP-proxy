function [TCPU] = GetProductiondata2(filename, N)

fid = fopen([filename '.RSM'], 'r');
if fid == -1
    disp('File open not successful')
else
    
%% Read the RSM file
for i = 1:10
    aline{i} = fgetl(fid);
end
Index_row = [30;42;54];
Index_cf  = [1,1,1;1000,1000,1000];     % conversion factor
temp      = split(aline{7},' ');
column    = find(strcmp(temp,'*10**3')==1)';
N         = length(column);
if N > 0
    fgetl(fid);
end

C = textscan(fid, '%f %f %f %f %f %f %f %f');
%% Select the values

TCPU = C{1,end};
%% File close
fclose(fid);

end
end