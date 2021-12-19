function [FOPT, FWPT, FWIT, TCPU] = GetProductiondata(filename, N)

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

C = textscan(fid, '%f %f %f %f %f %f %f');
%% Select the values

FOPT = C{1,3} * Index_cf(~isempty(find(Index_row(1,:)> column, 1))+1,1);
FWPT = C{1,4} * Index_cf(~isempty(find(and(Index_row(1,:)< column, column<Index_row(2,:)), 1))+1,2);
FWIT = C{1,5} * Index_cf(~isempty(find(and(Index_row(2,:)<=column, column<=Index_row(3,:)), 1))+1,2);
TCPU = C{1,end};
%% File close
fclose(fid);

end
end