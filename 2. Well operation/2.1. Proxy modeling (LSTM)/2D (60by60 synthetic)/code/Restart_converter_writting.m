function Restart_converter_writting(dataname,tstep)
fid=fopen('Restart_converter.log', 'w');
    count=fprintf(fid, 'U \n');
    count=fprintf(fid, '%s \n',dataname);
    count=fprintf(fid, '1 \n');
    count=fprintf(fid, '1 \n');
    count=fprintf(fid, '%d \n',tstep);
    count=fprintf(fid, '%d \n',tstep);  
    if exist([dataname, '.F000',int2str(tstep)],'file')==2
        count=fprintf(fid, 'Y \n');
        count=fprintf(fid, 'N \n');    
    else
        count=fprintf(fid, 'N \n');    
    end
 fclose('all');
 
 
%  U 
% SIMPLE_CASE 
% 1 
% 14 
% 1 
% 1 
% Y 
% N