function Restart_converter_writting_batch(dataname,tstep,Nidx)


fid=fopen('Restart_converter.log', 'w');
    count=fprintf(fid, 'U \n');
    count=fprintf(fid, '%s \n',[dataname '_' int2str(Nidx(1))]);
    count=fprintf(fid, '1 \n');
    count=fprintf(fid, '1 \n');
    count=fprintf(fid, '%d \n',tstep);
    count=fprintf(fid, '%d \n',tstep);  
    
    for i = Nidx(1)+1 : Nidx(end)
        
        if ~exist([dataname '_' int2str(i), '.F000',int2str(tstep)],'file')
            count=fprintf(fid, 'Y \n');
            count=fprintf(fid, 'U \n');
            count=fprintf(fid, '%s \n',[dataname '_' int2str(i)]);
            count=fprintf(fid, '1 \n');
            count=fprintf(fid, '1 \n');
            count=fprintf(fid, '%d \n',tstep);
            count=fprintf(fid, '%d \n',tstep);
        else
            count=fprintf(fid, 'N \n');
        end
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