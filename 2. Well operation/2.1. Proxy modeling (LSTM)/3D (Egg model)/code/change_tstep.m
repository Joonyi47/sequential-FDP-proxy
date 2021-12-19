function [] = change_tstep(filename, ntstep, dtstep)
global Np ...


fid = fopen([filename '.DATA'], 'r');
    total=100;
    if fid == -1
        disp('File open not successful')
    end
    Ecl=cell(300,300);
    i=1;
    while ~feof(fid)
        Ecl{i} = fgetl(fid);
        i=i+1;
    end
    row=find(strcmp(Ecl,'TSTEP')==1,total);
    idx = zeros(1,200);
    idx(1,row+1) = 1;
    
%     for j = 1:Np
%         fid = fopen([filename '_' int2str(j) '.DATA'], 'w');
        fid = fopen([filename '.DATA'], 'w');
        for i = 1:200
            if idx(i) == 1
                fprintf(fid, '%s\n', [int2str(ntstep) '*' int2str(dtstep) '/']);
            else
                fprintf(fid, '%s\n', Ecl{i});
            end
        end
        fclose(fid);
%     end
    
    fclose('all');

end