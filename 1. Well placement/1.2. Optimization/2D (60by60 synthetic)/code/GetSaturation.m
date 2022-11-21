function [SOIL] = GetSaturation(pos, type, contr, directory, datafile, permfile, net)
global  N_ens ...
        pmax tstep...
        posfile constfile ...


SOIL = [];

Np = size(pos,1);

%% for parallel computing %%
[home] =cd(directory);
for p = 1:N_ens
    fid = fopen([datafile '.DATA'], 'r');
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
    row=find(strcmp(Ecl,'INCLUDE')==1,total);
    for j = 1:Np
        fid = fopen([datafile '_' int2str(j) '.DATA'], 'w');
        for i = 1:i
            if i == row(1,1)+1
                fprintf(fid, '%s\n', [char(39)  permfile '_' int2str(p) '.DATA'  char(39) ' /']);
            elseif i == row(2,1)+1
                fprintf(fid, '%s\n', [char(39)  posfile '_' int2str(j) '.DATA'  char(39) ' /']);
            elseif i == row(3,1)+1
                fprintf(fid, '%s\n', [char(39)  constfile '_' int2str(j) '.DATA' char(39) ' /']);
            else
                fprintf(fid, '%s\n', Ecl{i});
            end
        end
        fclose(fid);
    end
    
    fclose('all');
    
    %%  Evaluate
    
    for j = 1:Np
        Setwellpos([pos(j,:), type], [posfile '_' int2str(j) '.DATA']);
        Setwellcontrol([pos(j,:), type, contr], [constfile '_' int2str(j) '.DATA']);
        %         Setwellcontrol([pos(j,:), type], ['CONSTRAINT_' int2str(j) '.DATA']);
        fclose('all');
    end
%         tic
        parfor j = 1:Np
            % get simulation result & evaluate fitness and violation
            dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' datafile '_' int2str(j) ' > NUL']);
        end
%         toc;
        delete('*.F*');
        Restart_converter_writting([datafile '_1'], pmax/tstep, pmax/tstep);
        [~,~] = dos('$convert < Restart_converter.log > NUL');
        soil = 1 - [GetTOF([datafile '_1'], pmax/tstep, 'SWAT')];
        delete('*.X*', '*.S*', '*.F*');
        %         fitness = [fitness;-FOPT(end)];
        
    SOIL= [SOIL, soil];
    
end

% delete all file except *.DATA, *.RSM
delete('*.PRT', '*.EGRID','*.UNSMRY', '*.SMSPEC','*.ECLEND', '*.GRID', '*INIT');
delete('*.S*', '*.F*', '*.X*');
fclose('all');
cd(home);

end