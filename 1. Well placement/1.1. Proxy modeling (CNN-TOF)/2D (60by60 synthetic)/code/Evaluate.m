function [pos_fit, pos_vio, fitness] = Evaluate(pos, type, wset, directory, datafile, permfile, net)
    global N ...
           discount_rate observed_term discount_term ...
           Cw Po Cpw Ciw ...
           posfile constfile 
           
    Np = size(pos,1);
    
    pos_fit = [];
    pos_vio = [];
    
    %% for parallel computing %%
    [home] =cd(directory);
    
    fid = fopen([datafile '.DATA'], 'r');
    total=100;
    if fid == -1
        disp('File open not successful')
    end
    Ecl=[];
    i=1;
    while ~feof(fid)
        Ecl{i} = fgetl(fid);
        i=i+1;
    end
    row=find(strcmp(Ecl,'INCLUDE')==1,total);
    for j = 1:Np
        fid = fopen([datafile '_' int2str(j) '.DATA'], 'w');
        for i = 1:size(Ecl,2)
            if i == row(1)+1
                fprintf(fid, '%s\n', [char(39)  permfile '.DATA'  char(39) ' /']);
            elseif i == row(2)+1
                fprintf(fid, '%s\n', [char(39)  posfile '_' int2str(j) '.DATA'  char(39) ' /']);
            elseif i == row(3)+1
                fprintf(fid, '%s\n', [char(39)  constfile '_' int2str(j) '.DATA' char(39) ' /']);
            else
                fprintf(fid, '%s\n', Ecl{i});
            end
        end
        fclose(fid);
    end
    
    fclose('all');
    
    %%  Evaluate
    fitness   = [];
    violation = [];
    for j = 1:Np
        Setwellpos([pos(j,:), type], [posfile '_' int2str(j) '.DATA']);
        Setwellcontrol([pos(j,:), type, wset], [constfile '_' int2str(j) '.DATA']);
        fclose('all');
    end
    if isempty(net)
        parfor j = 1:Np
            % get simulation result & evaluate fitness and violation
            dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' datafile '_' int2str(j) ' > NUL']);
        end
        for j = 1:Np
            [FOPT, FWPT, FWIT, ~] = GetProductiondata([datafile '_' int2str(j)], 3);
            g=(Po*FOPT(1)-Cpw*FWPT(1)-Ciw*FWIT(1))/((1+discount_rate)^(observed_term/discount_term));
            for k=1:size(FOPT,1)-1
                g=g+(Po*(FOPT(k+1)-FOPT(k))-Cpw*(FWPT(k+1)-FWPT(k))-Ciw*(FWIT(k+1)-FWIT(k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
            end
            g=g-Cw*sum(pos(j,2*N+1:end) <= -1/3 | pos(j,2*N+1:end) >= 1/3);
            [fit] = -g;
            fitness = [fitness;fit];
        end
        %         fitness = [fitness;-FOPT(end)];
    else
        fitness = net(pos');
        fitness = fitness';
    end
    for j = 1:Np
        [vio] = constViolation(pos(j,1:2*N), pos(j,2*N+1:end));
        violation = [violation;vio];
    end
    pos_fit   = fitness;
    pos_vio   = violation;
    % delete all file except *.DATA, *.RSM
    delete('*.PRT', '*.EGRID','*.UNSMRY', '*.SMSPEC','*.ECLEND');
    fclose('all');
    cd(home);
    
    


end