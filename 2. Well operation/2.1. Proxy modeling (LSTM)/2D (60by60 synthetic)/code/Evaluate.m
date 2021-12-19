function [pos_fit, pos_vio, fitness, rawdata] = Evaluate(pos, type, wset, directory, datafile, permfile, net)
    global N Nstep pmax...
           discount_rate observed_term discount_term ...
           Cw Po Cpw Ciw ...
           posfile constfile ...
           opt_type
    
      
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
                fprintf(fid, '%s\n', [char(39)  permfile '.DATA'  char(39) ' /']);
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
    fitness   = [];
    violation = [];
    rawdata   = cell(1,3);
    
    for j = 1:Np
        Setwellpos([pos(j,:), opt_type], [posfile '_' int2str(j) '.DATA']);
        Setwellcontrol([pos(j,:), type(j,:), wset], [constfile '_' int2str(j) '.DATA']);
        fclose('all');
    end
    if isempty(net)
        parfor j = 1:Np
            % get simulation result & evaluate fitness and violation
            dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' datafile '_' int2str(j) ' > NUL']);
        end
        for j = 1:Np
            [FOPT, FWPT, FWIT, ~] = GetProductiondata([datafile '_' int2str(j)], 3);
            rawdata{1} = [rawdata{1}, FOPT]; 
            rawdata{2} = [rawdata{2}, FWPT];
            rawdata{3} = [rawdata{3}, FWIT];
            
            g=(Po*FOPT(1)-Cpw*FWPT(1)-Ciw*FWIT(1))/((1+discount_rate)^(observed_term/discount_term));
            for k=1:size(FOPT,1)-1
                g=g+(Po*(FOPT(k+1)-FOPT(k))-Cpw*(FWPT(k+1)-FWPT(k))-Ciw*(FWIT(k+1)-FWIT(k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
            end
            
            tp = reshape(type(j,:), 4, N);            
            nw = zeros(N ,Nstep);
            for i = 1:N
                nw(i, ceil((Nstep) * tp(1,i))) = 1;
            end            
            nw = sum(nw,1);
            for k =1:Nstep
                g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
            end
            [fit] = -g;
            fitness = [fitness;fit];            
        end
    else
        fitness = net(pos');
        fitness = fitness';
    end
    for j = 1:Np
        [vio] = constViolation(type(j,:));
        violation = [violation;vio];
    end
    pos_fit   = fitness;
    pos_vio   = violation;
    % delete all file except *.DATA, *.RSM
    delete('*.PRT', '*.EGRID','*.UNSMRY', '*.SMSPEC','*.ECLEND');
    fclose('all');
    cd(home);
    
    


end