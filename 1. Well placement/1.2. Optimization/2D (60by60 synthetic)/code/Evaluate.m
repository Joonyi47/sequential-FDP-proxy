function [pos_fit, pos_vio, pos_tcpu, fitall, pos_tof, pos_p] = Evaluate(pos, type, contr, directory, datafile, permfile, net)
global  N ...
        N_ens ...
        nx ny ...
        discount_rate observed_term discount_term ...
        Cw Po Cpw Ciw ...
        samp_mean samp_std ...
        pmax dtstep ...
        slstep nsteps ...
        posfile constfile ...
        maxtof maxp


pos_fit  = [];
pos_vio  = [];
pos_tcpu = [];
pos_tof  = cell(1,2);
pos_p    = cell(1,1);

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
    fitness   = [];
    violation = [];
    tcpu      = [];
    
    for j = 1:Np
        Setwellpos([pos(j,:), type], [posfile '_' int2str(j) '.DATA']);
        Setwellcontrol([pos(j,:), type, contr], [constfile '_' int2str(j) '.DATA']);
        fclose('all');
    end
    if isempty(net)
        parfor j = 1:Np
            % get simulation result & evaluate fitness and violation
            dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' datafile '_' int2str(j) ' > NUL']);
        end
        for j = 1:Np
            [FOPT, FWPT, FWIT, TCPU] = GetProductiondata([datafile '_' int2str(j)], 3);
            g=(Po*FOPT(1)-Cpw*FWPT(1)-Ciw*FWIT(1))/((1+discount_rate)^(observed_term/discount_term));
            for k=1:size(FOPT,1)-1
                g=g+(Po*(FOPT(k+1)-FOPT(k))-Cpw*(FWPT(k+1)-FWPT(k))-Ciw*(FWIT(k+1)-FWIT(k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
            end
            g=g-Cw*sum(pos(j,2*N+1:end)<=-1/3 | pos(j,2*N+1:end)>=1/3);
            [fit] = -g;
            fitness = [fitness;fit];
            tcpu    = [tcpu;TCPU(end)];
        end
    else
        
        %%%%%%%%%%%%%%%%%%% streamline simulation %%%%%%%%%%%%%%%%%%%%%
        
        %%%% frontsim %%%%
        parfor j = 1:Np
            change_tstep([constfile '_' int2str(j)], nsteps, slstep/nsteps);
        end
%         tic
        parfor j = 1:Np
            [~,~] = dos(['C:\ecl\2009.1\bin\pc\frontsim.exe ' datafile '_' int2str(j) ' > NUL']);
        end
%         toc;
        for j = 1:Np
            [TCPU] = GetProductiondata2([datafile '_' int2str(j)], 3);
            tcpu   = [tcpu;TCPU(end)];
        end
        parfor j = 1:Np
            change_tstep([constfile '_' int2str(j)], pmax/dtstep, dtstep);
        end
        
        Restart_converter_writting_batch(datafile, nsteps, 1:Np);
        [~,~] = dos('$convert < Restart_converter.log > NUL');                              % º¯È¯½ÃÅ°±â.
        
        TOF_beg = [];
        TOF_end = [];
        PRESS   = [];
        
        parfor j= 1:Np
            
            TOF_beg = [TOF_beg, GetTOF([datafile '_' int2str(j)], nsteps, 'TIME_BEG')];
            TOF_end = [TOF_end, GetTOF([datafile '_' int2str(j)], nsteps, 'TIME_END')];
            PRESS   = [PRESS, GetTOF([datafile '_' int2str(j)], nsteps, 'PRESSURE')];
            
        end
        
        total_TOF = {TOF_beg, TOF_end};
         
        TOF_beg = TOF_beg;
        TOF_end = TOF_end;
        
        %%%% impes & ss %%%%
%         params = parameters;
%         
%         perm = GetPermx([permfile '.DATA']);
%         
%         P = [];
%         So = [];
%         
%         parfor j = 1:Np
%             
%             [p, so] = IMPES(pos(j,:), type(j,:), perm, params.swof, params);
%             P = [P, p(:,end)];
%             So = [So, so(:, end)];
%             
%         end
%         
%         TOF_beg = [];
%         TOF_end = [];
%         
%         parfor j = 1:Np
%             
%             [tbeg, tend] = StreamSim(pos(j,:), type(j,:), perm, params.swof, P(:,j), So(:,j), params);
%             TOF_beg = [TOF_beg, tbeg];
%             TOF_end = [TOF_end, tend];
%             
%         end
%         
%         parfor j = 1:Np
%             
%             TOF_beg(:,j) = reshape(reshape(TOF_beg(:,j), nx, ny)', nx*ny, 1);
%             TOF_end(:,j) = reshape(reshape(TOF_end(:,j), nx, ny)', nx*ny, 1);
%             
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        INPUT = [];
        
        for j = 1:Np
            
            INPUT(:,:,1,j) = (reshape(TOF_beg(:,j),nx,ny)'-maxtof/2)/maxtof;
            INPUT(:,:,2,j) = (reshape(TOF_end(:,j),nx,ny)'-maxtof/2)/maxtof;
            INPUT(:,:,3,j) = (reshape(PRESS(:,j),nx,ny)'-maxp)/maxp;
            INPUT(:,:,4,j) = 0*ones(nx,ny);
            
            if isempty(type)
                tp = pos(:,2*N+1:end);
            else
                tp = type;
            end

            for k = 1:N
                if tp(j,k) >= 1/3
                    INPUT(pos(j,2*k), pos(j,2*(k-1)+1), 4, j) = 1;
                elseif tp(j,k) <= -1/3
                    INPUT(pos(j,2*k), pos(j,2*(k-1)+1), 4, j) = -1;
                end
            end
            
        end
        
        fitness = predict(net, INPUT);
        fitness = -(fitness .* samp_std + samp_mean);
        fitness = -(-fitness - Cw * sum(tp <= -1/3 | tp >= 1/3, 2));
        
        pos_tof   = {[pos_tof{1}, total_TOF{1}], [pos_tof{2}, total_TOF{2}]};
        pos_p     = {[pos_p{1}, PRESS]};
        
        delete('*.S*', '*.F*', '*.X*');

    end

    for j = 1:Np
        [vio] = constViolation(pos(j,1:2*N), pos(j,2*N+1:end));
        violation = [violation;vio];
    end
    
    pos_fit   = [pos_fit, fitness];
    pos_vio   = [pos_vio, violation];
    pos_tcpu  = [pos_tcpu; tcpu];
end
fitall   = reshape(pos_fit, Np*N_ens, 1);
pos_fit  = mean(pos_fit,2);
pos_vio  = mean(pos_vio,2);
pos_tcpu = mean(pos_tcpu,2);

% delete all file except *.DATA, *.RSM
delete('*.PRT', '*.EGRID','*.UNSMRY', '*.SMSPEC','*.ECLEND', '*.GRID', '*INIT');
delete('*.S*', '*.F*', '*.X*');
fclose('all');
cd(home);




end
