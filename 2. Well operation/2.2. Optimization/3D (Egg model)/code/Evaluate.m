function [pos_fit, pos_vio, pos_tcpu, fitall, pos_raw] = Evaluate(pos, type, contr, directory, datafile, permfile, net)
global  N Nstep ...
        N_ens ...
        nx ny ...
        discount_rate observed_term discount_term ...
        Cw Po Cpw Ciw ...
        samp_mean samp_std ...
        pmax dstep tstep...
        nsteps slstep ...
        posfile constfile ...
        ACTIVE ...
        opt_pos opt_type ...
        parameters


pos_fit  = [];
pos_vio  = [];
pos_tcpu = [];
pos_raw  = cell(1,3);

Np = size(pos,1);

%% for parallel computing %%
[home] =cd(directory);
if isempty(net)
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
            if i == row(2,1)+1
                fprintf(fid, '%s\n', [char(39)  permfile '_' int2str(p) '.DATA'  char(39) ' /']);
            elseif i == row(3,1)+1
                fprintf(fid, '%s\n', [char(39)  posfile '_' int2str(j) '.DATA'  char(39) ' /']);
            elseif i == row(4,1)+1
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
        Setwellpos([pos(j,:), opt_type], [posfile '_' int2str(j) '.DATA']);
        Setwellcontrol([pos(j,:),  type(j,:), contr], [constfile '_' int2str(j) '.DATA']);
        %         Setwellcontrol([pos(j,:), type], ['CONSTRAINT_' int2str(j) '.DATA']);
        fclose('all');
    end
        parfor j = 1:Np
            % get simulation result & evaluate fitness and violation
            dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' datafile '_' int2str(j) ' > NUL']);
        end
        for j = 1:Np
            [FOPT, FWPT, FWIT, TCPU] = GetProductiondata([datafile '_' int2str(j)], 3);
            pos_raw{1} = [pos_raw{1}, FOPT];
            pos_raw{2} = [pos_raw{2}, FWPT];
            pos_raw{3} = [pos_raw{3}, FWIT];
            
            g=(Po*FOPT(1)-Cpw*FWPT(1)-Ciw*FWIT(1))/((1+discount_rate)^(observed_term/discount_term));
            for k=1:size(FOPT,1)-1
                g=g+(Po*(FOPT(k+1)-FOPT(k))-Cpw*(FWPT(k+1)-FWPT(k))-Ciw*(FWIT(k+1)-FWIT(k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
            end
            tp = reshape(type(j,:), 4, N);
            nw = zeros(N ,Nstep);
            for i = 1:N
                nw(i, (tp(1,i) == 0) + ceil(Nstep * tp(1,i))) = 1;
            end            
            nw = sum(nw,1);
            for k =1:Nstep
                g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
            end
            [fit] = -g;
            fitness = [fitness;fit];
            tcpu    = [tcpu;TCPU(end)];
        end
        
        for j = 1:Np
            [vio] = constViolation(type(j,:));
            violation = [violation;vio];
        end
        
        pos_fit   = [pos_fit, fitness];
        pos_vio   = [pos_vio, violation];
        pos_tcpu  = [pos_tcpu; tcpu];
end

else
fitness   = [];
    violation = [];
    
    fun = @(x,y,z,b,t)x*(1+cos(y*(1-t)+z))+b;
    
    type_ = zeros(Np, N*Nstep);
    for i = 1:Np
        tmp = type(i,:);
        tmp = reshape(tmp, 4, N);
        schd = zeros(N, Nstep);
        for k = 1:N
            schd(k, (tmp(1,k) == 0)+ceil(Nstep * tmp(1,k)):end) = 1;
        end
        
        tmp3 = [];
        for j = 1:N
            if opt_type(j) >= 1/3
                tmp2 = feval(fun, (contr(2) - contr(1))/2*tmp(2,j), tmp(3,j), tmp(4,j), contr(1), [1/Nstep:1/Nstep:1]);
                tmp2 = tmp2.*schd(j,:);
%                 tmp2(tmp2 == 0) = contr(2)+100;
                tmp2(tmp2 == 0) = 500;
            else
                tmp2 = feval(fun, (contr(4) - contr(3))/2*tmp(2,j), tmp(3,j), tmp(4,j), contr(3), [1/Nstep:1/Nstep:1]);
                tmp2 = tmp2.*schd(j,:);
%                 tmp2(tmp2 == 0) = contr(3)-100;
                tmp2(tmp2 == 0) = 300;
            end
            tmp3 = [tmp3; tmp2];
        end
%         tmp3 = [tmp3; sum(schd,1)];
        
        type_(i,:) = reshape(tmp3', 1, numel(tmp3));
               
    end
    
    idx  = opt_type;
    idx  = (idx >= 1/3);   %% prod == 1, inj == 0
    
    data = cell(Np, 1);
    for i = 1:Np
        tmp = reshape(type_(i,:), Nstep, N)';
        data{i} = [tmp(idx == 1,:); tmp(idx == 0,:)];  %% sort (prod -> inj)
%         data{i} = tmp;  %% sort (prod -> inj)
    end
    
    %%% nomalize the data %%%
    X_  = data;
    for i = 1:numel(X_)
        X_{i} = (X_{i} - samp_mean{1}) ./ samp_std{1};
    end
    
    %%% predict FT by rnn %%%
    outp = predict(net, X_);
    for i = 1:numel(outp)
        outp{i} = outp{i} .* samp_std{2} + samp_mean{2};
    end
    Y_p = interpFT(outp); %%% interpolate output data for npv calculation
    
    for j = 1:Np
        g=(Po*Y_p{j}(1,1)-Cpw*Y_p{j}(2,1)-Ciw*Y_p{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
        for k=1:size(Y_p{1},2)-1
            g=g+(Po*(Y_p{j}(1,k+1)-Y_p{j}(1,k)) ...
                -Cpw*(Y_p{j}(2,k+1)-Y_p{j}(2,k)) ...
                -Ciw*(Y_p{j}(3,k+1)-Y_p{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
        end
        
        tp = reshape(type(j,:), 4, N);
        nw = zeros(N ,Nstep);
        for i = 1:N
            nw(i, (tp(1,i) == 0) + ceil((Nstep) * tp(1,i))) = 1;
        end
        nw = sum(nw,1);
        for k =1:Nstep
            g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
        end
        [fit] = -g;
        fitness = [fitness;fit];
    end
    
    for j = 1:Np
        [vio] = constViolation(type(j,:));
        violation = [violation;vio];
    end
    
    pos_fit  =  repmat(fitness, 1,N_ens);
    pos_vio   = repmat(violation, 1,N_ens);
end
fitall   = reshape(pos_fit,Np*N_ens,1);
pos_fit  = mean(pos_fit,2);
pos_vio  = mean(pos_vio,2);
pos_tcpu = mean(pos_tcpu,2);

% delete all file except *.DATA, *.RSM
delete('*.PRT', '*.EGRID','*.UNSMRY', '*.SMSPEC','*.ECLEND', '*.GRID', '*INIT');
delete('*.S*', '*.F*', '*.X*');
fclose('all');
cd(home);




end