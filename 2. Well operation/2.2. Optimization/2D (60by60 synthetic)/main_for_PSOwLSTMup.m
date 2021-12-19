%%  Main for well operation opt. (PSOwRNNup) 
%------------------------------------------%
%---------$   by Joonyi Kim   $------------%
%---------$      200120       $------------%
%---------$   for 2D(60by60)  $------------%
%------------------------------------------%

% based on the well placement opt. result

%%% Run simulation %%%

clear; addpath(pwd); addpath([pwd '/code']);
copyfile('data/*.*', pwd);
close all; warning off;
close (findall(groot, 'Type', 'Figure'));

global Po Ciw Cpw discount_rate discount_term observed_term ... 
       Cw Nmax N Nstep N_ens N_iter N_scale ...
       nx ny dx dy...
       range ...
       area ...
       dstep tstep pmax ...
       direc_var  direc_fig ...
       posfile constfile ...
       opt_pos opt_type ...
       retrain iter
   

Po  = 60;
Ciw = 5;
Cpw = 3;
discount_rate = 0.1;
discount_term = 365;
observed_term = 30;
Cw = 2E+06;
Nmax = 14;  % N. of max wells
N_ens  = 10;
N_iter = 50;
N_scale = 3;    % upscaling ratio 
nx = 60;   % original scale
ny = 60;
dx = 120;
dy = 120;
area = 40;
range = [-4 4];

dstep = 90;
tstep = 30;
pmax = 7200;
Nstep = pmax/dstep;

load 'PERMX5.mat';
load 'PERMX5_selected_idx.mat';
load 'pos_opt(5_5_cnn_up(13)).mat'; % result of well placement opt.
load 'trainedNet(5_5_cnn_up(13)).mat'; % LSTM proxy model

type_idx = (pos_opt(1,2*Nmax+1:end) >= 1/3 | pos_opt(1,2*Nmax+1:end) <= -1/3);
N = sum(type_idx);  % N. of drilled wells
opt_pos  = pos_opt(1, repelem(type_idx, 1, 2));
opt_type = pos_opt(1, logical([zeros(1, 2*Nmax), type_idx]));
opt_type(opt_type >= 1/3) = 1; opt_type(opt_type <= -1/3) = -1;

ref_type = zeros(N,4);
ref_type(opt_type == 1, 2) = 1/3; ref_type(opt_type == -1, 2) = 1-1/3;
ref_type = reshape(ref_type', 1, 4*N);
%% PSO
%%%% PSO parameter %%%%%
global Np w c1 c2 maxgen ...
       Nvar

Np = 40;
w  = 0.729;
c1 = 1.494;
c2 = 1.494;

maxgen = 70;

casename  = '(perm5_5_cnn_up(13)_rnn_up(test1))';

directory = ['reference' casename];
direc_var = ['variables' casename];
direc_fig = ['Figures' casename];
ecldata   = '2D_JY_Eclrun';
frsdata   = '2D_JY_Frsrun';
permfile  = '2D_PERMX';
posfile   = '2D_POSITION';
constfile = '2D_CONSTRAINT';

%%% RNN parameters %%%
global samp_std samp_mean Nsamp
samp_std = temp_std; samp_mean = temp_mean; Nsamp = N_ens*size(pos_samp,1);
retrain = 1;
rtgen = [15, 30, 45, 60];
% rtgen = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now initialization');
mkdir(directory); copyfile(fullfile('*.DATA'), directory);
mkdir(direc_var); mkdir(direc_fig);   
copyfile('$convert.bat', directory);
    
gen = 1;

pos = Initialize(opt_type);
wset = [1000,2500, 4500,6000];  % prod: pressure, injec: pressure

for j = 1:N_ens
    MakePermxFile(PERMX.original(:,selected(j)), [permfile '_' int2str(j) '.DATA']);
    copyfile([permfile '_' int2str(j) '.DATA'], directory);
end

Nvar = size(pos,2);
bound = repmat([0, 1; 0, 1; 0, pi; 0, pi], N, 1);

vel = zeros(Np,Nvar); 
var_max = bound(:,2);
var_min = bound(:,1);
dVar = var_max - var_min;
maxvel = dVar/5;

[opt_fit, ~, ~, ~, ~] = Evaluate(opt_pos, ref_type, wset, directory, ecldata, permfile, []);
display(['Opt NPV: ' num2str(-opt_fit)]);

[ecl_fit,       ~,       ~, eclfitall, raw] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, ecldata, permfile, []);
[pos_fit, pos_vio, pos_cpu, rnnfitall,   ~] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, frsdata, permfile, trainedNet);
figure; h_reg = plotregression2(-ecl_fit, -pos_fit);  
saveas(h_reg, [direc_fig '/regression #' num2str(gen-1) '.png']); close;

pbest = pos;
pbest_fit   = pos_fit;
pbest_vio   = pos_vio;

dominated   = checkDomination(pos_fit, pos_vio);
REP.pos     = pos(~dominated, :);
REP.pos_fit = pos_fit(~dominated,:);
REP.pos_vio = pos_vio(~dominated,:);


    %%%%% save initial variables %%%%%
    a=who;
    save([direc_var '/generation1.mat'], a{:});
    % number of infeasible solutions
    disp(['infeasible solutions: ' num2str(numel(find(pos_vio)>0))]);

    %} Plotting and verbose
    h_fig = figure();
    h_par = plot(Np*N_ens*repmat(gen,1,length(find(pos_fit(:,1)<0))), -pos_fit((pos_fit(:,1)<0)),'or'); hold on;
    h_rep = plot(Np*N_ens*gen, -REP.pos_fit(:,1),'ok'); hold on;
    grid on; xlabel('iteration'); ylabel('NPV');
    axis([0 Np*N_ens*maxgen 3e+8 5e+8]);
    axis square;
    saveas(h_rep,[direc_fig '/generation #' num2str(gen-1) '.png']);

    display(['Generation #0 - best NPV: ' num2str(-REP.pos_fit(:,1))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen = 2;
stopCondition = false;
REP_POS = REP.pos; REP_fit = REP.pos_fit; 
POS_fit = pos_fit; POS_vio = pos_vio; POS = pos; POS_cpu = pos_cpu; 
ECL_fit = ecl_fit; ECLFIT = eclfitall; RNNFIT = rnnfitall; 
ECLTIME = []; RNNTIME = [];

while ~stopCondition
    display(['now Generation #' num2str(gen-1)]);
        %} Select leader
        try
        h = selectNeighbor(pbest_fit(:,1), pbest_vio, REP);
        catch ME
        end
        
        %} Update speeds and positions
        vel = w.*vel + c1*rand(Np,Nvar).*(pbest-pos) ...
                     + c2*rand(Np,Nvar).*(repmat(REP.pos,Np,1)-pos);

        %} Check boundaries
        pos = pos + vel;
   
        %} Check boundaries
        [pos,vel] = checkBoundaries(pos,vel,maxvel,var_max,var_min);
                
        %} Evaluate the population
        tic; [ecl_fit,       ~,       ~, eclfitall, raw] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, ecldata, permfile, []);         t1 = toc;
        tic; [pos_fit, pos_vio, pos_cpu, rnnfitall,   ~] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, frsdata, permfile, trainedNet); t2 = toc;
        figure; h_reg = plotregression2(-ecl_fit, -pos_fit); 
        saveas(h_reg, [direc_fig '/regression #' num2str(gen-1) '.png']); close;

        disp(['infeasible solutions: ' num2str(numel(find(pos_vio)>0))]);
        
        %} Update the repository
        REP = updateRepository(REP,pos,pos_fit,pos_vio);
        REP.gen = gen;
        
        %} Update the best positions found so far for each particle
        pos_best = dominates([pos_fit(:,1), pos_vio], [pbest_fit(:,1), pbest_vio]);
        if(sum(pos_best)>1)
            pbest_fit(logical(pos_best),:) = pos_fit(logical(pos_best),:);
            pbest_vio(logical(pos_best),:) = pos_vio(logical(pos_best),:);
            pbest(logical(pos_best),:) = pos(logical(pos_best),:);
        end

        %} Plotting and verbose
        figure(h_fig);
        h_par = plot(Np*N_ens*repmat(gen,1,length(find(pos_fit(:,1)<0))), -pos_fit((pos_fit(:,1)<0)),'or'); hold on;
        h_rep = plot(Np*N_ens*gen, -REP.pos_fit(:,1),'ok'); hold on;
        if gen == maxgen
            saveas(h_rep,[direc_fig '/generation #' num2str(gen-1) '.png']);
            saveas(h_rep,[direc_fig '/generation #' num2str(gen-1) '.fig']);
        end
        display(['Generation #' num2str(gen-1) '- best NPV: ' num2str(-REP.pos_fit(:,1))]);
        
        %} Update generation and check for termination
        gen = gen + 1;
        if(gen>maxgen), stopCondition = true; end
        fclose('all');
                
        REP_POS = [REP_POS; REP.pos];
        REP_fit = [REP_fit; REP.pos_fit];
        POS_fit = [POS_fit, pos_fit];
        POS_vio = [POS_vio, pos_vio];
        ECL_fit = [ECL_fit, ecl_fit];
        POS_cpu = [POS_cpu; pos_cpu];
        POS     = [POS; pos];
        ECLFIT  = [ECLFIT, eclfitall];
        RNNFIT  = [RNNFIT, rnnfitall];
        ECLTIME = [ECLTIME; t1];
        RNNTIME = [RNNTIME; t2];

        %%%%% update proxy model %%%%%
        
        if ismember(gen-1, rtgen)

            disp('now retraining proxy model'); iter = 1;
            
            pos_samp     = [pos_samp; pos];
            pos_fit_samp = [pos_fit_samp; eclfitall];
            raw_m = cell(1,3);
            for i = 1:Np
                raw_m{1} = [raw_m{1}, mean(raw{1}(:,i:Np:end), 2)];
                raw_m{2} = [raw_m{2}, mean(raw{2}(:,i:Np:end), 2)];
                raw_m{3} = [raw_m{3}, mean(raw{3}(:,i:Np:end), 2)];
            end            
            RAW_samp{1}  = [RAW_samp{1}, raw_m{1}];
            RAW_samp{2}  = [RAW_samp{2}, raw_m{2}];
            RAW_samp{3}  = [RAW_samp{3}, raw_m{3}];
            figure; h_reg = plotregression2(-ecl_fit, -pos_fit);
            print('-r600', '-dpng', [direc_fig '/regression before training #' num2str(retrain) '.png']); close;
            
            tic; [trainedNet, traininfo, samp_mean, samp_std] =  trainRNN(RAW_samp, repmat(opt_pos, Np, 1), pos_samp, wset, pos_fit_samp); t3=toc;
            [pos_fit_, pos_vio_, ~, ~, ~] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, frsdata, permfile, trainedNet);
            figure; h_reg = plotregression2(-ecl_fit, -pos_fit_);
            print('-r600', '-dpng', [direc_fig '/regression after training #' num2str(retrain) '.png']); close;
            
            best = find(checkDomination(pos_fit_, pos_vio_)==0);
            
            rel_error = abs((ecl_fit(best) - pos_fit_(best))/ecl_fit(best));
            stopcondition = rel_error < 0.0015;
            disp(['rel. error: ' num2str(rel_error(1))]);
            
            while ~stopcondition
                
                iter = iter + 1;  disp(['iteration ' num2str(iter)]);
                
                %} Update speeds and positions
                vel_ = w.*vel + c1*rand(Np,Nvar).*(pbest-pos) ...
                    + c2*rand(Np,Nvar).*(repmat(REP.pos,Np,1)-pos);
                
                %} Check boundaries
                pos_ = pos + vel_;
                
                %} Check boundaries
                [pos_,vel_] = checkBoundaries(pos_,vel_,maxvel,var_max,var_min);
                
                %} Evaluate the population
                [~, ~, ~, eclfitall_, raw] = Evaluate(repmat(opt_pos, Np, 1), pos_, wset, directory, ecldata, permfile, []);         
                [~, ~, ~,          ~,   ~] = Evaluate(repmat(opt_pos, Np, 1), pos_, wset, directory, frsdata, permfile, trainedNet); 
                
                pos_samp     = [pos_samp; pos_];
                pos_fit_samp = [pos_fit_samp; eclfitall_];
                raw_m = cell(1,3);
                for i = 1:Np
                    raw_m{1} = [raw_m{1}, mean(raw{1}(:,i:Np:end), 2)];
                    raw_m{2} = [raw_m{2}, mean(raw{2}(:,i:Np:end), 2)];
                    raw_m{3} = [raw_m{3}, mean(raw{3}(:,i:Np:end), 2)];
                end                
                RAW_samp{1}  = [RAW_samp{1}, raw_m{1}];
                RAW_samp{2}  = [RAW_samp{2}, raw_m{2}];
                RAW_samp{3}  = [RAW_samp{3}, raw_m{3}];

                tic; [trainedNet, traininfo, samp_mean, samp_std] =  trainRNN(RAW_samp, repmat(opt_pos, Np, 1), pos_samp, wset, pos_fit_samp); t3=t3+toc;
                [pos_fit_, pos_vio_, ~, ~, ~] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, frsdata, permfile, trainedNet); 
                figure; h_reg = plotregression2(-ecl_fit, -pos_fit_);
                print('-r600', '-dpng', [direc_fig '/regression after training #' num2str(retrain) '_' num2str(iter) '.png']); close;

                best = find(checkDomination(pos_fit_, pos_vio_)==0);
                               
                rel_error = abs((ecl_fit(best) - pos_fit_(best))/ecl_fit(best));
                stopcondition = rel_error < 0.0015;
                disp(['rel. error: ' num2str(rel_error)]);
            end
            
            [REP.pos_fit, ~, ~, ~, ~] = Evaluate(opt_pos, REP.pos, wset, directory, frsdata, permfile, trainedNet);
            
            retrain = retrain + 1;
            
            save([direc_var '/trainedNet' num2str(retrain) '.mat'], 'RAW_samp', 'pos_fit_samp', 'pos_samp', 'temp_mean', 'temp_std', 'trainedNet', 'traininfo');
            
            disp(['validation rmse = ' num2str(traininfo.ValidationRMSE(end)) ', training time: ' num2str(t3) ' sec.']);
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        a{1,1} = whos;
        tosave = cellfun(@isempty, regexp({a{1,1}.class}, '^matlab\.(ui|graphics)\.'));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'train[ei]')));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'PERMX')));
        save([direc_var '/generation' num2str(gen-1) '.mat'], a{1,1}(tosave).name);

        disp(['Ecl time:' num2str(t1) ', Rnn time:' num2str(t2)]);
        mdl = fitlm( -ecl_fit, -pos_fit);
        disp(['R square: ' num2str(mdl.Rsquared.Ordinary)]);

end

[ecl_fit_final,~,~, ~,~] = Evaluate(opt_pos, REP.pos, wset, directory, ecldata, permfile, []);
pos_opt = REP.pos; pos_opt_fit = REP.pos_fit;
pos_opt_fit_real = ecl_fit_final;
save(['pos_opt' casename '.mat'], 'pos_opt', 'pos_opt_fit', 'pos_opt_fit_real');

rmpath(pwd);
warning on

