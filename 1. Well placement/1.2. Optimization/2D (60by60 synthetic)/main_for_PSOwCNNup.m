w%%  Main for well placement opt. (PSOwCNNup)
%------------------------------------------%
%---------$   by Joonyi Kim   $------------%
%---------$      190529       $------------%
%---------$   for 2D(60by60)  $------------%
%------------------------------------------%

%%% Run simulation %%%
clear; addpath(pwd); addpath([pwd '/code']);
copyfile('data/*.*', pwd);
close all; warning off;

global Po Ciw Cpw discount_rate discount_term observed_term ... 
       Cw N N_ens N_iter N_scale ...
       nx ny dx dy...
       range ...
       area ...
       dtstep tstep pmax ...
       slstep nsteps ...
       maxtof maxp ...
       retrain iter ...
       direc_var  direc_fig ...
       posfile constfile

Po  = 60;
Ciw = 5;
Cpw = 3;
discount_rate = 0.1;
discount_term = 365;
observed_term = 30;
Cw = 2E+6;
N = 14;  % N. of wells
N_ens  = 10;
N_iter = 50;
N_scale = 3;    % upscaling ratio 
nx = 60;   % original scale
ny = 60;
dx = 120;
dy = 120;
area = 40;
range = [-4 4];

slstep = 30;            
nsteps = 1;
dtstep = 30;
tstep = 30;
pmax = 7200;

maxtof = 10000;
maxp   = 3500;

load 'PERMX5.mat';
load 'PERMX5_selected_idx.mat';

load 'trainedNet2218.mat'; % proxy model data

%% PSO
%%%% PSO parameter %%%%%
global Np w c1 c2 maxgen ...
       Nvar

Np = 40;
w  = 0.729;
c1 = 1.494;
c2 = 1.494;

maxgen = 200;
%%%%%%%%%%%%%%%%%%%%%%%%

casename  = '(perm5_5_cnn_up(test1))';

directory = ['reference' casename];
direc_var = ['variables' casename];
direc_fig = ['Figures' casename];
ecldata   = '2D_JY_Eclrun';
frsdata   = '2D_JY_Frsrun';
permfile  = '2D_PERMX';
posfile   = '2D_POSITION';
constfile = '2D_CONSTRAINT';

%%% CNN parameters %%%
global samp_std samp_mean Nsamp
samp_std = temp_std; samp_mean = temp_mean; Nsamp = size(pos_samp,1);
retrain = 1; 
% rtgen = maxgen + 1;  %% no retraining
rtgen = [40, 80, 120, 160];  %% retraining at certain generations

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now initialization');
mkdir(directory); copyfile(fullfile('*.DATA'), directory);
mkdir(direc_var); mkdir(direc_fig);   
copyfile('$convert.bat', directory);
    
gen = 1;

type = [];   % 1 - Prod., 0 - Injec., Empty - Random
pos = Initialize(type);
wset = [1500,1500, 5500,5500];  % Prod: BHP(psi), Injec: BHP(psi)

for j = 1:N_ens
    MakePermxFile(PERMX.original(:,selected(j)), [permfile '_' int2str(j) '.DATA']);
    copyfile([permfile '_' int2str(j) '.DATA'], directory);
end

Nvar = size(pos,2);

vel = zeros(Np, Nvar);
var_max = [repmat([nx;ny],N, 1); ones(N,1)];
var_min = [repmat([1;1], N, 1); 1/3*ones(1,1); -1*ones(N-1,1)];
dVar = var_max - var_min;
maxvel = dVar;

[ecl_fit,       ~,       ~, eclfitall, ~, ~] = Evaluate(pos, type, wset, directory, ecldata, permfile, []);
[pos_fit, pos_vio, pos_cpu, cnnfitall, ~, ~] = Evaluate(pos, type, wset, directory, frsdata, permfile, trainedNet);
figure; plotregression2(-eclfitall, -cnnfitall); 
print('-r600', '-dpng', [direc_fig '/regression_all #' num2str(gen-1) '.png']); close;
figure; plotregression2(-ecl_fit, -pos_fit);
print('-r600', '-dpng', [direc_fig '/regression #' num2str(gen-1) '.png']); close;

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
    h_fig = figure(2);
    h_par = plot(Np*N_ens*repmat(gen,1,length(find(pos_fit(:,1)<0))), -pos_fit((pos_fit(:,1)<0)),'or'); hold on;
    h_rep = plot(Np*N_ens*gen, -REP.pos_fit(:,1),'ok'); hold on;
    grid on; xlabel('iteration'); ylabel('NPV');
    axis([0 Np*N_ens*maxgen 0e+8 4.5e+8]);
    axis square;
    saveas(h_rep,[direc_fig '/generation #' num2str(gen-1) '.png']);

    display(['Generation #0 - best NPV: ' num2str(-REP.pos_fit(:,1))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen = 2;
stopCondition = false;
REP_POS = REP.pos; REP_fit = REP.pos_fit; 
POS_fit = pos_fit; POS_vio = pos_vio; POS = pos; POS_cpu = pos_cpu; 
ECL_fit = ecl_fit; ECLFIT = eclfitall; CNNFIT = cnnfitall; 
ECLTIME = []; CNNTIME = [];

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
        tic; [ecl_fit,       ~,       ~, eclfitall,   ~, ~] = Evaluate(pos, type, wset, directory, ecldata, permfile, []);         t1 = toc;
        tic; [pos_fit, pos_vio, pos_cpu, cnnfitall, TOF, P] = Evaluate(pos, type, wset, directory, frsdata, permfile, trainedNet); t2 = toc;
%         figure; h_reg = plotregression2(-eclfitall, -cnnfitall); 
        figure; h_reg = plotregression2(-ecl_fit, -pos_fit);
        print('-r600', '-dpng', [direc_fig '/regression #' num2str(gen-1) '.png']); close;

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
        end
        display(['Generation #' num2str(gen-1) '- best NPV: ' num2str(-REP.pos_fit(:,1))]);
        
        %} Update generation and check for termination
        gen = gen + 1;
        if(gen>maxgen), stopCondition = true; end
        fclose('all');
        try
        a{1,1} = whos;
        catch ME
            keyboard();
        end
        [REP]=SolutionVisualization(REP, PERMX.original(:,selected), ['Final solution ' casename '_for journal']);

        
        REP_POS = [REP_POS; REP.pos];
        REP_fit = [REP_fit; REP.pos_fit];
        POS_fit = [POS_fit, pos_fit];
        POS_vio = [POS_vio, pos_vio];
        ECL_fit = [ECL_fit, ecl_fit];
        POS_cpu = [POS_cpu; pos_cpu];
        POS     = [POS; pos];
        ECLFIT  = [ECLFIT, eclfitall];
        CNNFIT  = [CNNFIT, cnnfitall];
        ECLTIME = [ECLTIME; t1];
        CNNTIME = [CNNTIME; t2];
        
        %%%%% update proxy model %%%%%
        
        if ismember(gen-1, rtgen)

            disp('now retraining proxy model'); iter = 1;
            
            tp = repmat(pos(:,2*N+1:3*N), N_ens, 1);
            eclfitall_ = -(-eclfitall + Cw * sum(tp <= -1/3 | tp >= 1/3, 2));
            
            pos_samp     = [pos_samp; repmat(pos, N_ens,1)];
            pos_fit_samp = [pos_fit_samp; eclfitall_];
            TOF_samp{1}  = [TOF_samp{1}, TOF{1}];
            TOF_samp{2}  = [TOF_samp{2}, TOF{2}];
            P_samp{1}    = [P_samp{1}, P{1}];
            figure; h_reg = plotregression2(-ecl_fit, -pos_fit);
            print('-r600', '-dpng', [direc_fig '/regression before training #' num2str(retrain) '.png']); close;

            tic; [trainedNet, traininfo, samp_mean, samp_std] =  trainCNN(TOF_samp, P_samp, pos_samp, type, pos_fit_samp); t3=toc;
            [pos_fit_, pos_vio_, ~, ~, ~, ~] = Evaluate(pos, type, wset, directory, frsdata, permfile, trainedNet);
            figure; h_reg = plotregression2(-ecl_fit, -pos_fit_);
            print('-r600', '-dpng', [direc_fig '/regression after training #' num2str(retrain) '.png']); close;
            best = find(checkDomination(pos_fit_, pos_vio_)==0);
            
            rel_error = abs((ecl_fit(best) - pos_fit_(best))/ecl_fit(best));
            stopcondition = rel_error < 0.01;
            disp(['rel. error: ' num2str(rel_error)]);
            
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
                [~, ~, ~, eclfitall_,   ~, ~] = Evaluate(pos_, type, wset, directory, ecldata, permfile, []);         
                [~, ~, ~,          ~, TOF, P] = Evaluate(pos_, type, wset, directory, frsdata, permfile, trainedNet); 
                tp = repmat(pos_(:,2*N+1:3*N), N_ens, 1);
                eclfitall_ = -(-eclfitall_ + Cw * sum(tp <= -1/3 | tp >= 1/3, 2));
                
                pos_samp     = [pos_samp; repmat(pos_, N_ens,1)];
                pos_fit_samp = [pos_fit_samp; eclfitall_];
                TOF_samp{1}  = [TOF_samp{1}, TOF{1}];
                TOF_samp{2}  = [TOF_samp{2}, TOF{2}];
                P_samp{1}    = [P_samp{1}, P{1}];
                
                tic; [trainedNet, traininfo, samp_mean, samp_std] =  trainCNN(TOF_samp, P_samp, pos_samp, type, pos_fit_samp); t3=t3+toc;
                [pos_fit_, pos_vio_, ~, ~, ~, ~] = Evaluate(pos, type, wset, directory, frsdata, permfile, trainedNet); 
                figure; h_reg = plotregression2(-ecl_fit, -pos_fit_);
                print('-r600', '-dpng', [direc_fig '/regression after training #' num2str(retrain) '_' num2str(iter) '.png']); close;

                best = find(checkDomination(pos_fit_, pos_vio_)==0);
                               
                rel_error = abs((ecl_fit(best) - pos_fit_(best))/ecl_fit(best));
                stopcondition = rel_error < 0.01;
                disp(['rel. error: ' num2str(rel_error)]);
            end
            
            [REP.pos_fit, ~, ~, ~, ~, ~] = Evaluate(REP.pos, type, wset, directory, frsdata, permfile, trainedNet); 
            
            retrain = retrain + 1;
            
            save([direc_var '/trainedNet' num2str(retrain) '.mat'], 'TOF_samp', 'pos_fit_samp', 'pos_samp', 'temp_mean', 'temp_std', 'trainedNet', 'traininfo');
            
            disp(['validation rmse = ' num2str(traininfo.ValidationRMSE(end)) ', training time: ' num2str(t3) ' sec.']);
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tosave = cellfun(@isempty, regexp({a{1,1}.class}, '^matlab\.(ui|graphics)\.'));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'train[ei]')));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'TOF')));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'TOF_samp')));
        tosave = and(tosave, ~strcmp({a{1,1}.name}, 'P'));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'P_samp')));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'pos_samp')));
        tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'PERMX')));
        save([direc_var '/generation' num2str(gen-1) '.mat'], a{1,1}(tosave).name);

        disp(['Ecl time:' num2str(t1) ', Frs time:' num2str(t2)]);
        mdl = fitlm( -ecl_fit, -pos_fit);
        disp(['R square: ' num2str(mdl.Rsquared.Ordinary)]);


end

[ecl_fit_final,~,~,~,~,~] = Evaluate(REP.pos, type, wset, directory, ecldata, permfile, []);
pos_opt = REP.pos; pos_opt_fit = REP.pos_fit;
pos_opt_fit_real = ecl_fit_final;
save(['pos_opt' casename '.mat'], 'pos_opt', 'pos_opt_fit', 'pos_opt_fit_real');

rmpath(pwd);
warning on

