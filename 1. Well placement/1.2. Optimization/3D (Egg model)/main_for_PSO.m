%% Main for PSO (egg model) %%
%%% Run simulation %%%
addpath(pwd); close all;
global Po Ciw Cpw discount_rate discount_term observed_term ... 
       Cw N Nvar N_ens N_iter N_scale ...
       nx ny dx dy...
       range ...
       area ...
       dstep tstep pmax ...
       direc_var  direc_fig ...
       posfile constfile
       
Po  = 60;
Ciw = 5;
Cpw = 3;
discount_rate = 0.1;
discount_term = 365;
observed_term = 30;
Cw = 2E+6;
N = 12;  % N. of wells
Nvar = 3*N;
N_ens  = 10;
N_iter = 50;
N_scale = 3;    % upscaling ratio 
nx = 60;   % original scale
ny = 60;
dx = 36/0.3048;   % m
dy = 36/0.3048;
area = 40;
range = [-4 4];

dstep = 30;
tstep = 30;
pmax = 3600;

load 'ACTIVE.mat';
load 'PERM1.mat';
load 'PERM1_selected_idx.mat';

%% PSO
%%%% PSO parameter %%%%%
global Np w c1 c2 maxgen

Np = 40;
w  = 0.729;
c1 = 1.494;
c2 = 1.494;

maxgen = 200;

casename  = '(egg_perm1_6(4)_(7))'; 

directory = ['reference' casename];
direc_var = ['variables' casename];
direc_fig = ['Figures' casename];
datafile  = 'Egg_Model_Eclrun';
permfile  = 'Egg_Model_mDARCY';
posfile   = 'Egg_Model_POSITION';
constfile = 'Egg_Model_SCHEDULE';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now initialization');
mkdir(directory); copyfile(fullfile('*.DATA'), directory);
mkdir(direc_var); mkdir(direc_fig);   
copyfile('$convert.bat', directory);
    
gen = 1;

type = [];   % 1 - Prod., 0 - Injec.
pos = Initialize(type);
load('variables(egg_perm1_6(4)_cnn_up(15))/generation1.mat', 'pos');
wset = [300,300 500,500];  % Prod: BHP(barsa), Injec: BHP(barsa)

for j = 1:N_ens
    MakePermxFile(PERMX(:,selected(j)), [permfile '_' int2str(j) '.DATA']);
    copyfile([permfile '_' int2str(j) '.DATA'], directory);
end

vel = zeros(Np, 3*N);
var_max = [repmat([nx-1;ny-2],N, 1); ones(N,1)];
var_min = [repmat([2;2], N, 1); 1/3*ones(1,1); -1*ones(N-1,1)];
dVar = var_max - var_min;
maxvel = dVar;

[pos_fit, pos_vio, pos_cpu, eclfitall, ~,~] = Evaluate(pos, type, wset, directory, datafile, permfile, []);

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
    h_par = plot(Np*repmat(gen,1,length(find(pos_fit(:,1)<0))), -pos_fit((pos_fit(:,1)<0)),'or'); hold on;
    h_rep = plot(Np*gen, -REP.pos_fit(:,1),'ok'); hold on;
    grid on; xlabel('Simulations'); ylabel('NPV');
    axis([0 Np*maxgen 0e+8 1.2e+8]);
    axis square;
    saveas(h_rep,[direc_fig '/generation #' num2str(gen-1) '.png']);

    display(['Generation #0 - best NPV: ' num2str(-REP.pos_fit(:,1))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen = 2;
stopCondition = false;
REP_POS = REP.pos; REP_fit = REP.pos_fit; POS_fit = pos_fit; POS = pos;
POS_vio = pos_vio;
ECLFIT = eclfitall; 
POS_cpu = pos_cpu; TIME = [];

while ~stopCondition
    display(['now Generation #' num2str(gen-1)]);
        tic
        %} Select leader
        try
        h = selectNeighbor(pbest_fit(:,1), pbest_vio, REP);
        catch ME
        end
        
        %} Update speeds and positions
        vel = w.*vel + c1*rand(Np,3*N).*(pbest-pos) ...
                     + c2*rand(Np,3*N).*(repmat(REP.pos,Np,1)-pos);

        %} Update positions
        pos = pos + vel;
   
        %} Check boundaries
        [pos,vel] = checkBoundaries(pos,vel,maxvel,var_max,var_min);      
              
        %} Evaluate the population
        tic; [pos_fit, pos_vio, pos_cpu, eclfitall, ~,~] = Evaluate(pos, type, wset, directory, datafile, permfile, []); t=toc;
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
        h_par = plot(Np*repmat(gen,1,length(find(pos_fit(:,1)<0))), -pos_fit((pos_fit(:,1)<0)),'or'); hold on;
        h_rep = plot(Np*gen, -REP.pos_fit(:,1),'ok'); hold on;

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
        [REP]=SolutionVisualization(REP, PERMX(:,selected), ['Final solution ' casename]);
     
                
        REP_POS = [REP_POS; REP.pos];
        REP_fit = [REP_fit; REP.pos_fit];
        POS_fit = [POS_fit; pos_fit];
        POS_vio = [POS_vio, pos_vio];
        POS_cpu = [POS_cpu; pos_cpu];
        POS     = [POS; pos];
        ECLFIT  = [ECLFIT, eclfitall];
        TIME    = [TIME; t];
        
         tosave = cellfun(@isempty, regexp({a{1,1}.class}, '^matlab\.(ui|graphics)\.'));
         tosave = and(tosave, cellfun(@isempty, regexp({a{1,1}.name}, 'PERMX')));
        save([direc_var '/generation' num2str(gen-1) '.mat'], a{1,1}(tosave).name);
        
        disp(['Ecl time:' num2str(t)]);
end

pos_opt = REP.pos; pos_opt_fit = REP.pos_fit;
save(['pos_opt' casename '.mat'], 'pos_opt', 'pos_opt_fit');

rmpath(pwd);


