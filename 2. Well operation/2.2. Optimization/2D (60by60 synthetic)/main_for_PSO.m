%% Main for PSO %%
%%% Run simulation %%%

clear; close gcf
addpath(pwd); 
global Po Ciw Cpw discount_rate discount_term observed_term ... 
       Cw Nmax N Nstep N_ens N_iter N_scale ...
       nx ny dx dy...
       range ...
       area ...
       dstep tstep pmax ...
       direc_var  direc_fig ...
       posfile constfile ...
       opt_pos opt_type

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
load 'PERMX5_selected(narrow).mat';
load 'pos_opt(5_5_cnn_up(13)).mat';
% load 'pos_opt(base_case).mat';

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

casename  = '(perm5_5_cnn_up_(13)_ecl(test1))';   

directory = ['reference' casename];
direc_var = ['variables' casename];
direc_fig = ['Figures' casename];
datafile  = '2D_JY_Eclrun1';
permfile  = '2D_PERMX1';
posfile   = '2D_POSITION1';
constfile = '2D_CONSTRAINT1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now initialization');
mkdir(directory); copyfile(fullfile('*.DATA'), directory);
mkdir(direc_var); mkdir(direc_fig);   
copyfile('$convert.bat', directory);
    
gen = 1;

pos = Initialize(opt_type);
wset = [1000,2500, 4500,6000];  % prod: pressure, injec: pressure
% wset = [1500, 1000];  % prod: 2000 psi, injec: 2500 std/d

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

[opt_fit, ~, ~, ~, ~] = Evaluate(opt_pos, ref_type, wset, directory, datafile, permfile, []);
display(['Opt NPV: ' num2str(-opt_fit)]);
    
[pos_fit, pos_vio, pos_cpu, eclfitall, raw] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, datafile, permfile, []);

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
%     axis([0 Np*N_ens*maxgen 3.5e+8 4.5e+8]);
    axis([0 Np*N_ens*maxgen 3e+8 5e+8]);
    axis square;
    saveas(h_rep,[direc_fig '/generation #' num2str(gen-1) '.png']);

    display(['Generation #0 - best NPV: ' num2str(-REP.pos_fit(:,1))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen = 2;
stopCondition = false;
REP_POS = REP.pos; REP_fit = REP.pos_fit; POS_fit = pos_fit; POS_vio = pos_vio;
ECLFIT = eclfitall; POS = pos; POS_cpu = pos_cpu; TIME = [];

while ~stopCondition
    display(['now Generation #' num2str(gen-1)]);
        %} Select leader
        try
        h = selectNeighbor(pbest_fit(:,1), pbest_vio, REP);
        catch ME
        end
        
        %} Reinitialize
%         [pos, vel] = Reinitialize(pos, vel);
        
        %} Update speeds and positions
        vel = w.*vel + c1*rand(Np,Nvar).*(pbest-pos) ...
                     + c2*rand(Np,Nvar).*(repmat(REP.pos,Np,1)-pos);

        %} Check boundaries
        pos = pos + vel;
   
        %} Check boundaries
        [pos,vel] = checkBoundaries(pos,vel,maxvel,var_max,var_min);                
        
        
        %} Evaluate the population
        tic
        [pos_fit, pos_vio, pos_cpu, eclfitall, raw] = Evaluate(repmat(opt_pos, Np, 1), pos, wset, directory, datafile, permfile, []);
        t=toc;

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
        
        REP_POS = [REP_POS; REP.pos];
        REP_fit = [REP_fit; REP.pos_fit];
        POS_fit = [POS_fit, pos_fit];
        POS_vio = [POS_vio, pos_vio];
        POS_cpu = [POS_cpu; pos_cpu];
        POS     = [POS; pos];
        ECLFIT  = [ECLFIT, eclfitall];
        TIME    = [TIME; t];

        tosave = cellfun(@isempty, regexp({a{1,1}.class}, '^matlab\.(ui|graphics)\.'));
        save([direc_var '/generation' num2str(gen-1) '.mat'], a{1,1}(tosave).name);

        disp(['Ecl time:' num2str(t)]);


end

rmpath(pwd);


