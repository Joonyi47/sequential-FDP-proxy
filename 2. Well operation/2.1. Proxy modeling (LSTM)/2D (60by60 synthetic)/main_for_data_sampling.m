%%  Main for LSTM proxy 
%------------------------------------------%
%---------$   by Joonyi Kim   $------------%
%---------$      200120       $------------%
%---------$   for 2D(60by60)  $------------%
%------------------------------------------%

% based on the well placement opt. result

%%% Run simulation %%%
clear; addpath(pwd); addpath([pwd '/code']);
copyfile('data/*.*', pwd);
global Po Ciw Cpw discount_rate discount_term observed_term ... 
       Cw N Np Nmax Nstep N_ens N_iter N_scale ...
       nx ny dx dy...
       range ...
       area ...
       dstep tstep pmax ...
       opt_pos opt_type

%%% sample parameters %%%
Np = 500;
N_ens = 10;             % # of perm.fields

%%% field parameters %%%
Nmax = 14;              % # of wells
nx = 60;                % # of grids. x-direction
ny = 60;
dx = 120;
dy = 120;
area = 40;

%%% NPV parameters %%%
Po  = 60;
Ciw = 5;
Cpw = 3;
discount_rate = 0.1;
discount_term = 365;
observed_term = 30;
Cw = 2E+06;             % drilling cost

%%% simulation parameters %%%
tstep = 30;
dstep = 90;             % drilling term (step)
pmax = 7200;
Nstep = pmax/dstep;

%%% scaling parameters%%%
N_iter = 50;
N_scale = 3;    % upscaling ratio 
range = [-4 4];

load 'PERMX5.mat';
load 'PERMX5_selected_idx.mat';
load 'pos_opt(5_5_cnn_up(13)).mat';   % result of well placement opt.

type_idx = (pos_opt(1,2*Nmax+1:end) >= 1/3 | pos_opt(1,2*Nmax+1:end) <= -1/3);
N = sum(type_idx);  % N. of drilled wells
opt_pos  = pos_opt(1, repelem(type_idx, 1, 2));
opt_type = pos_opt(1, logical([zeros(1, 2*Nmax), type_idx]));
opt_type(opt_type >= 1/3) = 1; opt_type(opt_type <= -1/3) = -1;

%% make training set
%%%% parameters %%%%%
global Np ...
       direc_var  direc_fig ...
       posfile constfile

Np = 500;

directory = 'simulation';
direc_var = 'variables';
direc_fig = 'Figures';
ecldata   = '2D_JY_Eclrun';
permfile  = '2D_PERMX';
posfile   = '2D_POSITION';
constfile = '2D_CONSTRAINT';


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now initialization');
mkdir(directory); copyfile(fullfile('*.DATA'), directory);
mkdir(direc_var); mkdir(direc_fig);   
copyfile('$convert.bat', directory);
    
gen = 1;

for j = 1:1
   tic;  [pos{j}] = Initialize(opt_type); toc;
end
pos = repmat(pos, N_ens, 1)';
wset = [1000, 2500, 4500, 6000];  % prod: 1000~2500 psi, injec: 4500~6000 psi

pos_fit = []; pos_vio = []; pos_mat = []; pos_mat2 = []; tof_mat = []; total_TOF = cell(1,2); t_ecl = []; t_frs = []; total_raw = cell(1,3);
for j = 1:N_ens
    disp(['now iteration ' int2str(j)]);
    
    MakePermxFile(PERMX.original(:,selected(j)), [permfile '.DATA']);
    copyfile([permfile '.DATA'], directory);
    
    tic; [fit, vio, ~, raw] = Evaluate(repmat(opt_pos, Np, 1), pos{j}, wset, directory, ecldata, permfile, []); t1 = toc;
       
    pos_fit  = [pos_fit; fit];
    pos_vio  = [pos_vio; vio];
    t_ecl    = [t_ecl, t1];
    total_raw{1} = [total_raw{1}, raw{1}]; % FOPT
    total_raw{2} = [total_raw{2}, raw{2}]; % FWPT
    total_raw{3} = [total_raw{3}, raw{3}]; % FWIT

    disp(['Ecl time:' num2str(t1)]);
end

pos = cell2mat(pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath(pwd);

%%
t = datestr(now);
save(['Result_(pos_opt(5_5_cnn_up(13)))' t(1:11) ' ' t(13:14) t(16:17)]);
