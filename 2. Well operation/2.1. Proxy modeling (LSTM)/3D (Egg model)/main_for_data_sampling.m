%%  Main for LSTM proxy 
%------------------------------------------%
%---------$   by Joonyi Kim   $------------%
%---------$      200120       $------------%
%---------$   for Egg model   $------------%
%------------------------------------------%

% based on the well placement opt. result

%%% Run simulation %%%
clear; 
addpath(pwd);
global Po Ciw Cpw discount_rate discount_term observed_term ... 
       Cw N Nmax Nstep N_ens N_iter N_scale ...
       nx ny nz dx dy...
       range ...
       area ...
       tstep dstep pmax ...
       opt_pos opt_type
       
Po  = 60;
Ciw = 5;
Cpw = 3;
discount_rate = 0.1;
discount_term = 365;
observed_term = 30;
Cw = 2E+06;
Nmax = 12;  % # of wells (4 prod., 8 injec.)
N_ens  = 10;  % # of perm. fields
N_iter = 50;
N_scale = 3;    % upscaling ratio 
nx = 60;   % original scale
ny = 60;
nz = 2;
dx = 36/0.3048;
dy = 36/0.3048;
area = 40;
range = [-4 4];

tstep = 30;
dstep = 90;
pmax = 3600;
Nstep = pmax/dstep;

load 'ACTIVE.mat';
load 'PERM1.mat';
load 'PERM1_selected_idx.mat';
load 'pos_opt(1_6(4)_(7)).mat'; % result of well placement opt.

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
ecldata   = 'Egg_Model_Eclrun';   
permfile  = 'Egg_Model_mDARCY';
posfile   = 'Egg_Model_POSITION';
constfile = 'Egg_Model_SCHEDULE';


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now initialization');
mkdir(directory); copyfile(fullfile('*.DATA'), directory);
mkdir(direc_var); mkdir(direc_fig);   
copyfile('$convert.bat', directory);
    
gen = 1;
for j = 1:1
    [pos{j}] = Initialize(opt_type);
end
pos = repmat(pos, N_ens, 1)';
wset = [250, 350, 450, 550];  % Prod: BHP(barsa), Injec:BHP(barsa)


%%
pos_fit = []; pos_vio = []; pos_tcpu = []; pos_mat = []; pos_mat2 = []; tof_mat = []; total_TOF = cell(1,2); t_ecl = []; t_frs = []; total_P = cell(1,1);
total_raw = cell(1,3);
for j = 1:N_ens 
    disp(['now iteration ' int2str(j)]);
    
    MakePermxFile(PERMX(:,selected(j)), [permfile '.DATA']);
    copyfile([permfile '.DATA'], directory);
    tic;    [fit, vio, tcpu, raw] = Evaluate(repmat(opt_pos, Np, 1), pos{j}, wset, directory, ecldata, permfile, []); t1 = toc;
        
    pos_fit  = [pos_fit; fit];
    pos_vio  = [pos_vio; vio];
    pos_tcpu = [pos_tcpu; tcpu];
    t_ecl    = [t_ecl, t1];
    total_raw{1} = [total_raw{1}, raw{1}];
    total_raw{2} = [total_raw{2}, raw{2}];
    total_raw{3} = [total_raw{3}, raw{3}];
    
    disp(['Ecl time:' num2str(t1)]);
end

pos = cell2mat(pos');

rmpath(pwd);
%%
t = datestr(now);
save(['Result_' t(1:11) ' ' t(13:14) t(16:17)]);