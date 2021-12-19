%%  Main for  CNN proxy 
%------------------------------------------%
%---------$   by Joonyi Kim   $------------%
%---------$      190529       $------------%
%---------$   for Egg models  $------------%
%------------------------------------------%

%%% Run simulation %%%
clear; addpath(pwd); addpath([pwd '/code']);
copyfile('data/*.*', pwd);
global Po Ciw Cpw discount_rate discount_term observed_term ... 
       Cw N N_ens N_iter ...
       nx ny nz dx dy...
       area ...
       dtstep pmax ...
       slstep nsteps ...
       ACTIVE
       
Po  = 60;
Ciw = 5;
Cpw = 3;
discount_rate = 0.1;
discount_term = 365;
observed_term = 30;
Cw = 0E+06;
N = 12;  % # of wells (4 prod., 8 injec.)
N_ens  = 10;  % # of perm. fields
N_iter = 50;
nx = 60;   % original scale
ny = 60;
nz = 2;
dx = 36/0.3048;
dy = 36/0.3048;
area = 40;

slstep = 30;
nsteps = 1;
dtstep = 30;
pmax = 3600;

load 'ACTIVE.mat';
load 'PERM1.mat';
load 'PERM1_selected_idx.mat';

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
frsdata   = 'Egg_Model_Frsrun';   % frontsim datafile name
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
for j = 1:N_ens
    pos{j} = Initialize(type);
end
wset = [300, 500];  % Prod. BHP(barsa), Inj. BHP(barsa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Simulate sample data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
pos_fit = []; pos_vio = []; pos_tcpu = []; pos_mat = []; pos_mat2 = []; tof_mat = []; total_TOF = cell(1,2); t_ecl = []; t_frs = []; total_P = cell(1,1);
for j = 1:N_ens 
    disp(['now iteration ' int2str(j)]);
    
    MakePermxFile(PERMX(:,selected(j)), [permfile '.DATA']);
    copyfile([permfile '.DATA'], directory);
    
    tic; [fit, vio, tcpu]  = Evaluate(pos{j}, type, wset, directory, ecldata, permfile, []); t1 = toc;
    tic; [mat1, tof, p, ~] = SLsimulate(pos{j}(:,1:2*N), pos{j}(:,2*N+1:3*N), directory, frsdata, permfile, j); t2 = toc;
    
    pos_fit  = [pos_fit; fit];
    pos_vio  = [pos_vio; vio];
    pos_tcpu = [pos_tcpu; tcpu];
    pos_mat  = [pos_mat, mat1];
    t_ecl    = [t_ecl, t1];
    t_frs    = [t_frs, t2];

    total_TOF{1} = [total_TOF{1}, tof{1}];
    total_TOF{2} = [total_TOF{2}, tof{2}];
    total_P{1}   = [total_P{1}, p{1}];
    
    disp(['Ecl time:' num2str(t1) ', Frs time:' num2str(t2)]);
end

nump = length(find(type == 1));
numi = length(find(type == 0));

Idx_act = find(ACTIVE == 1);
TOF_beg = 30000*ones(nx*ny*2, Np*N_ens);
TOF_beg(Idx_act,:) = total_TOF{1};
TOF_end = 30000*ones(nx*ny*2, Np*N_ens);
TOF_end(Idx_act,:) = total_TOF{2};
PRESS = 400*ones(nx*ny*2, Np*N_ens);
PRESS(Idx_act,:) = total_P{1};
total_TOF{1} = TOF_beg; total_TOF{2} = TOF_end; 
total_P{1} = PRESS;
pos = cell2mat(pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmpath(pwd);
%%
t = datestr(now);
save(['Result_' t(1:11) ' ' t(13:14) t(16:17)]);