addpath(pwd);

datafile  = 'Egg_Model_ECL';
directory = 'sim_for_NPVrange';
copyfile('*.DATA', directory);
copyfile('*.mat', directory);

[home] = cd(directory);
% load PERM_toLy2.mat;
% PERMX = PERMX_toLy2;
load PERM1.mat;

fitness   = [];
tcpu      = [];

for ii = 1:100
    MakePermxFile(PERMX(:,ii), 'mDARCY.DATA');
    dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' datafile ' > NUL']);

        [FOPT, FWPT, FWIT, TCPU] = GetProductiondata(datafile, 3);
        g=(Po*FOPT(1)-Cpw*FWPT(1)-Ciw*FWIT(1))/((1+discount_rate)^(observed_term/discount_term));
        for k=1:size(FOPT,1)-1
            g=g+(Po*(FOPT(k+1)-FOPT(k))-Cpw*(FWPT(k+1)-FWPT(k))-Ciw*(FWIT(k+1)-FWIT(k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
        end
        g=g-Cw*N;
        [fit] = -g;
        fitness = [fitness;fit];
        tcpu    = [tcpu;TCPU(end)];
  
    
end

selected = [];
fit_sorted = sort(-fitness, 'ascend');
for ii = 1:10
   [~,I] = min(abs(fit_sorted - quantile(fit_sorted, 0.05+0.1*(ii-1))));
   selected = [selected, I];
end

cd(home);
rmpath(pwd);