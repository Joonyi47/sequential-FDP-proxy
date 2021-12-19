%% Run simulations
Po  = 65;
Ciw = 8;
Cpw = 5;
discount_rate = 0;
discount_term = 365;
observed_term = 30;
Cw = 20E+06;
N = 5;  % N. of wells
N_ens = 1000;

load 'PERMX.mat';

%% original
fit = []; FOPT = []; FWPT = []; FWIT = []; TCPU = []; time = [];
for i = 1:N_ens
    filename = '2D_JY_eclrun1';
    MakePermxFile(PERMX.original(:,i), '2D_PERMX1.DATA');
    tic
    dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' filename ' > NUL']);
    t = toc
    i
    [FOPT(:,i), FWPT(:,i), FWIT(:,i), TCPU(:,i)] = GetProductiondata(filename, 3);
    g=(Po*FOPT(1)-Cpw*FWPT(1)-Ciw*FWIT(1))/((1+discount_rate)^(observed_term/discount_term));
    for j=1:size(FOPT,1)-1
      g=g+(Po*(FOPT(j+1,i)-FOPT(j,i))-Cpw*(FWPT(j+1,i)-FWPT(j,i))-Ciw*(FWIT(j+1,i)-FWIT(j,i)))/((1+discount_rate)^(observed_term*(j+1)/discount_term));
    end
    g=g-Cw*N;
    
    fit = [fit, g];
    time = [time, t];
end

%% upscale 1
fit1 = []; FOPT1 = []; FWPT1 = []; FWIT1 = []; TCPU1 = []; time1 = [];
for i = 1:N_ens
    filename = '2D_JY_eclrun2';
    MakePermxFile(PERMX.ups1(:,i), '2D_PERMX2.DATA');
    tic
    dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' filename ' > NUL']);
    t1 = toc;
    [FOPT1(:,i), FWPT1(:,i), FWIT1(:,i), TCPU1(:,i)] = GetProductiondata(filename, 3);
    g=(Po*FOPT1(1)-Cpw*FWPT1(1)-Ciw*FWIT1(1))/((1+discount_rate)^(observed_term/discount_term));
    for j=1:size(FOPT1,1)-1
      g=g+(Po*(FOPT1(j+1,i)-FOPT1(j,i))-Cpw*(FWPT1(j+1,i)-FWPT1(j,i))-Ciw*(FWIT1(j+1,i)-FWIT1(j,i)))/((1+discount_rate)^(observed_term*(j+1)/discount_term));
    end
    g=g-Cw*N;
    
    fit1 = [fit1, g];
    time1 = [time1, t1];
end
%% upscale 2
fit2 = []; FOPT2 = []; FWPT2 = []; FWIT2 = []; TCPU2 = []; time2 = [];
for i = 1:N_ens
    filename = '2D_JY_eclrun3';
    MakePermxFile(PERMX.ups2(:,i), '2D_PERMX3.DATA');
    tic
    dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' filename ' > NUL']);
    t2 = toc;
    [FOPT2(:,i), FWPT2(:,i), FWIT2(:,i), TCPU2(:,i)] = GetProductiondata(filename, 3);
    g=(Po*FOPT2(1)-Cpw*FWPT2(1)-Ciw*FWIT2(1))/((1+discount_rate)^(observed_term/discount_term));
    for j=1:size(FOPT2,1)-1
      g=g+(Po*(FOPT2(j+1,i)-FOPT2(j,i))-Cpw*(FWPT2(j+1,i)-FWPT2(j,i))-Ciw*(FWIT2(j+1,i)-FWIT2(j,i)))/((1+discount_rate)^(observed_term*(j+1)/discount_term));
    end
    g=g-Cw*N;
    
    fit2 = [fit2, g];
    time2 = [time2, t2];
end



