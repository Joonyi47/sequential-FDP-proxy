addpath([pwd '/code']);
%% Data preprocessing

fitness = [];
for j = 1:size(pos_fit,1)
    g=(Po*total_raw{1}(1,j)-Cpw*total_raw{2}(1,j)-Ciw*total_raw{3}(3,j))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(total_raw{1},1)-1
        g=g+(Po*(total_raw{1}(k+1,j)-total_raw{1}(k,j)) ...
            -Cpw*(total_raw{2}(k+1,j)-total_raw{2}(k,j)) ...
            -Ciw*(total_raw{3}(k+1,j)-total_raw{3}(k,j)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness = [fitness;fit];
end

fitness_m = mean(reshape(fitness, Np, N_ens),2);

%%% 평균 production (or injection) total

total_raw_m = cell(1,3);
for j = 1:Np
    total_raw_m{1} = [total_raw_m{1}, mean(total_raw{1}(:,j:Np:end), 2)];
    total_raw_m{2} = [total_raw_m{2}, mean(total_raw{2}(:,j:Np:end), 2)];
    total_raw_m{3} = [total_raw_m{3}, mean(total_raw{3}(:,j:Np:end), 2)];
end

fitness_m_ = [];
for j = 1:size(total_raw_m{1},2)
    g=(Po*total_raw_m{1}(1,j)-Cpw*total_raw_m{2}(1,j)-Ciw*total_raw_m{3}(3,j))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(total_raw_m{1},1)-1
        g=g+(Po*(total_raw_m{1}(k+1,j)-total_raw_m{1}(k,j)) ...
            -Cpw*(total_raw_m{2}(k+1,j)-total_raw_m{2}(k,j)) ...
            -Ciw*(total_raw_m{3}(k+1,j)-total_raw_m{3}(k,j)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_m_ = [fitness_m_;fit];
end


RAW = cell(1,3);
for i = 1:Np
    RAW{1} = [RAW{1}, mean(total_raw{1}(:,i:Np:end), 2)];
    RAW{2} = [RAW{2}, mean(total_raw{2}(:,i:Np:end), 2)];
    RAW{3} = [RAW{3}, mean(total_raw{3}(:,i:Np:end), 2)];
end


fun = @(x,y,z,b,t)x*(1+cos(y*(1-t)+z))+b;

pos2 = zeros(Np, N*Nstep);
for i = 1:Np
    tmp = pos(i,:);
    tmp = reshape(tmp, 4, N);
    schd = zeros(N, Nstep);
    for k = 1:N
        schd(k, ceil(Nstep * tmp(1,k)):end) = 1;
    end
    
    tmp3 = []; tmp2= [];
   for j = 1:N
       if opt_type(j) >= 1/3
           tmp2 = feval(fun, (wset(2) - wset(1))/2*tmp(2,j), tmp(3,j), tmp(4,j), wset(1), [1/Nstep:1/Nstep:1]);
           tmp2 = tmp2.*schd(j,:);
           tmp2(tmp2 == 0) = 400; % null value (initial reservoir pressure in this case)
       else
           tmp2 = feval(fun, (wset(4) - wset(3))/2*tmp(2,j), tmp(3,j), tmp(4,j), wset(3), [1/Nstep:1/Nstep:1]);
           tmp2 = tmp2.*schd(j,:);
           tmp2(tmp2 == 0) = 400;
       end
       tmp3 = [tmp3; tmp2];
   end
   
   pos2(i,:) = reshape(tmp3', 1, N*Nstep);
   
end


%% Initial conditions
case_sample = 1;
case_sensitivity = 1;

Nall = 500;
Ntrain = ceil(Nall*0.7);
Ntest = Nall * 0.15;
Nsample = Nall - Ntest;

selected = randperm(Nall, Nsample);
selected = sort(selected, 'ascend');

aa = randperm(Nsample,Ntrain);
tr = selected(aa);
val = [];
for iter = selected
   if find(iter==tr)
      continue
   end
    val = [val,iter];
end

nonselected = [];
for j = 1:Nall
   if ~ismember(j, selected)
       nonselected = [nonselected, j];
   end
end

DATA   = cell(Nall, 1);
idx    = opt_type;
idx    = (idx >= 1/3);

a = 1/3+(1-1/3)/4;
for i = 1:Nall   
    tmp = reshape(pos2(i,:), Nstep, N)';
%     tmp = reshape(pos(i,:), Nstep, N)';
    DATA{i} = [tmp(idx == 1,:); tmp(idx ==0,:)];
end

yDATA        = cell(Nall, 1);
pstep        = Nstep;

for i = 1:Nall
   yDATA{i} =  [total_raw_m{1}(pmax/tstep/pstep:pmax/tstep/pstep:end,i)'; ...
                total_raw_m{2}(pmax/tstep/pstep:pmax/tstep/pstep:end,i)'; ...
                total_raw_m{3}(pmax/tstep/pstep:pmax/tstep/pstep:end,i)'];
end

varname = ['Result_(PERM1_6(4)_(7))_' t(1:11) ' ' t(13:14) t(16:17)];   %%% streamline simulation 결과 변수
dirname = [varname '/' 'Nsample ' int2str(Nsample) '/sample case ' int2str(case_sample) '/sensitivity case ' int2str(case_sensitivity)];
mkdir(dirname);

%% Input data
case_parameter = '1';
mkdir([dirname '/' case_parameter]);

X_total = DATA([tr, val]);
Y_total = yDATA([tr, val]);

X_train = DATA(tr);
mu = mean([X_total{:}],2);
sig = std([X_total{:}],0,2);
tmp2 = [X_total{:}];

for i = 1:numel(X_train)
    X_train{i} = (X_train{i} - mu) ./ sig;
end
X_validation = DATA(val);
for i = 1:numel(X_validation)
    X_validation{i} = (X_validation{i} - mu) ./ sig;
end

Y_total_ = [];
for i = 1:numel(Y_total)
   Y_total_(:,:,i) = Y_total{i};
end
mu2 = mean(Y_total_,3);
sig2 = std(Y_total_,0,3);

Y_train = yDATA(tr);
for i = 1:numel(Y_train)
    Y_train{i} = (Y_train{i} - mu2) ./ sig2;
    Y_train{i}(isnan(Y_train{i})) = 0;
end
Y_validation = yDATA(val);
for i = 1:numel(Y_validation)
    Y_validation{i} = (Y_validation{i} - mu2) ./ sig2;
    Y_validation{i}(isnan(Y_validation{i})) = 0;
end

%% LSTM model
numFeatures = N;
numHiddenUnits1 = 125;
numHiddenUnits2 = 100;
numHiddenUnits3 = 75;
numResponses = size(yDATA{1},1);
layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits1,'OutputMode','sequence')
    dropoutLayer(0.3)
    lstmLayer(numHiddenUnits2,'OutputMode','sequence')
    dropoutLayer(0.3)
    fullyConnectedLayer(numResponses)
    regressionLayer];

maxEpochs = 300;
miniBatchSize = 30;
validationFrequency = floor(Ntrain/miniBatchSize)*maxEpochs/100;

options = trainingOptions('adam', ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',false, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{X_validation,Y_validation}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','training-progress');

[trainedNet, traininfo] = trainNetwork(X_train,Y_train,layers, options);

%% train data
OUTPUT3 = predict(trainedNet,X_train);
for i = 1:numel(OUTPUT3)
   OUTPUT3{i} = OUTPUT3{i}.* sig2 + mu2;
end
Y_train_p = interpFT(OUTPUT3);

OUTPUT3_real = Y_train;
for i = 1:numel(OUTPUT3_real)
   OUTPUT3_real{i} = OUTPUT3_real{i}.* sig2 + mu2;
end
Y_train_real = interpFT(OUTPUT3_real);

fitness_tr = []; pos_ = pos(tr,:);
for j = 1:length(tr)
    g=(Po*Y_train_p{j}(1,1)-Cpw*Y_train_p{j}(2,1)-Ciw*Y_train_p{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_train_p{1},2)-1
        g=g+(Po*(Y_train_p{j}(1,k+1)-Y_train_p{j}(1,k)) ...
            -Cpw*(Y_train_p{j}(2,k+1)-Y_train_p{j}(2,k)) ...
            -Ciw*(Y_train_p{j}(3,k+1)-Y_train_p{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_tr = [fitness_tr;fit];
end

fitness_tr_real = [];
for j = 1:length(tr)
    g=(Po*Y_train_real{j}(1,1)-Cpw*Y_train_real{j}(2,1)-Ciw*Y_train_real{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_train_real{1},2)-1
        g=g+(Po*(Y_train_real{j}(1,k+1)-Y_train_real{j}(1,k)) ...
            -Cpw*(Y_train_real{j}(2,k+1)-Y_train_real{j}(2,k)) ...
            -Ciw*(Y_train_real{j}(3,k+1)-Y_train_real{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_tr_real = [fitness_tr_real;fit];
end
figure; h1=plotregression2(-fitness_tr_real, -fitness_tr);
% fig = gcf;
% fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' case_parameter '/training.png']);
close; 
%% validation data
OUTPUT = predict(trainedNet,X_validation);
for i = 1:numel(OUTPUT)
   OUTPUT{i} = OUTPUT{i}.* sig2 + mu2;
   OUTPUT{i}(OUTPUT{i} < 0) = 0;
end
Y_val_p = interpFT(OUTPUT);

OUTPUT_real = Y_validation;
for i = 1:numel(OUTPUT_real)
   OUTPUT_real{i} = OUTPUT_real{i}.* sig2 + mu2;
end
Y_val_real = interpFT(OUTPUT_real);

fitness_val = []; pos_ = pos(val,:);
for j = 1:length(val)
    g=(Po*Y_val_p{j}(1,1)-Cpw*Y_val_p{j}(2,1)-Ciw*Y_val_p{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_val_p{1},2)-1
        g=g+(Po*(Y_val_p{j}(1,k+1)-Y_val_p{j}(1,k)) ...
            -Cpw*(Y_val_p{j}(2,k+1)-Y_val_p{j}(2,k)) ...
            -Ciw*(Y_val_p{j}(3,k+1)-Y_val_p{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_val = [fitness_val;fit];
end

fitness_val_real = [];
for j = 1:length(val)
    g=(Po*Y_val_real{j}(1,1)-Cpw*Y_val_real{j}(2,1)-Ciw*Y_val_real{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_val_real{1},2)-1
        g=g+(Po*(Y_val_real{j}(1,k+1)-Y_val_real{j}(1,k)) ...
            -Cpw*(Y_val_real{j}(2,k+1)-Y_val_real{j}(2,k)) ...
            -Ciw*(Y_val_real{j}(3,k+1)-Y_val_real{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_val_real = [fitness_val_real;fit];
end
% fitness_val_real = fitness_m(val);
figure; h2=plotregression2(-fitness_val_real, -fitness_val);
% fig = gcf;
% fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' case_parameter '/validation.png']);
close; 
%% test data
bb = nonselected(1, randperm(Nall - Nsample, Ntest));
X_test = DATA(bb);
mu_ = mean([X_test{:}],2);
sig_ = std([X_test{:}],0,2);
mu_ = mu;
sig_ = sig;
 
for i = 1:numel(X_test)
    X_test{i} = (X_test{i} - mu_) ./ sig_;
end
Y_test = yDATA(bb);

OUTPUT2 = predict(trainedNet,X_test);
for i = 1:numel(OUTPUT2)
   OUTPUT2{i} = OUTPUT2{i}.* sig2 + mu2;
   OUTPUT2{i}(OUTPUT2{i} < 0) = 0;
end
Y_test_p = interpFT(OUTPUT2);

OUTPUT2_real = Y_test;
Y_test_real = interpFT(OUTPUT2_real);

fitness_test = []; pos_ = pos(bb,:);
for j = 1:length(bb)
    g=(Po*Y_test_p{j}(1,1)-Cpw*Y_test_p{j}(2,1)-Ciw*Y_test_p{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_test_p{1},2)-1
        g=g+(Po*(Y_test_p{j}(1,k+1)-Y_test_p{j}(1,k)) ...
            -Cpw*(Y_test_p{j}(2,k+1)-Y_test_p{j}(2,k)) ...
            -Ciw*(Y_test_p{j}(3,k+1)-Y_test_p{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_test = [fitness_test;fit];
end

fitness_test_real = [];
for j = 1:length(bb)
    g=(Po*Y_test_real{j}(1,1)-Cpw*Y_test_real{j}(2,1)-Ciw*Y_test_real{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_test_real{1},2)-1
        g=g+(Po*(Y_test_real{j}(1,k+1)-Y_test_real{j}(1,k)) ...
            -Cpw*(Y_test_real{j}(2,k+1)-Y_test_real{j}(2,k)) ...
            -Ciw*(Y_test_real{j}(3,k+1)-Y_test_real{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(pos_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    a = 1/3+(1-1/3)/4;
    for i = 1:N
        nw(i, ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_test_real = [fitness_test_real;fit];
end
figure; h3=plotregression2(-fitness_test_real, -fitness_test);
% fig = gcf;
% fig.PaperPosition = [5 5 8 8];
% ax = gca; ax.YTick = ax.XTick;
print('-r600', '-dpng', [dirname '/' case_parameter '/other samples.png']);
close;
% figure; h3=plotregression2(-fitness_m(bb), -fitness_test);

%% save
temp_mean{1} = mu; temp_mean{2} = mu2; temp_std{1} = sig; temp_std{2} = sig2;
pos_samp = pos(selected, :); pos_fit_samp = pos_fit;
RAW_samp{1} = RAW{1}(:,selected); RAW_samp{2} = RAW{2}(:,selected); RAW_samp{3} = RAW{3}(:,selected);
% saveas(h1, [dirname '/validation.png']); saveas(h1, [dirname '/validation.fig']);
% saveas(h2, [dirname '/training.png']); saveas(h2, [dirname '/training.fig']);
% saveas(h3, [dirname '/other samples.png']); saveas(h3, [dirname '/other samples.fig']);
save([dirname '/' case_parameter '/trainedNet.mat'], 'trainedNet', 'traininfo', 'temp_mean', 'temp_std', 'pos_samp', 'pos_fit_samp', 'RAW_samp', 'DATA', 'yDATA');
fid = fopen([dirname '/' case_parameter '/rmse ' num2str(traininfo.ValidationRMSE(end)) '.txt'], 'w'); fclose(fid); close all
close (findall(groot, 'Type', 'Figure'));
%% save data
save([dirname '/' case_parameter '/predicted.mat'], 'Y_test_real','Y_test_p','Y_train_p','Y_train_real','Y_val_p','Y_val_real',...
    'fitness_test_real', 'fitness_test','fitness_val_real', 'fitness_val','fitness_tr_real', 'fitness_tr');

rmpath([pwd '/code']);