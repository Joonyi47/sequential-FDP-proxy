function [trainedNet, traininfo, temp_mean, temp_std] = trainRNN(total_RAW_m, pos, type, wset, pos_fit)
global N Np ...
       opt_type ...
       pmax Nstep dstep tstep ...
       discount_rate observed_term discount_term...
       Cw Po Cpw Ciw ...
       direc_fig retrain iter
       
fun = @(x,y,z,b,t)x*(1+cos(y*(1-t)+z))+b;

type_ = zeros(size(type,1), N*Nstep);
for i = 1:size(type,1)
    tmp = type(i,:);
    tmp = reshape(tmp, 4, N);
    schd = zeros(N, Nstep);
    for k = 1:N
        schd(k, (tmp(1,k) == 0)+ceil(Nstep * tmp(1,k)):end) = 1;
    end
    
    tmp3 = [];
   for j = 1:N
       if opt_type(j) >= 1/3
           tmp2 = feval(fun, (wset(2) - wset(1))/2*tmp(2,j), tmp(3,j), tmp(4,j), wset(1), [1/Nstep:1/Nstep:1]);
           tmp2 = tmp2.*schd(j,:);
%            tmp2(tmp2 == 0) = wset(2);
           tmp2(tmp2 == 0) = 3500;
       else
           tmp2 = feval(fun, (wset(4) - wset(3))/2*tmp(2,j), tmp(3,j), tmp(4,j), wset(3), [1/Nstep:1/Nstep:1]);
           tmp2 = tmp2.*schd(j,:);
%            tmp2(tmp2 == 0) = wset(3);
           tmp2(tmp2 == 0) = 3500;
       end
       tmp3 = [tmp3; tmp2];
   end
   
   type_(i,:) = reshape(tmp3', 1, numel(tmp3));
   
   
end

Nsample = size(type,1);
Ntrain  = ceil(Nsample*0.8);
selected = 1:Nsample;

aa = randperm(Nsample,Ntrain);
tr = selected(aa);
val = [];
for s = selected
   if find(s==tr)
      continue
   end
   val = [val,s];
end

idx  = opt_type;
idx  = (idx >= 1/3);   %% prod == 1, inj == 0

%%% transform pos data of n-th gen. for rnn network %%%
%%%% -- input -- %%%%
data = cell(Nsample, 1);
for i = 1:Nsample
    tmp = reshape(type_(i,:), Nstep, N)';
    data{i} = [tmp(idx == 1,:); tmp(idx == 0,:)];  %% sort (prod -> inj
end
%%%% -- output -- %%%%
ydata = cell(Nsample, 1);
for i = 1:Nsample
    ydata{i} = [total_RAW_m{1}(pmax/tstep/Nstep:pmax/tstep/Nstep:end,i)'; ...
                total_RAW_m{2}(pmax/tstep/Nstep:pmax/tstep/Nstep:end,i)'; ...
                total_RAW_m{3}(pmax/tstep/Nstep:pmax/tstep/Nstep:end,i)'];
end


%%% nomalize the data %%%
X_  = data;
temp_mean{1} = mean([X_{:}],2);
temp_std{1} = std([X_{:}],0,2);

for i = 1:numel(ydata)
   Y_(:,:,i) = ydata{i};
end
temp_mean{2} = mean(Y_,3);
temp_std{2} = std(Y_,0,3);

X_train = data(tr);
X_validation = data(val);
for i = 1:numel(X_train)
    X_train{i} = (X_train{i} - temp_mean{1}) ./ temp_std{1};
end
for i = 1:numel(X_validation)
    X_validation{i} = (X_validation{i} - temp_mean{1}) ./ temp_std{1};
end

Y_train = ydata(tr);
Y_validation = ydata(val);
for i = 1:numel(Y_train)
    Y_train{i} = (Y_train{i} - temp_mean{2}) ./ temp_std{2};
end
for i = 1:numel(Y_validation)
    Y_validation{i} = (Y_validation{i} - temp_mean{2}) ./ temp_std{2};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFeatures = N;
numHiddenUnits1 = 125;
numHiddenUnits2 = 100;
numHiddenUnits3 = 75;
numResponses = size(ydata{1},1);
layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits1,'OutputMode','sequence')
    dropoutLayer(0.3)
    lstmLayer(numHiddenUnits2,'OutputMode','sequence')
    dropoutLayer(0.3)
%     lstmLayer(numHiddenUnits3,'OutputMode','sequence')
%     dropoutLayer(0.2)
    fullyConnectedLayer(numResponses)
    regressionLayer];

maxEpochs = 300;
miniBatchSize = 30;

options = trainingOptions('adam', ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',false, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{X_validation,Y_validation}, ...
    'Plots','training-progress');

[trainedNet, traininfo] = trainNetwork(X_train,Y_train,layers, options);

%% visualization
%%%% train data %%%%
OUTPUT3 = predict(trainedNet,X_train);
for i = 1:numel(OUTPUT3)
   OUTPUT3{i} = OUTPUT3{i}.* temp_std{2} + temp_mean{2};
end
Y_train_p = interpFT(OUTPUT3);

OUTPUT3_real = Y_train;
for i = 1:numel(OUTPUT3_real)
   OUTPUT3_real{i} = OUTPUT3_real{i}.* temp_std{2} + temp_mean{2};
end
Y_train_real = interpFT(OUTPUT3_real);

fitness_tr = []; type_ = type(tr,:);
for j = 1:length(tr)
    g=(Po*Y_train_p{j}(1,1)-Cpw*Y_train_p{j}(2,1)-Ciw*Y_train_p{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_train_p{1},2)-1
        g=g+(Po*(Y_train_p{j}(1,k+1)-Y_train_p{j}(1,k)) ...
            -Cpw*(Y_train_p{j}(2,k+1)-Y_train_p{j}(2,k)) ...
            -Ciw*(Y_train_p{j}(3,k+1)-Y_train_p{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(type_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    for i = 1:N
        nw(i, (tp(1,i) == 0) + ceil((Nstep) * tp(1,i))) = 1;
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
    tp = reshape(type_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    for i = 1:N
        nw(i, (tp(1,i) == 0) + ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_tr_real = [fitness_tr_real;fit];
end
% fitness_val_real = fitness_m(val);
figure; h1=plotregression2(-fitness_tr_real, -fitness_tr);
print('-r600', '-dpng', [direc_fig '/training data #' num2str(retrain) '_' num2str(iter) '.png']); close;

%% validation
OUTPUT = predict(trainedNet,X_validation);
for i = 1:numel(OUTPUT)
   OUTPUT{i} = OUTPUT{i}.* temp_std{2} + temp_mean{2};
end
Y_val_p = interpFT(OUTPUT);

OUTPUT_real = Y_validation;
for i = 1:numel(OUTPUT_real)
   OUTPUT_real{i} = OUTPUT_real{i}.* temp_std{2} + temp_mean{2};
end
Y_val_real = interpFT(OUTPUT_real);

fitness_val = []; type_ = type(val,:);
for j = 1:length(val)
    g=(Po*Y_val_p{j}(1,1)-Cpw*Y_val_p{j}(2,1)-Ciw*Y_val_p{j}(3,1))/((1+discount_rate)^(observed_term/discount_term));
    for k=1:size(Y_val_p{1},2)-1
        g=g+(Po*(Y_val_p{j}(1,k+1)-Y_val_p{j}(1,k)) ...
            -Cpw*(Y_val_p{j}(2,k+1)-Y_val_p{j}(2,k)) ...
            -Ciw*(Y_val_p{j}(3,k+1)-Y_val_p{j}(3,k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
    end
    tp = reshape(type_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    for i = 1:N
        nw(i, (tp(1,i) == 0) + ceil((Nstep) * tp(1,i))) = 1;
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
    tp = reshape(type_(j,:), 4, N);
    nw = zeros(N ,Nstep);
    for i = 1:N
        nw(i, (tp(1,i) == 0) + ceil((Nstep) * tp(1,i))) = 1;
    end
    nw = sum(nw,1);
    for k =1:Nstep
        g=g-Cw*nw(k)/((1+discount_rate)^((pmax/Nstep*(k-1)/discount_term)));
    end
    [fit] = -g;
    fitness_val_real = [fitness_val_real;fit];
end
figure; h2=plotregression2(-fitness_val_real, -fitness_val);
print('-r600', '-dpng', [direc_fig '/validation data #' num2str(retrain) '_' num2str(iter) '.png']); close;

end