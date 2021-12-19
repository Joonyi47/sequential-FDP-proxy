addpath([pwd '/code']);
%% Initial conditions
case_sample = 1;
case_sensitivity = 3;

Nall = 5000;
Ntrain = ceil(Nall*0.7);  % 70% used for training
Ntest = Nall * 0.15;      % 15% used for test
Nsample = Nall - Ntest;  % traing + validation

selected = randperm(Nall, Nsample);
selected = sort(selected, 'ascend');

pos_samp     = pos(selected, :);
pos_fit_samp = pos_fit(selected, :);
TOF_samp{1}  = total_TOF{1}(:,selected);
TOF_samp{2}  = total_TOF{2}(:,selected);
% P_samp{1}  = total_P{1}(:,selected);

maxtof = 10000; maxp = 3500;
% total_TOF_{1} = (total_TOF{1} - mean(mean(total_TOF{1}))) ./ mean(std(total_TOF{1},0,1));
% total_TOF_{2} = (total_TOF{2} - mean(mean(total_TOF{2}))) ./ mean(std(total_TOF{2},0,1));
total_TOF_{1} = (total_TOF{1}-maxtof/2)/maxtof;
total_TOF_{2} = (total_TOF{2}-maxtof/2)/maxtof;
% total_P_{1} = (total_P{1}-maxp)/maxp;

varname = ['Result_' t(1:11) ' ' t(13:14) t(16:17)];   %%% sample data 결과 변수이름으로 설정
dirname = [varname '/' 'Nsample ' int2str(Nsample) '/sample case ' int2str(case_sample) '/sensitivity case ' int2str(case_sensitivity)];
mkdir(dirname);
%% Input data
case_parameter = '1';

INPUT = [];
for j = 1:Nsample

if case_sensitivity == 1  
    INPUT(:,:,1,j) = (maxtof-reshape((total_TOF{1}(1:nx*ny,selected(j))+total_TOF{2}(1:nx*ny,selected(j))),nx,ny)')/maxtof;
    INPUT(:,:,2,j) = (maxtof-reshape((total_TOF{1}(nx*ny+1:2*nx*ny,selected(j))+total_TOF{2}(nx*ny+1:2*nx*ny,selected(j))),nx,ny)')/maxtof;
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    INPUT(:,:,3,j) = 0  *ones(ny,nx);
    for k = 1:N
        if tp(selected(j),k) == 1
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 3, j) = 1;
        else
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 3, j) = -1;
        end
    end
elseif case_sensitivity == 2 % only TOF maps
    INPUT(:,:,1,j) = (maxtof - reshape(total_TOF{1}(1:nx*ny,selected(j)),nx,ny)')/maxtof;
    INPUT(:,:,2,j) = (maxtof - reshape(total_TOF{2}(1:nx*ny,selected(j)),nx,ny)')/maxtof;
elseif case_sensitivity == 3 % TOF maps + well configuration map
    INPUT(:,:,1,j) = reshape(total_TOF_{1}(1:nx*ny,selected(j)),nx,ny)';
    INPUT(:,:,2,j) = reshape(total_TOF_{2}(1:nx*ny,selected(j)),nx,ny)';
    INPUT(:,:,3,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(selected(j),k) >= 1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 3, j) = 1;
        elseif tp(selected(j),k) <= -1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 3, j) = -1;
        end             
    end
else % TOF maps + Pressure map + well configuration map
    INPUT(:,:,1,j) = reshape(total_TOF_{1}(1:nx*ny,selected(j)),nx,ny)';
    INPUT(:,:,2,j) = reshape(total_TOF_{2}(1:nx*ny,selected(j)),nx,ny)';
    INPUT(:,:,3,j) = reshape(total_P_{1}(1:nx*ny,selected(j)),nx,ny)';
    INPUT(:,:,4,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(selected(j),k) >= 1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 4, j) = 1;
        elseif tp(selected(j),k) <= -1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 4, j) = -1;
        end             
    end
end
end


temp_mean = mean(-pos_fit(selected,1));
temp_std = std(-pos_fit(selected,1), 1);
INPUT_NPV = ( - pos_fit(selected,1) - temp_mean ) / temp_std;

tr = randperm(Nsample,Ntrain);
val = [];
for iter = 1:1:Nsample
   if find(iter==tr)
      continue
   end
    val = [val,iter];
end


%% CNN model
miniBatchSize  = 120;
validationFrequency = floor((Nsample-Ntrain)/miniBatchSize);
options = trainingOptions('adam', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',60, ...
    'InitialLearnRate',3e-3, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{INPUT(:,:,:,val),INPUT_NPV(val,:)}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','training-progress', ...
    'Verbose',false, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',20);
%     'Shuffle','every-epoch', ...

subdirname = ['BS ' int2str(miniBatchSize) ', ME ' int2str(options.MaxEpochs) ', ILR ' num2str(options.InitialLearnRate) '_' case_parameter]; 
home = cd(dirname);    
mkdir(subdirname);
cd(home)

layers = [
    imageInputLayer([ny nx size(INPUT,3)])

    convolution2dLayer(3,32,'Padding','same','WeightsInitializer','glorot')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2,'Stride',2)
  
    convolution2dLayer(3,64,'Padding','same','WeightsInitializer','glorot')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2,'Stride',2)

    convolution2dLayer(3,64,'Padding','same','WeightsInitializer','glorot')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,128,'Padding','same','WeightsInitializer','glorot')
    batchNormalizationLayer
    reluLayer

    averagePooling2dLayer(2,'Stride',2)

    fullyConnectedLayer(128)
    batchNormalizationLayer
    reluLayer
        
    dropoutLayer(0.4)

    fullyConnectedLayer(1)
    batchNormalizationLayer

    regressionLayer];


[trainedNet, traininfo] = trainNetwork(INPUT(:,:,:,tr), INPUT_NPV(tr,1),layers,options);

%% for training data and validation data
% modelfile = 'best_model.h5';
% trainedNet = importKerasNetwork(modelfile);
OUTPUT = predict(trainedNet,INPUT(:,:,:,val));
OUTPUT = [OUTPUT, INPUT_NPV(val)] .* temp_std + temp_mean;
figure; h1 = plotregression2(OUTPUT(:,2), OUTPUT(:,1)); ax = gca; ax.YLim = ax.XLim;
fig = gcf;
fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' subdirname '/validation.png']); 
p=polyfit(OUTPUT(:,1), OUTPUT(:,2),1); p(1)
[~, I1] = sort(OUTPUT(:,1), 'ascend');
[~, I2] = sort(OUTPUT(:,2), 'ascend');
sorted_out = [I1, I2];
as = ismember(sorted_out(1:(Nsample - Ntrain)*0.1,1), sorted_out(1:(Nsample - Ntrain)*0.1,2));
asd1=length(find(as == 1))/((Nsample - Ntrain)*0.1);
as = ismember(sorted_out(1:30,2), sorted_out(1:30,1));
asd2=length(find(as == 1))/30

% 
OUTPUT2 = predict(trainedNet,INPUT(:,:,:,tr));
OUTPUT2 = [OUTPUT2, INPUT_NPV(tr)].* temp_std + temp_mean;
figure; h2 = plotregression2(OUTPUT2(:,2), OUTPUT2(:,1)); ax = gca; ax.YLim = ax.XLim;
fig = gcf;
fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' subdirname '/training.png']); 
p=polyfit(OUTPUT2(:,1), OUTPUT2(:,2),1); p(1)

%% for test data
nonselected = [];
for j = 1:Nall
   if ~ismember(j, selected)
       nonselected = [nonselected, j];
   end
end

INPUT = [];
for j = 1:Ntest
if case_sensitivity == 1
    INPUT(:,:,1,j) = (maxtof-reshape((total_TOF{1}(1:nx*ny,nonselected(j))+total_TOF{2}(1:nx*ny,nonselected(j))),nx,ny)')/maxtof;
    INPUT(:,:,2,j) = (maxtof-reshape((total_TOF{1}(nx*ny+1:2*nx*ny,nonselected(j))+total_TOF{2}(nx*ny+1:2*nx*ny,nonselected(j))),nx,ny)')/maxtof;
    tp = pos(:,2*N+1:3*N);
    INPUT(:,:,3,j) = 0*ones(ny,nx);
    for k = 1:N
        if tp(nonselected(j),k) == 1
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 3, j) = 1;
        else
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 3, j) = -1;
        end
    end
elseif case_sensitivity == 2
    INPUT(:,:,1,j) = (maxtof - reshape(total_TOF{1}(:,nonselected(j)),nx,ny)')/maxtof;
    INPUT(:,:,2,j) = (maxtof - reshape(total_TOF{2}(:,nonselected(j)),nx,ny)')/maxtof;
elseif case_sensitivity == 3
    INPUT(:,:,1,j) = reshape(total_TOF_{1}(1:nx*ny,nonselected(j)),nx,ny)';
    INPUT(:,:,2,j) = reshape(total_TOF_{2}(1:nx*ny,nonselected(j)),nx,ny)';
    INPUT(:,:,3,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(nonselected(j),k) >= 1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 3, j) = 1;
        elseif tp(nonselected(j),k) <= -1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 3, j) = -1;
        end
    end
else
    INPUT(:,:,1,j) = reshape(total_TOF_{1}(1:nx*ny,nonselected(j)),nx,ny)';
    INPUT(:,:,2,j) = reshape(total_TOF_{2}(1:nx*ny,nonselected(j)),nx,ny)';
    INPUT(:,:,3,j) = reshape(total_P_{1}(1:nx*ny,nonselected(j)),nx,ny)';
    INPUT(:,:,4,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(nonselected(j),k) >= 1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 4, j) = 1;
        elseif tp(nonselected(j),k) <= -1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 4, j) = -1;
        end             
    end
end
  
end

temp_mean = mean(-pos_fit(nonselected,1));
temp_std = std(-pos_fit(nonselected,1), 1);
INPUT_NPV = ( - pos_fit(nonselected,1) - temp_mean ) / temp_std;

OUTPUT3 = predict(trainedNet, INPUT(:,:,:,:));
OUTPUT3 = [OUTPUT3, INPUT_NPV(:)].* temp_std + temp_mean;
figure; h3 = plotregression2(OUTPUT3(:,2), OUTPUT3(:,1)); ax = gca; ax.YLim = ax.XLim;
fig = gcf;
fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' subdirname '/test.png']); 
p=polyfit(OUTPUT3(:,1), OUTPUT3(:,2),1); p(1)
pos_fit_test = pos_fit(nonselected,1);
[~, I3] = sort(-pos_fit_test, 'descend');
OUTPUT4 = predict(trainedNet, INPUT(:,:,:,I3(1:40,1)));
OUTPUT4 = [OUTPUT4, INPUT_NPV(I3(1:40,1))].* temp_std + temp_mean;
figure; h4 = plotregression2(OUTPUT4(:,2), OUTPUT4(:,1));
[~, I1] = sort(OUTPUT3(:,1), 'ascend');
[~, I2] = sort(OUTPUT3(:,2), 'ascend');
sorted_out = [I1, I2];
as = ismember(sorted_out(1:(Nsample - Ntrain)*0.1,1), sorted_out(1:(Nsample - Ntrain)*0.1,2));
asd3=length(find(as == 1))/((Nsample - Ntrain)*0.1)
as = ismember(sorted_out(1:30,2), sorted_out(1:30,1));
asd4=length(find(as == 1))/30

%% save
saveas(h1, [dirname '/' subdirname '/validation' '(' num2str(asd1) ',' num2str(asd2) ')' '.fig']);
saveas(h2, [dirname '/' subdirname '/training.fig']);
saveas(h3, [dirname '/' subdirname '/other samples' '(' num2str(asd3) ',' num2str(asd4) ')' '.fig']);
% save([dirname '/' subdirname '/trainedNet.mat'], 'trainedNet', 'traininfo', 'temp_std', 'temp_mean', 'pos_samp', 'pos_fit_samp', 'TOF_samp', 'P_samp');
save([dirname '/' subdirname '/trainedNet.mat'], 'trainedNet', 'traininfo', 'temp_std', 'temp_mean', 'pos_samp', 'pos_fit_samp', 'TOF_samp');
fid = fopen([dirname '/' subdirname '/rmse ' num2str(traininfo.ValidationRMSE(end)) '.txt'], 'w'); fclose(fid); close all
close (findall(groot, 'Type', 'Figure'));
%% save data
save([dirname '/predicted.mat'], 'OUTPUT', 'OUTPUT2', 'OUTPUT3', 'OUTPUT4');

rmpath([pwd '/code']);