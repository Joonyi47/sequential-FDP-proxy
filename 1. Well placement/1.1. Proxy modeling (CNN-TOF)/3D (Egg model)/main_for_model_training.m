addpath([pwd '/code']);
%% Initial conditions
case_sample = 1;
case_sensitivity = 3;

Nsample = Np*N_ens;
Ntrain = ceil(Nsample*0.7);
Ntest = Nsample * 0.15;
Ntrainle = Nsample - Ntest;  % traing + validation

selected = randperm(Nsample, Ntrainle);
selected = sort(selected, 'ascend');

pos_samp     = pos(selected, :);
pos_fit_samp = pos_fit(selected, :);
TOF_samp{1}  = total_TOF{1}(:,selected);
TOF_samp{2}  = total_TOF{2}(:,selected);
% P_samp{1}  = total_P{1}(:,selected);

maxtof = 30000; maxp = 400;
% total_TOF{1}(total_TOF{1} >= maxtof) = maxtof;
% total_TOF{2}(total_TOF{2} >= maxtof) = maxtof;
total_TOF_{1} = (total_TOF{1}-maxtof/2)/maxtof;
total_TOF_{2} = (total_TOF{2}-maxtof/2)/maxtof;
% total_P_{1} = (total_P{1}-maxp/2)/maxp;

varname = ['Result_' t(1:11) ' ' t(13:14) t(16:17)];   %%% sample data °á°ú º¯¼öÀÌ¸§À¸·Î ¼³Á¤
dirname = [varname '/' 'Ntrainle ' int2str(Ntrainle) '/sample case ' int2str(case_sample) '/sensitivity case ' int2str(case_sensitivity)];
mkdir(dirname);
%% Input data
case_parameter = '1';
Nlayer = 2;

INPUT = [];
for j = 1:Ntrainle

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
elseif case_sensitivity == 2
    INPUT(:,:,1,j) = (maxtof - reshape(total_TOF{1}(1:nx*ny,selected(j)),nx,ny)')/maxtof;
    INPUT(:,:,2,j) = (maxtof - reshape(total_TOF{2}(1:nx*ny,selected(j)),nx,ny)')/maxtof;
elseif case_sensitivity == 3 % TOF data
    for k = 1:Nlayer
        INPUT(:,:,k,j) = reshape(total_TOF_{1}(nx*ny*(k-1)+1:nx*ny*k,selected(j)),nx,ny)';
    end
    for k = 1:Nlayer
        INPUT(:,:,Nlayer+k,j) = reshape(total_TOF_{2}(nx*ny*(k-1)+1:nx*ny*k,selected(j)),nx,ny)';
    end
    INPUT(:,:,2*Nlayer+1,j) = 0*ones(ny,nx);
        
    tp = pos(:,2*N+1:3*N);

    for k = 1:N
        if tp(selected(j),k) >= 1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 2*Nlayer+1, j) = 1;
        elseif tp(selected(j),k) <= -1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 2*Nlayer+1, j) = -1;
        end
    end
elseif case_sensitivity == 4  % incorporate pressure data additioNsampley
    for k = 1:Nlayer
        INPUT(:,:,k,j) = reshape(total_TOF_{1}(nx*ny*(k-1)+1:nx*ny*k,selected(j)),nx,ny)';
    end
    for k = 1:Nlayer
        INPUT(:,:,Nlayer+k,j) = reshape(total_TOF_{2}(nx*ny*(k-1)+1:nx*ny*k,selected(j)),nx,ny)';
    end
    for k = 1:Nlayer
        INPUT(:,:,Nlayer*2+k,j) = reshape(total_P_{1}(nx*ny*(k-1)+1:nx*ny*k,selected(j)),nx,ny)';
    end
    
    INPUT(:,:,3*Nlayer+1,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(selected(j),k) >= 1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 3*Nlayer+1, j) = 1;
        elseif tp(selected(j),k) <= -1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), 3*Nlayer+1, j) = -1;
        end
    end
else % only pressure data
    for k = 1:Nlayer
        INPUT(:,:,k,j) = reshape(total_P_{1}(nx*ny*(k-1)+1:nx*ny*k,selected(j)),nx,ny)';
    end
    
    INPUT(:,:,Nlayer+1,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(selected(j),k) >= 1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), Nlayer+1, j) = 1;
        elseif tp(selected(j),k) <= -1/3
            INPUT(pos(selected(j),2*k), pos(selected(j),2*(k-1)+1), Nlayer+1, j) = -1;
        end
    end
end
end



temp_mean = mean(-pos_fit(selected,1));
temp_std = std(-pos_fit(selected,1), 1);
INPUT_NPV = ( - pos_fit(selected,1) - temp_mean ) / temp_std;

tr = randperm(Nsample,Ntrain);
val = [];
for iter = 1:Nsample
   if find(iter==tr)
      continue
   end
    val = [val,iter];
end

%% CNN model
miniBatchSize  = 80;
validationFrequency = floor((Nsample-Ntrain)/miniBatchSize);
options = trainingOptions('adam', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',50, ...
    'InitialLearnRate',3e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',20, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{INPUT(:,:,:,val),INPUT_NPV(val,:)}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','training-progress', ...
    'Verbose',false);
%     'Shuffle','every-epoch', ...

subdirname = ['BS ' int2str(miniBatchSize) ', MaxEpoch ' int2str(options.MaxEpochs) ', IniLR ' num2str(options.InitialLearnRate) ... 
                ', Ly ' num2str(Nlayer) '_' case_parameter]; 
home = cd(dirname);    
mkdir(subdirname);
cd(home)

layers = [
    imageInputLayer([ny nx size(INPUT,3)])

    convolution2dLayer(3,128,'Padding','same','WeightsInitializer','glorot')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2,'Stride',2)
  
    convolution2dLayer(3,128,'Padding','same','WeightsInitializer','glorot')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2,'Stride',2)

    convolution2dLayer(3,256,'Padding','same','WeightsInitializer','glorot')
    batchNormalizationLayer
    reluLayer
%     
    averagePooling2dLayer(2,'Stride',2)
%     
    convolution2dLayer(3,256,'Padding','same','WeightsInitializer','glorot')
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
        
    dropoutLayer(0.3)

    fullyConnectedLayer(1)
    batchNormalizationLayer    regressionLayer];


[trainedNet, traininfo] = trainNetwork(INPUT(:,:,:,tr), INPUT_NPV(tr,1),layers,options);

%% for training data and validation data
% modelfile = 'best_model.h5';
% trainedNet = importKerasNetwork(modelfile);
OUTPUT = predict(trainedNet,INPUT(:,:,:,val));
OUTPUT = [OUTPUT, INPUT_NPV(val)] .* temp_std + temp_mean;
figure; h1 = plotregression2(OUTPUT((OUTPUT(:,2)>0),2), OUTPUT((OUTPUT(:,2)>0),1)); ax = gca; ax.YLim = ax.XLim;
% p=polyfit(OUTPUT(:,1), OUTPUT(:,2),1);
fig = gcf;
fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' subdirname '/validation.png']); 
[~, I1] = sort(OUTPUT(:,1), 'ascend');
[~, I2] = sort(OUTPUT(:,2), 'ascend');
sorted_out = [I1, I2];
as = ismember(sorted_out(1:(Nsample - Ntrain)*0.1,1), sorted_out(1:(Nsample - Ntrain)*0.1,2));
asd1=length(find(as == 1))/((Nsample - Ntrain)*0.1)
as = ismember(sorted_out(1:30,2), sorted_out(1:30,1));
asd2=length(find(as == 1))/30

% 
OUTPUT2 = predict(trainedNet,INPUT(:,:,:,tr));
OUTPUT2 = [OUTPUT2, INPUT_NPV(tr)].* temp_std + temp_mean;
figure; h2 = plotregression2(OUTPUT2((OUTPUT2(:,2)>0),2), OUTPUT2((OUTPUT2(:,2)>0),1)); ax = gca; ax.YLim = ax.XLim;
fig = gcf;
fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' subdirname '/training.png']); 

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
    for k = 1:Nlayer
        INPUT(:,:,k,j) = reshape(total_TOF_{1}(nx*ny*(k-1)+1:nx*ny*k,nonselected(j)),nx,ny)';
    end
    for k = 1:Nlayer
        INPUT(:,:,Nlayer+k,j) = reshape(total_TOF_{2}(nx*ny*(k-1)+1:nx*ny*k,nonselected(j)),nx,ny)';
    end
    INPUT(:,:,2*Nlayer+1,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(nonselected(j),k) >= 1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 2*Nlayer+1, j) = 1;
        elseif tp(nonselected(j),k) <= -1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 2*Nlayer+1, j) = -1;
        end
    end
    
elseif case_sensitivity == 4    
    for k = 1:Nlayer
        INPUT(:,:,k,j) = reshape(total_TOF_{1}(nx*ny*(k-1)+1:nx*ny*k,nonselected(j)),nx,ny)';
    end
    for k = 1:Nlayer
        INPUT(:,:,Nlayer+k,j) = reshape(total_TOF_{2}(nx*ny*(k-1)+1:nx*ny*k,nonselected(j)),nx,ny)';
    end
    for k = 1:Nlayer
        INPUT(:,:,Nlayer*2+k,j) = reshape(total_P_{1}(nx*ny*(k-1)+1:nx*ny*k,nonselected(j)),nx,ny)';
    end
    
    INPUT(:,:,3*Nlayer+1,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(nonselected(j),k) >= 1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 3*Nlayer+1, j) = 1;
        elseif tp(nonselected(j),k) <= -1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), 3*Nlayer+1, j) = -1;
        end
    end
else
    for k = 1:Nlayer
        INPUT(:,:,k,j) = reshape(total_P_{1}(nx*ny*(k-1)+1:nx*ny*k,nonselected(j)),nx,ny)';
    end
    
    INPUT(:,:,Nlayer+1,j) = 0*ones(ny,nx);
    if isempty(type)
        tp = pos(:,2*N+1:3*N);
    else
        tp = repmat(type,N_ens*Np,1);
    end
    for k = 1:N
        if tp(nonselected(j),k) >= 1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), Nlayer+1, j) = 1;
        elseif tp(nonselected(j),k) <= -1/3
            INPUT(pos(nonselected(j),2*k), pos(nonselected(j),2*(k-1)+1), Nlayer+1, j) = -1;
        end
    end
end
end

temp_mean = mean(-pos_fit(nonselected,1));
temp_std = std(-pos_fit(nonselected,1), 1);
INPUT_NPV = ( - pos_fit(nonselected,1) - temp_mean ) / temp_std;

OUTPUT3 = predict(trainedNet, INPUT(:,:,:,:));
OUTPUT3 = [OUTPUT3, INPUT_NPV(:)].* temp_std + temp_mean;
figure; h3 = plotregression2(OUTPUT3((OUTPUT3(:,2)>0),2), OUTPUT3((OUTPUT3(:,2)>0),1)); ax = gca; ax.YLim = ax.XLim;
fig = gcf;
fig.PaperPosition = [5 5 8 8];
print('-r600', '-dpng', [dirname '/' subdirname '/test.png']); 
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
save([dirname '/' subdirname '/trainedNet.mat'], 'trainedNet', 'traininfo', 'temp_std', 'temp_mean', 'pos_samp', 'pos_fit_samp', 'TOF_samp', 'P_samp');
fid = fopen([dirname '/' subdirname '/rmse ' num2str(traininfo.ValidationRMSE(end)) '.txt'], 'w'); fclose(fid); close all
close (findall(groot, 'Type', 'Figure'));
%% save data
save([dirname '/predicted.mat'], 'OUTPUT', 'OUTPUT2', 'OUTPUT3', 'OUTPUT4');

rmpath([pwd '/code']);
