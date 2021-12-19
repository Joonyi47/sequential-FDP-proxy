function [XTrain,YTrain] = prepareDataTrain(filenamePredictors)

dataTrain = dlmread(filenamePredictors);

numObservations = max(dataTrain(:,1));

XTrain = cell(numObservations,1);
YTrain = cell(numObservations,1);
for i = 1:numObservations
    idx = dataTrain(:,1) == i;
    
    X = dataTrain(idx,3:end)';
    XTrain{i} = X;
    
    timeSteps = dataTrain(idx,2)';
    Y = fliplr(timeSteps);
    YTrain{i} = Y;
end

end