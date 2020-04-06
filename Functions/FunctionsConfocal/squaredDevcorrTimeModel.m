function [squaredDev]=squaredDevcorrTimeModel(params,k1,minTimes,maxTimes,Dxs,corrTimeMeas)
    numTimeScales=numel(minTimes);
    squaredDev=0;
    
    for i = 1:numTimeScales
        modelvals=corrTimeModel([params(1) params(i*2) params(i*2+1)],k1,minTimes(i),maxTimes(i),Dxs{i});
        numDatapoints=numel(Dxs{i});
        squaredDev=squaredDev+sum((abs((modelvals-corrTimeMeas{i})))/(numDatapoints*(maxTimes(i)-minTimes(i))));
    end
end