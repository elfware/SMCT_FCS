function [vec]=squaredDevcorrTimeModelvec(params,k1,minTimes,maxTimes,Dxs,corrTimeMeas)
    numTimeScales=numel(minTimes);
    squaredDev=0;
    vec=zeros(1,numTimeScales);
    for i = 1:numTimeScales
        modelvals=corrTimeModel([params(1) params(i*2) params(i*2+1)],k1,minTimes(i),maxTimes(i),Dxs{i});
        numDatapoints=numel(Dxs{i});
        if (params(1)<0||params(i*2)<0||params(i*2)>params(i*2+1)*(maxTimes(i)-minTimes(i))||params(i*2+1)<0)
            vec(i)=10^5;
        else    
            vec(i)=sum((abs((modelvals-corrTimeMeas{i})))/(numDatapoints*std(corrTimeMeas{i})));
        end
    end
end