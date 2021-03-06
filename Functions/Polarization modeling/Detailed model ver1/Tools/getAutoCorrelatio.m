function [covIhExactMean covIvExactMean covIhDiscreteMean covIvDiscreteMean TimeLag]=getAutoCorrelatio(ParamAutoCorr,I)

IhExact=I.IhExact;
IvExact=I.IvExact;
IhDiscrete=I.IhDiscrete;
IvDiscrete=I.IvDiscrete;

numlags=ParamAutoCorr.numlags;
ScaleAutocorr=ParamAutoCorr.ScaleAutocorr;
dt=ParamAutoCorr.dt;
%%Autocorralation and plot
numTrajs=size(IhExact,1);
flag=not(isempty(IhDiscrete))&not(isempty(IvDiscrete));





covIvExact=zeros(numTrajs,2*numlags+1);
covIhExact=zeros(numTrajs,2*numlags+1);
covIvDiscrete=zeros(numTrajs,2*numlags+1);
covIhDiscrete=zeros(numTrajs,2*numlags+1);
for trajInd =1:numTrajs   
    covIhExact(trajInd,:)=xcov(IhExact(trajInd,:),numlags,'none');
    covIvExact(trajInd,:)=xcov(IvExact(trajInd,:),numlags,'none');
    if flag
    covIhDiscrete(trajInd,:)=xcov(IhDiscrete(trajInd,:),numlags,'none');
    covIvDiscrete(trajInd,:)=xcov(IvDiscrete(trajInd,:),numlags,'none');
    end
    covIhExact(trajInd,:)=covIhExact(trajInd,:)/covIhExact(trajInd,numlags+1);
    covIvExact(trajInd,:)=covIvExact(trajInd,:)/covIvExact(trajInd,numlags+1);
    covIhDiscrete(trajInd,:)=covIvDiscrete(trajInd,:)/covIhDiscrete(trajInd,numlags+1);
    covIvDiscrete(trajInd,:)=covIvDiscrete(trajInd,:)/covIvDiscrete(trajInd,numlags+1);
end
Lag=[-numlags:numlags];


covIhExactMean=mean(covIhExact(:,numlags+1:end),1);
covIvExactMean=mean(covIvExact(:,numlags+1:end),1);

if flag
  covIhDiscreteMean=mean(covIhDiscrete(:,numlags+1:end),1);
  covIvDiscreteMean=mean(covIvDiscrete(:,numlags+1:end),1);
else 
  covIhDiscreteMean=[];
  covIvDiscreteMean=[];
end

TimeLag=dt*Lag(numlags+1:end);
covIhExactMean(2:end)=ScaleAutocorr(1)*covIhExactMean(2:end);
covIvExactMean(2:end)=ScaleAutocorr(2)*covIvExactMean(2:end);
