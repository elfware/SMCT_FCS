function [covIhExactMean,covIvExactMean,covIhDiscreteMean,covIvDiscreteMean,covTotExactMean,covTotDiscreteMean,covPExactMean,covPDiscreteMean,covCrossExactMean,covCrossDiscreteMean,TimeLag]=getAutoCorrelatio_correctnorm(ParamAutoCorr,I)

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

covTotExact=zeros(numTrajs,2*numlags+1);
covTotDiscrete=zeros(numTrajs,2*numlags+1);


covPExact=zeros(numTrajs,2*numlags+1);
covPDiscrete=zeros(numTrajs,2*numlags+1);

covCrossExact=zeros(numTrajs,2*numlags+1);
covCrossDiscrete=zeros(numTrajs,2*numlags+1);
for trajInd =1:numTrajs   
    covIhExact(trajInd,:)=xcov(IhExact(trajInd,:),numlags,'none');
    covIvExact(trajInd,:)=xcov(IvExact(trajInd,:),numlags,'none');
    covTotExact=xcov(IvExact(trajInd,:)+IhExact(trajInd,:),numlags,'none');
    PExact=(IhExact(trajInd,:)-IvExact(trajInd,:))./(IhExact(trajInd,:)+IvExact(trajInd,:));
    PExact(isnan(PExact))=0;
    covPExact(trajInd,:)=xcov(PExact,numlags,'none');
    covCrossExact(trajInd,:)=xcov(IhExact(trajInd,:),IvExact(trajInd,:),numlags,'none');
    if flag
    covIhDiscrete(trajInd,:)=xcov(IhDiscrete(trajInd,:),numlags,'none');
    covIvDiscrete(trajInd,:)=xcov(IvDiscrete(trajInd,:),numlags,'none');
    covTotDiscrete(trajInd,:)=xcov(IvDiscrete(trajInd,:)+IhDiscrete(trajInd,:),numlags,'none');
    PDiscrete=(IhDiscrete(trajInd,:)-IvDiscrete(trajInd,:))./(IhDiscrete(trajInd,:)+IvDiscrete(trajInd,:));
    PDiscrete(isnan(PDiscrete))=0;
    covPDiscrete(trajInd,:)=xcov(PDiscrete,numlags,'none');
    covCrossDiscrete(trajInd,:)=xcov(IhDiscrete(trajInd,:),IvDiscrete(trajInd,:),numlags,'none');
    end


end
Lag=[-numlags:numlags];


covIhExactMean=sum(covIhExact(:,numlags+1:end),1);
covIvExactMean=sum(covIvExact(:,numlags+1:end),1);
covTotExactMean=sum(covTotExact(:,numlags+1:end),1);
covPExactMean=sum(covPExact(:,numlags+1:end),1);
covCrossExactMean=sum(covCrossExact(:,numlags+1:end),1);

covIhExactMean=covIhExactMean/covIhExactMean(1);
covIvExactMean=covIvExactMean/covIvExactMean(1);
covTotExactMean=covTotExactMean/covTotExactMean(1);
covPExactMean=covPExactMean/covPExactMean(1);
covCrossExactMean=covCrossExactMean/covCrossExactMean(1);


if flag
  covIhDiscreteMean=sum(covIhDiscrete(:,numlags+1:end),1);
  covIvDiscreteMean=sum(covIvDiscrete(:,numlags+1:end),1);
  covTotDiscreteMean=sum(covTotDiscrete(:,numlags+1:end),1);
  covPDiscreteMean=sum(covPDiscrete(:,numlags+1:end),1);
  covCrossDiscreteMean=sum(covCrossDiscrete(:,numlags+1:end),1);
  
  covIhDiscreteMean=covIhDiscreteMean/covIhDiscreteMean(1);
  covIvDiscreteMean=covIvDiscreteMean/covIvDiscreteMean(1);
  covTotDiscreteMean=covTotDiscreteMean/covTotDiscreteMean(1);
  covPDiscreteMean=covPDiscreteMean/covPDiscreteMean(1);
  covCrossDiscreteMean=covCrossDiscreteMean/covCrossDiscreteMean(1);

else 
  covIhDiscreteMean=[];
  covIvDiscreteMean=[];
  covTotDiscreteMean=[];
  covPDiscreteMean=[];
  covCrossDiscreteMean=[];
end

TimeLag=dt*Lag(numlags+1:end);
%covIhExactMean(2:end)=ScaleAutocorr(1)*covIhExactMean(2:end);
%covIvExactMean(2:end)=ScaleAutocorr(2)*covIvExactMean(2:end);
