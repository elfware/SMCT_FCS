function [covIhExact,covIvExact,covIhDiscrete,covIvDiscrete,covTotExact,covTotDiscrete,covPExact,covPDiscrete,covCrossExact,covCrossDiscrete,TimeLag]=getAutoCorrelatio_correctnorm_newinformat_single(ParamAutoCorr,IhExact,IvExact,IhDiscrete,IvDiscrete)

%IhExact=I.IhExact;
%IvExact=I.IvExact;
%IhDiscrete=I.IhDiscrete;
%IvDiscrete=I.IvDiscrete;

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
    covTotExact(trajInd,:)=xcov(IvExact(trajInd,:)+IhExact(trajInd,:),numlags,'none');
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


%covIhExactMean=sum(covIhExact(:,numlags+1:end),1);
%covIvExactMean=sum(covIvExact(:,numlags+1:end),1);
%covTotExactMean=sum(covTotExact(:,numlags+1:end),1);
%covPExactMean=sum(covPExact(:,numlags+1:end),1);
%covCrossExactMean=sum(covCrossExact(:,numlags+1:end),1);

covIhExact=covIhExact(:,numlags+1:end);
covIvExact=covIvExact(:,numlags+1:end);
covTotExact=covTotExact(:,numlags+1:end);
covPExact=covPExact(:,numlags+1:end);
covCrossExact=covCrossExact(:,numlags+1:end);

covIhExact=covIhExact./covIhExact(:,1)*ones(1,numel(numlags));
covIvExact=covIvExact./covIvExact(:,1)*ones(1,numel(numlags));
covTotExact=covTotExact./covTotExact(:,1)*ones(1,numel(numlags));
covPExact=covPExact./covPExact(:,1)*ones(1,numel(numlags));
covCrossExact=covCrossExact./covCrossExact(:,1)*ones(1,numel(numlags));


if flag
  %covIhDiscreteMean=sum(covIhDiscrete(:,numlags+1:end),1);
  %covIvDiscreteMean=sum(covIvDiscrete(:,numlags+1:end),1);
  %covTotDiscreteMean=sum(covTotDiscrete(:,numlags+1:end),1);
  %covPDiscreteMean=sum(covPDiscrete(:,numlags+1:end),1);
  %covCrossDiscreteMean=sum(covCrossDiscrete(:,numlags+1:end),1);
  
  covIhDiscrete=covIhDiscrete(:,numlags+1:end);
  covIvDiscrete=covIvDiscrete(:,numlags+1:end);
  covTotDiscrete=covTotDiscrete(:,numlags+1:end);
  covPDiscrete=covPDiscrete(:,numlags+1:end);
  covCrossDiscrete=covCrossDiscrete(:,numlags+1:end);
  
  
  covIhDiscrete=covIhDiscrete./covIhDiscrete(:,1)*ones(1,numel(numlags));
  covIvDiscrete=covIvDiscrete./covIvDiscrete(:,1)*ones(1,numel(numlags));
  covTotDiscrete=covTotDiscrete./covTotDiscrete(:,1)*ones(1,numel(numlags));
  covPDiscrete=covPDiscrete./covPDiscrete(:,1)*ones(1,numel(numlags));
  covCrossDiscrete=covCrossDiscrete./covCrossDiscrete(:,1)*ones(1,numel(numlags));

else 
  covIhDiscrete=[];
  covIvDiscrete=[];
  covTotDiscrete=[];
  covPDiscrete=[];
  covCrossDiscrete=[];
end

TimeLag=dt*Lag(numlags+1:end);
%covIhExactMean(2:end)=ScaleAutocorr(1)*covIhExactMean(2:end);
%covIvExactMean(2:end)=ScaleAutocorr(2)*covIvExactMean(2:end);
