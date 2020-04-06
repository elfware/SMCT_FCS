% Main function for fitting difference i autocorrelation
function [expFitCombDiff] = fitAutoCorrModelv2(timeWind,normNumber,kguessExpMax,kguessExpMin,numGuesses,theScales,compare,lb,ub,numBoots,numPreBoots,perc,scaleRes)
%UNTITLED4 Summary of this function goes here
% v2 Fixed so that sum squared deviation is minimzed

%theScales=[5:26];
numTerms=1;
numScales=numel(theScales);




%expsToFit=[2 3 4 5 6 7];
%bgToUse=[12 13 14 15 16 17];
%numFits=numel(expsToFit);






numFits=numel(compare);
bestparams=cell(1,numFits);
bestfval=zeros(1,numFits);
expFitCombDiff=cell(1,numFits);
totNumParams=2*numTerms;
numE=1;
startguess=zeros(1,totNumParams+numScales);
wrong=0;


for i=1:numFits
    disp(['Current Fit :' num2str(i)])
    theExp=i;
    lagTimes=zeros(numScales,timeWind);
    autoCorrsDiffBootPre=zeros(numScales,timeWind,numPreBoots);
    % Do pre bootstrapping to estimate the variance/ weights of the
    % datapoints which are use to weight the fitting
    parfor boot=1:numPreBoots
        autoCorrsDiff=zeros(numScales,timeWind);
        
        for j=1:numScales
            [numtrajs1,~]=size(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok);
            [numtrajs2,~]=size(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok);
            bootInds1=ceil(rand(1,numtrajs1)*numtrajs1);
            bootInds2=ceil(rand(1,numtrajs2)*numtrajs2);
            autoCorrsDiff(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok(bootInds1,:),perc)-trimmean(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok(bootInds2,:),perc);
          %  lagTimes(j,:)=scaleRes{theScales(j)}.lagTimes;
            
        end
        autoCorrsDiffBootPre(:,:,boot)=autoCorrsDiff;
    end
    
    
    stdDiffs=std(autoCorrsDiffBootPre,[],3);
    weights=1./stdDiffs;
    weights(weights==inf)=0;
    weights=weights';
    weights=weights(:);
    weights=weights';
    autoCorrsDiff=zeros(numScales,timeWind);
    bgAmpl=zeros(1,numScales);
    for j=1:numScales
        autoCorrsDiff(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok,perc)-trimmean(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok,perc);
        lagTimes(j,:)=scaleRes{theScales(j)}.lagTimes;
        bgAmpl(j)=1;%combUnNorm{i+3}.autoCorrScaled(j,1)-combUnNorm{i+3}.autoCorrScaled(j,end);
        
    end
    
    %bgAmpl=bgAmpl/sum(bgAmpl);
  
    %Function to minimize definition
    funDiff = @(x) (weights.*diffAutoCorrDiffMultScales3(ones(1,numTerms),x(totNumParams/2+1:totNumParams),x(end)*bgAmpl,x(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber)).^2;
    
    
    opt=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',10000,'Display','off');
    opt.MaxFunEvals=1000;
    %,'PlotFcns',@optimplotfval,'Display','iter');
    bestparams{i}=zeros(1,totNumParams+1);
    bestfval(i)=inf;
    startguess=zeros(1,totNumParams+1);
    wrong=0;
    
    for theG=1:numGuesses
        
        startguess(1:totNumParams/2)=1*abs(randn(1,totNumParams/2)); %kfactor
        startguess(totNumParams/2+1:totNumParams)=-10.^(kguessExpMin+rand(1,totNumParams/2)*(kguessExpMax-kguessExpMin)); %k
        %startguess(totNumParams+1:totNumParams+numScales)=10*abs(randn(1,numScales)); %bextra
        startguess(end)=10*abs(randn(1,1));
        try
        [params,fval]=lsqnonlin(funDiff,startguess,lb,ub,opt);
        if bestfval(i)>fval
            bestfval(i)=fval;
            bestparams{i}=params;
        end
        catch
        end
        
   end
    
    bgRes=diffAutoCorrDiffMultScales3BgRes(ones(1,numTerms),bestparams{i}(totNumParams/2+1:totNumParams),bestparams{i}(end)*bgAmpl,bestparams{i}(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber);
    
    
    % expFitRes{i}.aVals=bestparams{i}(1);
    % expFitRes{i}.kVals=bestparams{i}(2);
    expFitCombDiff{i}.bestparams=bestparams{i};
    expFitCombDiff{i}.bestfval=bestfval(i);
    
    expFitCombDiff{i}.autoCorrsDiff=autoCorrsDiff;
    expFitCombDiff{i}.lagTimes=lagTimes;
    expFitCombDiff{i}.bgAmpl=bgAmpl;
    %expFitCombDiff{i}.stdDiffs=stdDiff;
    expFitCombDiff{i}.weights=weights;
    expFitCombDiff{i}.bgRes=bgRes;
    %%Bootstrapping
    
    bestparamsBoot=zeros(numBoots,totNumParams+1);
    bestfvalBoot=zeros(1,numBoots);
    autoCorrsDiffBoot=zeros(numScales,timeWind,numBoots);
    bgResBoot=zeros(numScales,numBoots);
    parfor boot=1:numBoots
        autoCorrsDiff=zeros(numScales,timeWind);
        
        for j=1:numScales
            [numtrajs1,~]=size(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok);
            [numtrajs2,~]=size(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok);
            bootInds1=ceil(rand(1,numtrajs1)*numtrajs1);
            bootInds2=ceil(rand(1,numtrajs2)*numtrajs2);
            autoCorrsDiff(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok(bootInds1,:),perc)-trimmean(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok(bootInds2,:),perc);
          %  lagTimes(j,:)=scaleRes{theScales(j)}.lagTimes;
            
        end
        
        
        funDiff = @(x) (weights.*diffAutoCorrDiffMultScales3(ones(1,numTerms),x(totNumParams/2+1:totNumParams),x(end)*bgAmpl,x(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber)).^2;
        
        
        opt=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',10000,'Display','off');
        opt.MaxFunEvals=1000;
        %,'PlotFcns',@optimplotfval,'Display','iter');
        bestparams=zeros(1,totNumParams+1);
        bestfval=inf;
        startguess=zeros(1,totNumParams+1);
        wrong=0;
        
        for theG=1:numGuesses
            
            startguess(1:totNumParams/2)=1*abs(randn(1,totNumParams/2));
            startguess(totNumParams/2+1:totNumParams)=-10.^(kguessExpMin+rand(1,totNumParams/2)*(kguessExpMax-kguessExpMin));
            %startguess(totNumParams+1:totNumParams+numScales)=10*abs(randn(1,1));
            startguess(end)=10*abs(randn(1,1));
            try
            [params,fval]=lsqnonlin(funDiff,startguess,lb,ub,opt);
            if bestfval>fval
                bestfval=fval;
                bestparams=params;
            end
            catch
            end
        end
        
          bgResBoot(:,boot)=diffAutoCorrDiffMultScales3BgRes(ones(1,numTerms),bestparams(totNumParams/2+1:totNumParams),bestparams(end)*bgAmpl,bestparams(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber);
        
        % expFitRes{i}.aVals=bestparams{i}(1);
        % expFitRes{i}.kVals=bestparams{i}(2);
        bestparamsBoot(boot,:)=bestparams;
        bestfvalBoot(boot)=bestfval;
        
        autoCorrsDiffBoot(:,:,boot)=autoCorrsDiff;
    end
     expFitCombDiff{i}.boot.bestparamsBoot=bestparamsBoot;
      expFitCombDiff{i}.boot.bestfvalBoot=bestfvalBoot;
       expFitCombDiff{i}.boot.autoCorrsDiffBoot=autoCorrsDiffBoot;
            expFitCombDiff{i}.boot.bgResBoot=bgResBoot;
end

end

