% Main function for fitting difference i autocorrelation
function [expFitCombDiff] = fitAutoCorrModelv9_bgguess(timeWind,normNumber,kguessExpMax,kguessExpMin,numGuesses,theScales,compare,lb,ub,numBoots,numPreBoots,perc,scaleRes,combUnNorm,bgInds,aguessExpMax,aguessExpMin)
%UNTITLED4 Summary of this function goes here
% v2 Fixed so that sum squared deviation is minimzed
% v3 Fitting with background
% v4 fitting with background and background removed should look like
% difference of exponentials
%v9 removed so dthe fitting oesnt consider exact equation of autocorrelation (eq  22 in first Nature submission)
%but ony approximate form( eq 23 in first Nature submission)
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
    corr1=zeros(numScales,timeWind);
    corr2=zeros(numScales,timeWind);
%     bgAmpl=zeros(1,numScales);
     for j=1:numScales
        corr1(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok,perc);
        corr2(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok,perc);
        autoCorrsDiff(j,:)= corr1(j,:)-corr2(j,:);
         lagTimes(j,:)=scaleRes{theScales(j)}.lagTimes;
%         bgAmpl(j)=1;%combUnNorm{i+3}.autoCorrScaled(j,1)-combUnNorm{i+3}.autoCorrScaled(j,end);
         
    end
    
    %bgAmpl=bgAmpl/sum(bgAmpl);
    
    bgAmpl=ones(1,numScales);
    
    %Function to minimize definition
    
    autoCorrScaledBg1=combUnNorm{bgInds(1)}.autoCorrScaled(theScales-combUnNorm{bgInds(1)}.startScale+1,:);
    autoCorrScaledBg2=combUnNorm{bgInds(2)}.autoCorrScaled(theScales-combUnNorm{bgInds(2)}.startScale+1,:);
    funDiff = @(x) sqrt((weights.*diffAutoCorrDiffandExpMultScales9(ones(1,numTerms),x(totNumParams/2+1:totNumParams),x(end)*bgAmpl,x(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber,autoCorrScaledBg1,autoCorrScaledBg2)).^2);
    
    
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
        %startguess(end)=10*abs(randn(1,1));
        startguess(end)=10.^(aguessExpMin+rand(1,1)*(aguessExpMax-aguessExpMin));
      try
        [params,fval]=lsqnonlin(funDiff,startguess,lb,ub,opt);
        if bestfval(i)>fval
            bestfval(i)=fval;
            bestparams{i}=params;
        end
      catch
            wrong=wrong+1;
      end 
        
    end
    
       
 
    correctedCorr1=autoCorrBgRem([bestparams{i}(2) bestparams{i}(3)],corr1,autoCorrScaledBg1,normNumber,lagTimes);
    correctedCorr2=autoCorrBgRem([bestparams{i}(1)*bestparams{i}(2) bestparams{i}(3)],corr2,autoCorrScaledBg2,normNumber,lagTimes);
   
   % bgRes=diffAutoCorrDiffMultScales3BgRes(ones(1,numTerms),bestparams{i}(totNumParams/2+1:totNumParams),bestparams{i}(end)*bgAmpl,bestparams{i}(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber);
    
    
    % expFitRes{i}.aVals=bestparams{i}(1);
    % expFitRes{i}.kVals=bestparams{i}(2);
    expFitCombDiff{i}.bestparams=bestparams{i};
    expFitCombDiff{i}.bestfval=bestfval(i);
    expFitCombDiff{i}.corr1=corr1;
    expFitCombDiff{i}.corr2=corr2;
    expFitCombDiff{i}.autoCorrsDiff=autoCorrsDiff;
    expFitCombDiff{i}.lagTimes=lagTimes;
    expFitCombDiff{i}.bgAmpl=bgAmpl;
    %expFitCombDiff{i}.stdDiffs=stdDiff;
    expFitCombDiff{i}.weights=weights;
   % expFitCombDiff{i}.bgRes=bgRes;
   expFitCombDiff{i}.correctedCorr1=correctedCorr1;
   expFitCombDiff{i}.correctedCorr2=correctedCorr2;
   
    %%Bootstrapping
    
    bestparamsBoot=zeros(numBoots,totNumParams+1);
    bestfvalBoot=zeros(1,numBoots);
    autoCorrsDiffBoot=zeros(numScales,timeWind,numBoots);
    bgResBoot=zeros(numScales,numBoots);
    corr1Boot=zeros(numScales,timeWind,numBoots);
    corr2Boot=zeros(numScales,timeWind,numBoots);
    correctedCorr1Boot=zeros(numScales,timeWind,numBoots);
    correctedCorr2Boot=zeros(numScales,timeWind,numBoots);
    parfor boot=1:numBoots
        autoCorrsDiff=zeros(numScales,timeWind);
        corr1Temp=zeros(numScales,timeWind);
        corr2Temp=zeros(numScales,timeWind);
        for j=1:numScales
            [numtrajs1,~]=size(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok);
            [numtrajs2,~]=size(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok);
            bootInds1=ceil(rand(1,numtrajs1)*numtrajs1);
            bootInds2=ceil(rand(1,numtrajs2)*numtrajs2);
            
            corr1Temp(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok(bootInds1,:),perc);
            corr2Temp(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok(bootInds2,:),perc);
            
            autoCorrsDiff(j,:)=corr1Temp(j,:)-corr2Temp(j,:);
          %  lagTimes(j,:)=scaleRes{theScales(j)}.lagTimes;
            
        end
        
        
        funDiff = @(x) sqrt((weights.*diffAutoCorrDiffandExpMultScales9(ones(1,numTerms),x(totNumParams/2+1:totNumParams),x(end)*bgAmpl,x(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber,autoCorrScaledBg1,autoCorrScaledBg2)).^2);
        
        
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
            %startguess(end)=10*abs(randn(1,1));
            startguess(end)=10.^(aguessExpMin+rand(1,1)*(aguessExpMax-aguessExpMin));
            try
            [params,fval]=lsqnonlin(funDiff,startguess,lb,ub,opt);
            if bestfval>fval
                bestfval=fval;
                bestparams=params;
            end
            catch
                wrong=wrong+1;
            end
        end
        
        correctedCorr1=autoCorrBgRem([bestparams(2) bestparams(3)],corr1Temp,autoCorrScaledBg1,normNumber,lagTimes);
        correctedCorr2=autoCorrBgRem([bestparams(1)*bestparams(2) bestparams(3)],corr2Temp,autoCorrScaledBg2,normNumber,lagTimes);
   
        
       %   bgResBoot(:,boot)=diffAutoCorrDiffMultScales3BgRes(ones(1,numTerms),bestparams(totNumParams/2+1:totNumParams),bestparams(end)*bgAmpl,bestparams(1:totNumParams/2),autoCorrsDiff,lagTimes,normNumber);
        
        % expFitRes{i}.aVals=bestparams{i}(1);
        % expFitRes{i}.kVals=bestparams{i}(2);
        bestparamsBoot(boot,:)=bestparams;
        bestfvalBoot(boot)=bestfval;
        
        autoCorrsDiffBoot(:,:,boot)=autoCorrsDiff;
        corr1Boot(:,:,boot)=corr1Temp;
        corr2Boot(:,:,boot)=corr2Temp;
        correctedCorr1Boot(:,:,boot)=correctedCorr1;
        correctedCorr2Boot(:,:,boot)=correctedCorr2;
        
    end
     expFitCombDiff{i}.boot.bestparamsBoot=bestparamsBoot;
      expFitCombDiff{i}.boot.bestfvalBoot=bestfvalBoot;
       expFitCombDiff{i}.boot.autoCorrsDiffBoot=autoCorrsDiffBoot;
        expFitCombDiff{i}.boot.corr1Boot=corr1Boot;
        expFitCombDiff{i}.boot.corr2Boot=corr2Boot;
        expFitCombDiff{i}.boot.correctedCorr1Boot=correctedCorr1Boot;
        expFitCombDiff{i}.boot.correctedCorr2Boot=correctedCorr2Boot;
        
       %      expFitCombDiff{i}.boot.bgResBoot=bgResBoot;
end

end

