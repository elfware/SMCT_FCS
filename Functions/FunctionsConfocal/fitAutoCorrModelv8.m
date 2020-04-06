% Main function for fitting difference i autocorrelation
function [expFitCombDiff] = fitAutoCorrModelv8(theScales,compare,numBoots,numPreBoots,perc,scaleRes,tLagVals,pitchVals)
%UNTITLED4 Summary of this function goes here
% v2 Fixed so that sum squared deviation is minimzed
% v3 Fitting with background
% v4 fitting with background and background removed should look like
% difference of exponentials
%v8 only tyring to find where amplitude difference is biggest
%theScales=[5:26];
numTerms=1;
numScales=numel(theScales);
timeWind=size(scaleRes{theScales(1)}.combines{compare{1}(1)}.autoCorrValsAllok,2);



%expsToFit=[2 3 4 5 6 7];
%bgToUse=[12 13 14 15 16 17];
%numFits=numel(expsToFit);

numFits=numel(compare);

for i =1:numFits
    
    lagTimes=zeros(numScales,timeWind);
    corr1=zeros(numScales,timeWind);
    corr2=zeros(numScales,timeWind);
    autoCorrsDiff=zeros(numScales,timeWind);
    for j=1:numScales
        corr1(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(1)}.autoCorrValsAllok,perc);
        corr2(j,:)=trimmean(scaleRes{theScales(j)}.combines{compare{i}(2)}.autoCorrValsAllok,perc);
        autoCorrsDiff(j,:)= corr1(j,:)-corr2(j,:);
        lagTimes(j,:)=scaleRes{theScales(j)}.lagTimes;
    end
    [~,maxInd]=max(autoCorrsDiff);
    posAutoCorrsDiff=autoCorrsDiff.*(autoCorrsDiff>0)
    if maxInd==1
        theNorm=posAutoCorrsDiff(maxInd)+posAutoCorrsDiff(maxInd+1);
        theMaxAmplTLag=(lagTimes(maxInd)*posAutoCorrsDiff(maxInd)+lagTimes(maxInd+1)*posAutoCorrsDiff(maxInd+1))/theNorm;
    elseif maxInd==timeWind
        theNorm=posAutoCorrsDiff(maxInd)+posAutoCorrsDiff(maxInd-1);
        theMaxAmplTLag=(lagTimes(maxInd)*posAutoCorrsDiff(maxInd)+lagTimes(maxInd-1)*posAutoCorrsDiff(maxInd-1))/theNorm;
    else
        theNorm=posAutoCorrsDiff(maxInd-1)+posAutoCorrsDiff(maxInd)+posAutoCorrsDiff(maxInd+1);
        theMaxAmplTLag=(lagTimes(maxInd-1)*posAutoCorrsDiff(maxInd-1)+lagTimes(maxInd)*posAutoCorrsDiff(maxInd)+lagTimes(maxInd+1)*posAutoCorrsDiff(maxInd+1))/theNorm;
    end
    
    [~,tInd]=min(abs(tLagVals-theMaxAmplTLag));
    
    expFitCombDiff{i}.thePitch=pitchVals(tInd);
    
    thePitchBoot=zeros(1,numBoots);
    autoCorrsDiffBoot=cell(1,numBoots);
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
        
        [~,maxInd]=max(autoCorrsDiff);
        posAutoCorrsDiff=autoCorrsDiff.*(autoCorrsDiff>0)
        
        if maxInd==1
            theNorm=posAutoCorrsDiff(maxInd)+posAutoCorrsDiff(maxInd+1);
            theMaxAmplTLag=(lagTimes(maxInd)*posAutoCorrsDiff(maxInd)+lagTimes(maxInd+1)*posAutoCorrsDiff(maxInd+1))/theNorm;
        elseif maxInd==timeWind
            theNorm=posAutoCorrsDiff(maxInd)+posAutoCorrsDiff(maxInd-1);
            theMaxAmplTLag=(lagTimes(maxInd)*posAutoCorrsDiff(maxInd)+lagTimes(maxInd-1)*posAutoCorrsDiff(maxInd-1))/theNorm;
        else
            theNorm=posAutoCorrsDiff(maxInd-1)+posAutoCorrsDiff(maxInd)+posAutoCorrsDiff(maxInd+1);
            theMaxAmplTLag=(lagTimes(maxInd-1)*posAutoCorrsDiff(maxInd-1)+lagTimes(maxInd)*posAutoCorrsDiff(maxInd)+lagTimes(maxInd+1)*posAutoCorrsDiff(maxInd+1))/theNorm;
        end
        
        [~,tInd]=min(abs(tLagVals-theMaxAmplTLag));
        thePitchBoot(boot)=pitchVals(tInd);
        autoCorrsDiffBoot{boot}=autoCorrsDiff;
    end
    
    expFitCombDiff{i}.boot.thePitchBoot=thePitchBoot;
    expFitCombDiff{i}.boot.autoCorrsDiff=autoCorrsDiffBoot;
    
end




end

