function [fN] = fN2TermsMultScales(gUnNorm,fUnNorm,timeWind,normNumber)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
numScales=numel(fUnNorm)-timeWind+1;
inds=1:timeWind;
fN=zeros(1,numel(fUnNorm));
for i=1:numScales
    currInds=inds+(i-1);
   nom=(fUnNorm(currInds)-fUnNorm(currInds(end)))+(gUnNorm(currInds)-gUnNorm(currInds(end)));
   denom=mean(fUnNorm(currInds(1:normNumber))+gUnNorm(currInds(1:normNumber))-fUnNorm(currInds(end))-fUnNorm(currInds(end)));
    fN(1+(i-1)*timeWind:i*timeWind)=nom./denom;
end
end

