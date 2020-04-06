function [fN] = fN2TermsMultScalesver2(gUnNorm,fUnNorm,normNumber,lagTimesInds)

%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
[~,numT]=size(fUnNorm);

expVals=gUnNorm(lagTimesInds);
nom=fUnNorm-repmat(fUnNorm(:,end),[1 numT])+(expVals-repmat(expVals(:,end),[1 numT]));
denom=repmat(mean(fUnNorm(:,1:normNumber),2),[1 numT])-repmat(fUnNorm(:,end),[1 numT])+(repmat(mean(expVals(:,1:normNumber),2),[1 numT])-repmat(expVals(:,end),[1 numT]));

fN=nom./denom;
fN=fN';
fN=fN(:);

end

