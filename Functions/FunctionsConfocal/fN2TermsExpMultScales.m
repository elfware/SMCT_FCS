function [fN] = fN2TermsExpMultScales(params,fUnNorm,normNumber,lagTimes)

%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
[~,numT]=size(fUnNorm);

expVals=exp(params(2)*lagTimes);
nom=params(1)*(expVals-repmat(expVals(:,end),[1 numT]))+fUnNorm-repmat(fUnNorm(:,end),[1 numT]);
%themeanF=fUnNorm(:,1:normNumber)-repmat(fUnNorm(:,end),[normNumber numT]);

%denom=repmat(mean(fUnNorm(:,1:normNumber),2),[1 numT])-repmat(fUnNorm(:,end),[1 numT])+params(1)*(repmat(mean(expVals(:,1:normNumber),2),[1 numT])-repmat(expVals(:,end),[1 numT]));

denom=repmat(mean(nom(:,1:normNumber),2),[1 numT]);
fN=nom./denom;
fN=fN';
fN=fN(:);

end

