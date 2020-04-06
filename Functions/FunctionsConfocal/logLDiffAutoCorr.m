function [logL] = logLDiffAutoCorr(aVec,kVec,bExtra,kFactor,autoCorrsDiff,lagTimes,normNumber,sig)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
numTimes=numel(lagTimes);
Nt=numTimes-1;% -1 because last point is always zero defintion wise, so should not been take into account when calculating the likelihood
 TheoryDiff=corrFuncMultExpBgbExtra(aVec,kVec,bExtra,-inf,lagTimes,lagTimes(end),normNumber)-corrFuncMultExpBgbExtra(aVec,kVec.*kFactor,bExtra,-inf,lagTimes,lagTimes(end),normNumber);
 
 squaredDiff=(autoCorrsDiff-TheoryDiff).^2;
 
 logL=-(1/(sqrt(2)*sig))*sum(squaredDiff)-Nt*log(sig*sqrt(2*pi));
 
 
end

