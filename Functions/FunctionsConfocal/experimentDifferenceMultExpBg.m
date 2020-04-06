function [diff] = experimentDifferenceMultExpBg(aBgVec,kBgVec,kVec,t1Vec,t2Vec,t3Vec,vals)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


theoryVals1=normCorrTimeMultExpBgv2(aBgVec,kBgVec,kVec(1),t1Vec,t2Vec,t3Vec);
theoryVals2=normCorrTimeMultExpBgv2(aBgVec,kBgVec,kVec(2),t1Vec,t2Vec,t3Vec);
theoryDiff=theoryVals1-theoryVals2;

diff=theoryDiff-vals;

end

