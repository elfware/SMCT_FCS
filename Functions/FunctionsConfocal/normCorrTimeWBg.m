function [T] = normCorrTimeWBg(a,k,t1Vec,t2Vec,t3Vec,TBg,aVec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k2=k;
expT1=exp(k2*t1Vec);
vec=a*aVec;
numVec=TBg.*vec;
numerator = (1./(t2Vec-t1Vec)).*((exp(k2*t2Vec)-expT1)/k2)-exp(k2*t3Vec)+numVec;

denominator=expT1-exp(k2*t3Vec)+vec;

T= numerator./denominator;

end

