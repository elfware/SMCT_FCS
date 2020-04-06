function [T] = normCorrTime(aVec,bVec,k,t1Vec,t2Vec,t3Vec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%k2=-k*D/p^2;
k2=k;
expT1=exp(k2*t1Vec);
numerator = (1./(t2Vec-t1Vec)).*((exp(k2*t2Vec)-expT1)/k2)-exp(k2*t3Vec)+aVec;

denominator=expT1-exp(k2*t3Vec)+bVec;

T= numerator./denominator;

end

