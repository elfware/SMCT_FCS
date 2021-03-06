function [T] = normCorrTimeExpBg(aBg,kBg,k,t1Vec,t2Vec,t3Vec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%k2=-k*D/p^2;
k2=k;
expT1=exp(k2*t1Vec);
expBgT1=exp(kBg*t1Vec);

expT3=exp(k2*t3Vec);
expBgT3=exp(kBg*t3Vec);

aVec=(1./(t2Vec-t1Vec)).*(aBg*(exp(kBg*t2Vec)-expBgT1)/kBg)-aBg*expBgT3;
bVec=aBg*(expBgT1-expBgT3);
numerator = (1./(t2Vec-t1Vec)).*((exp(k2*t2Vec)-expT1)/k2)-expT3+aVec;

denominator=expT1-expT3+bVec;

T= numerator./denominator;

end

