function [T] = normCorrTimeMultExpBgv2(aBgVec,kBgVec,k,t1Vec,t2Vec,t3Vec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%k2=-k*D/p^2;
numBg=numel(aBgVec);
numT=numel(t1Vec);
k2=k;
expT1=exp(k2*t1Vec);
expT3=exp(k2*t3Vec);
aVec=zeros(1,numT);
bVec=zeros(1,numT);
for i=1:numBg
aBg=aBgVec(i);
kBg=kBgVec(i);

expBgT1=exp(kBg*t1Vec);
expBgT3=exp(kBg*t3Vec);

aVec=aVec+(1./(t2Vec-t1Vec)).*(aBg*(exp(kBg*t2Vec)-expBgT1)/kBg)-aBg*expBgT3;
bVec=bVec+aBg*(expBgT1-expBgT3);
end


numerator = (1./(t2Vec-t1Vec)).*((exp(k2*t2Vec)-expT1)/k2)-expT3;

denominator=expT1-expT3+bVec;

T= numerator./denominator+aVec./bVec;



end

