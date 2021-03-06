function [corrfunc] = corrFuncMultExpBg(aBgVec,kBgVec,k,tVec,t3,normNumber)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%k2=-k*D/p^2;
numBg=numel(aBgVec);
numT=numel(tVec);
k2=k;
expT1=exp(k2*tVec);
expT3=exp(k2*t3);
aVec=zeros(1,numT);
b=0;
for i=1:numBg
aBg=aBgVec(i);
kBg=kBgVec(i);

expBgT1=exp(kBg*tVec);
expBgT3=exp(kBg*t3);

%aVec=aVec+(1./(t2Vec-t1Vec)).*(aBg*(exp(kBg*t2Vec)-expBgT1)/kBg)-aBg*expBgT3;
%bVec=bVec+aBg*(expBgT1-expBgT3);

aVec=aVec+aBg*(expBgT1-expBgT3);
b=b+aBg*(mean(expBgT1(1:normNumber))-expBgT3);
end


numerator = expT1-expT3+aVec;

denominator=mean(expT1(1:normNumber))-expT3+b;

corrfunc= numerator./denominator;
corrfunc(isnan(corrfunc))=0;
end

