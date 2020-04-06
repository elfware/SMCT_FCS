function [diff] = multDiffCorrTimeMultExpBgv2(aBgVec,kBgVec,kVec,t1Vec,t2Vec,t3Vec,vals,inds)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

numE=numel(kVec);
numT=numel(t1Vec);
diff=zeros(1,numT);

for i=1:numE
  theExp=inds(i);  
  corrT=vals(theExp,:);
  
  %diff=diff+abs(normCorrTimeMultExpBgv2(aBgVec,kBgVec,kVec(i),t1Vec,t2Vec,t3Vec)-corrT);
  diff((i-1)*numT+1:(i-1)*numT+numT)=normCorrTimeMultExpBgv2(aBgVec,kBgVec,kVec(i),t1Vec,t2Vec,t3Vec)-corrT;
  
end


end

