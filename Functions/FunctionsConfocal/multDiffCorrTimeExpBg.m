function [diff] = multDiffCorrTimeExpBg(aBg,kBg,kVec,t1Vec,t2Vec,t3Vec,vals,inds)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

numE=numel(kVec);
numT=numel(t1Vec);
diff=zeros(1,numT);

for i=1:numE
  theExp=inds(i);  
  corrT=vals(theExp,:);
  
  diff=diff+abs(normCorrTimeExpBg(aBg,kBg,kVec(i),t1Vec,t2Vec,t3Vec)-corrT);
  
end


end

