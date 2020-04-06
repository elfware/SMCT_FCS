function [diff] = multDiffCorrTime(aVec,bVec,kVec,t1Vec,t2Vec,t3Vec,vals,inds)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

numE=numel(kVec);
numT=numel(t1Vec);
diff=zeros(1,numT);

for i=1:numE
  theExp=inds(i);  
  corrT=vals(theExp,:);
  
  diff=diff+abs(normCorrTime(aVec,bVec,kVec(theExp),t1Vec,t2Vec,t3Vec)-corrT);
  
end


end

