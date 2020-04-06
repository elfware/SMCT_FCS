function [thediff] = diffAutoCorrDiffMultScales3BgRes(aVec,kVec,bExtra,kFactor,autoCorrsDiff,lagTimes,normNumber)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[numScales,numTimes]=size(lagTimes);
numTerms=numel(aVec)/numScales;

thediff=zeros(1,numScales);

for i=1:numScales
  %diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBg(aVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),kVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),-inf,lagTimes(i,:),lagTimes(i,end))-autoCorrs(i,:);
  thediff(i)=corrFuncMultExpBgbExtra4BgRes(aVec,kVec,bExtra(i),-inf,lagTimes(i,:),lagTimes(i,end),normNumber)-corrFuncMultExpBgbExtra4BgRes(aVec,kVec.*kFactor,bExtra(i),-inf,lagTimes(i,:),lagTimes(i,end),normNumber);
  
 
end


end

