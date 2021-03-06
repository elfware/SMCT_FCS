function [diff] = diffAutoCorrDiffMultScales2(aVec1,kVec1,aVec2,kVec2,bExtra,kFactor,autoCorrsDiff,lagTimes,normNumber)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[numScales,numTimes]=size(lagTimes);
numTerms=numel(aVec1)/numScales;

diff=zeros(1,numScales*numTimes);

for i=1:numScales
  %diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBg(aVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),kVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),-inf,lagTimes(i,:),lagTimes(i,end))-autoCorrs(i,:);
  TheoryDiff=corrFuncMultExpBgbExtra(aVec1,kVec1,bExtra(i),-inf,lagTimes(i,:),lagTimes(i,end),normNumber)-corrFuncMultExpBgbExtra(aVec2,kVec2.*kFactor,bExtra(i),-inf,lagTimes(i,:),lagTimes(i,end),normNumber);
  
  diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=TheoryDiff-autoCorrsDiff(i,:);
end


end

