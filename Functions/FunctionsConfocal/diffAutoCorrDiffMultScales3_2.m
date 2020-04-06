function [diff] = diffAutoCorrDiffMultScales3_2(aVec,kVec,bExtra1,bExtra2,kFactor,autoCorrsDiff,lagTimes,normNumber)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[numScales,numTimes]=size(lagTimes);
numTerms=numel(aVec)/numScales;

diff=zeros(1,numScales*numTimes);

for i=1:numScales
  %diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBg(aVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),kVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),-inf,lagTimes(i,:),lagTimes(i,end))-autoCorrs(i,:);
  TheoryDiff=corrFuncMultExpBgbExtra4cleaned(aVec,kVec,bExtra1(i),lagTimes(i,:),normNumber)-corrFuncMultExpBgbExtra4cleaned(aVec,kVec.*kFactor,bExtra2(i),lagTimes(i,:),normNumber);
  
  diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=TheoryDiff-autoCorrsDiff(i,:);
end


end

