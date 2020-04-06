function [diff] = diffAutoCorrMultScales5(aVec,kVec,bExtra,autoCorrs,lagTimes,normNumber)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[numScales,numTimes]=size(lagTimes);
numTerms=numel(aVec)/numScales;

diff=zeros(1,numScales*numTimes);

for i=1:numScales
  %diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBg(aVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),kVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),-inf,lagTimes(i,:),lagTimes(i,end))-autoCorrs(i,:);
  %diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBg(aVec,kVec,-inf,lagTimes(i,:),lagTimes(i,end),normNumber)-autoCorrs(i,:);
  diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBgbExtra3(aVec,kVec,bExtra(i),-inf,lagTimes(i,:),lagTimes(i,end),normNumber)-autoCorrs(i,:);
  %diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=fN2TermsExpMultScales([aVec kVec],autoCorrScaledBg(i,:),normNumber,lagTimes(i,:))'-autoCorrs(i,:);
end


end

