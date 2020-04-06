function [diff] = diffAutoCorrDiffMultScales4(aVec,kVec,bExtra,kFactor,autoCorrsDiff,lagTimes,normNumber,autoCorrScaledBg1,autoCorrScaledBg2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% v4 fitting with backgrounds
[numScales,numTimes]=size(lagTimes);
numTerms=numel(aVec)/numScales;

diff=zeros(1,numScales*numTimes);

 TheoryDiff=difffN2TermsExpMultScales([1 kVec(1) 1 kVec(1).*kFactor],bExtra(1)*autoCorrScaledBg1,bExtra(1)*autoCorrScaledBg2,normNumber,lagTimes);
 autoCorrsDiffTrans=autoCorrsDiff';
 diff=abs(TheoryDiff-autoCorrsDiffTrans(:));   
 
%for i=1:numScales
%  diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBg(aVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),kVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),-inf,lagTimes(i,:),lagTimes(i,end))-autoCorrs(i,:);
%  TheoryDiff=corrFuncMultExpBgbExtra4cleaned(aVec,kVec,bExtra(i),lagTimes(i,:),normNumber)-corrFuncMultExpBgbExtra4cleaned(aVec,kVec.*kFactor,bExtra(i),lagTimes(i,:),normNumber);
%  
%  diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=TheoryDiff-autoCorrsDiff(i,:);
%end
%%



end

