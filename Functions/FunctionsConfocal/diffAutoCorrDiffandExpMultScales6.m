function [diff] = diffAutoCorrDiffandExpMultScales6(aVec,kVec,bExtra,kFactor,autoCorrsDiff,lagTimes,normNumber,autoCorrScaledBg1,autoCorrScaledBg2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% v4 fitting with backgrounds
[numScales,numTimes]=size(lagTimes);

diff=zeros(1,numScales*numTimes);

 TheoryDiff=difffN2TermsExpMultScales([1 kVec(1) 1 kVec(1).*kFactor],bExtra(1)*autoCorrScaledBg1,bExtra(1)*autoCorrScaledBg2,normNumber,lagTimes);
 autoCorrsDiffTrans=autoCorrsDiff';
 diff=TheoryDiff-autoCorrsDiffTrans(:);   
 
 %corrected=autoCorrBgRem([kVec(1) bExtra(1)],corr1,autoCorrScaledBg1,normNumber,lagTimes);
 %corrected2=autoCorrBgRem([kVec(1).*kFactor bExtra(1)],corr2,autoCorrScaledBg2,normNumber,lagTimes);
 
 %correctedDiff=corrected-corrected2;
 
 TheoryDiffCorrected=diffAutoCorrDiffMultScales3(1,kVec(1),bExtra,kFactor,autoCorrsDiff,lagTimes,normNumber);
 TheoryDiffCorrectedTrans=TheoryDiffCorrected';
 %correctedDiffTrans=correctedDiff';
 
 %diff=abs(diff)+abs(TheoryDiffCorrectedTrans-correctedDiffTrans(:));
  diff=abs(diff)+abs(TheoryDiffCorrectedTrans-autoCorrsDiffTrans(:));
 diff(isnan(diff))=inf;
 
%for i=1:numScales
%  diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=corrFuncMultExpBg(aVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),kVec((i-1)*numTerms+1:(i-1)*numTerms+numTerms),-inf,lagTimes(i,:),lagTimes(i,end))-autoCorrs(i,:);
%  TheoryDiff=corrFuncMultExpBgbExtra4cleaned(aVec,kVec,bExtra(i),lagTimes(i,:),normNumber)-corrFuncMultExpBgbExtra4cleaned(aVec,kVec.*kFactor,bExtra(i),lagTimes(i,:),normNumber);
%  
%  diff((i-1)*numTimes+1:(i-1)*numTimes+numTimes)=TheoryDiff-autoCorrsDiff(i,:);
%end
%%



end

