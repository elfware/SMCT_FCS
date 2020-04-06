function [diff] = diffAutoCorrDiffandExpMultScales8(aVec,kVec,bExtra,kFactor,autoCorrsDiff,lagTimes,normNumber,autoCorrScaledBg1,autoCorrScaledBg2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% v4 fitting with backgrounds
% v7 average of amplitued
% v8 indiviudal amplitudes
%v9 removed so doesnt consider exact equation (eq  22 in first Nature submission)
%but ony approximate form( eq 23 in first Nature submission)
[numScales,numTimes]=size(lagTimes);

diff=zeros(1,numScales*numTimes);

 %TheoryDiff=difffN2TermsExpMultScales([1 kVec(1) 1 kVec(1).*kFactor],bExtra(1)*autoCorrScaledBg1,bExtra(1)*autoCorrScaledBg2,normNumber,lagTimes);
 autoCorrsDiffTrans=autoCorrsDiff';
 %diff=TheoryDiff-autoCorrsDiffTrans(:);   
 
 %corrected=autoCorrBgRem([kVec(1) bExtra(1)],corr1,autoCorrScaledBg1,normNumber,lagTimes);
 %corrected2=autoCorrBgRem([kVec(1).*kFactor bExtra(1)],corr2,autoCorrScaledBg2,normNumber,lagTimes);
 
 %correctedDiff=corrected-corrected2;
% bNew=bExtra.*mean([(autoCorrScaledBg1(:,1)-autoCorrScaledBg1(:,end)) (autoCorrScaledBg2(:,1)-autoCorrScaledBg2(:,end))],2)';
 bNew1=bExtra.*[autoCorrScaledBg1(:,1)-autoCorrScaledBg1(:,end)];
 bNew2=bExtra.*[autoCorrScaledBg2(:,1)-autoCorrScaledBg2(:,end)];
 TheoryDiffCorrected=diffAutoCorrDiffMultScales3_2(1,kVec(1),bNew1,bNew2,kFactor,autoCorrsDiff,lagTimes,normNumber);
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

