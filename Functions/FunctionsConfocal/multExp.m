function [func] = multExp(aVec,kVec,tVec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numTerms=numel(aVec);
numTimes=numel(tVec);
func =zeros(1,numTimes);
for i=1:numTerms
    func=func+aVec(i)*exp(kVec(i)*tVec);
end

end

