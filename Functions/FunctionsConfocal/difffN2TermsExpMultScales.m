function [difference] = difffN2TermsExpMultScales(params,autoCorrScaledBg1,autoCorrScaledBg2,normNumber,lagTimes)
%params(1) = a firstt term
%params(2)= k firstterm
difference=fN2TermsExpMultScales([params(1) params(2)],autoCorrScaledBg1,normNumber,lagTimes)-fN2TermsExpMultScales([params(3) params(4)],autoCorrScaledBg2,normNumber,lagTimes);

end

