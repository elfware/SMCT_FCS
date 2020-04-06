function [difference] = diffExpBgAndModel(params,autoCorrScaledExp,autoCorrScaledBg,lagTimes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
difference=autoCorrScaledExp-params(3)*autoCorrScaledBg-params(1)*exp(params(2)*lagTimes)-params(4);

difference=difference(:);

%diffNorm=autoCorrScaledExp-params(3)*autoCorrScaledBg+params(4);
%diffNorm=params(3)*autoCorrScaledBg+params(4);
%difference=difference./abs(diffNorm(:));
    




end

