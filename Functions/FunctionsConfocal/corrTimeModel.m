function [corrTimes] = corrTimeModel(params,k1,minTime,maxTime,Dxs)
%params= [y a b]
y=params(1); %Pitch paramter, um/revolution (2*pi radians)
expBG1=exp(params(2)*minTime); %Backgrounparameter
expBGend=exp(params(2)*maxTime); %Backgroundparameter

a=params(2);
b=params(3);

%a = (params(3)/params(2))*(expBGend-expBG1)-params(3)*expBGend*(maxTime-minTime);
%b = params(3)*(expBG1 - expBGend);

ynew=y/(2*pi); % 

factors=ynew^2./(Dxs.*k1);
exponents1=exp(maxTime./factors);
exponents2=exp(minTime./factors);
corrTimesNumerator=factors.*(exponents1-exponents2)-(maxTime-minTime).*exponents1;
corrTimesDenominator=exponents2-exponents1+b;

corrTimes=(corrTimesNumerator./corrTimesDenominator)+a./corrTimesDenominator;

end

