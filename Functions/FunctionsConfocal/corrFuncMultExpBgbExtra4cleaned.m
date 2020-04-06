function [corrfunc] = corrFuncMultExpBgbExtra4cleaned(constants,kVec,bExtra,tVec,normNumber)
%This function calcutes one of the terms used for fitting of the
%autocorrelation function. An normalised exponential with a background
%constant in the denominator (bExtra)
%exponentials are on the form constants(i)*exp()
numTerms=numel(constants);
numT=numel(tVec);
aVec=zeros(1,numT);


%Loop over the number of terms/exponentials. Should be 1 term for LacI
%Orientation manuscript.
for i=1:numTerms

k=kVec(i);
expT=exp(k*tVec);


aVec=aVec+constants(i)*expT;
end


numerator=aVec-aVec(end);

denominator=mean(numerator(1:normNumber))+bExtra;


corrfunc= numerator./denominator;
corrfunc(isnan(corrfunc))=0;
end

