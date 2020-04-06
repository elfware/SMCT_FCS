function [P] = gaussianPolarization(x,GsquaredX)
%this function returns the polarasiation signals, P=[Phv;P45] defined by x=
%[meanAngle,sigmaAngle] for equal gaussian distributions of phi
%and theta. phi and theta defined according to Backer and Moerner 2014. 

meanAngle=x(1);
sigmaAngle=x(2);


a=GsquaredX(1,1)-GsquaredX(2,2);
b=GsquaredX(1,1)+GsquaredX(2,2);
c=2*GsquaredX(3,3);

factor=exp(-2*(sigmaAngle^2));
integral = factor*cos(2*meanAngle);
integralsin = factor*sin(2*meanAngle);

P=zeros(2,1);

k=(a-integral)/(b+c+(c-b)*integral);

P(1)=k*integral;
P(2)=k*integralsin;

end

