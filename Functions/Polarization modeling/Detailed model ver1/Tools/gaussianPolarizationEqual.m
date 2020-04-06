function [P] = gaussianPolarizationEqual(x,GsquaredX)
%this function returns the polarasiation signals, P=[Phv;P45] defined by x=
%[meanAngle,sigmaAngle] for gaussian distributions of phi
%and theta. sigma for both thta and phi are here assumed to be equal and 
% meanTheta=pi/2-meanSigma. phi and theta defined according to Backer and Moerner 2014.
%solution is without exciation scaling factor for my (same excitation in z as in x and y) 

meanAngle=abs(x(1)); % if negative angles are given, use the magnitude of them insted (phi and thet are defined positive
sigmaAngle=abs(x(2));


a=GsquaredX(1,1)-GsquaredX(2,2);
b=GsquaredX(1,1)+GsquaredX(2,2);
c=2*GsquaredX(3,3);

factor=exp(-2*(sigmaAngle^2));
integralPhicos = factor*cos(2*meanAngle);
integralPhisin = factor*sin(2*meanAngle);
integralTheta=factor*cos(2*(pi/2-meanAngle));
P=zeros(2,1);

k=(a-a*integralTheta)/(b+c+(c-b)*integralTheta);

P(1)=k*integralPhicos;
P(2)=k*integralPhisin;

end

