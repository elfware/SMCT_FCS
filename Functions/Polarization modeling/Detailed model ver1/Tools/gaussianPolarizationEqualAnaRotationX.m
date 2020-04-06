function [P,phiTot] = gaussianPolarizationEqualAnaRotationX(x,GsquaredX,numSamples)
%this function returns the polarasiation signals, P=[Phv;P45] defined by x=
%[meanAngle,sigmaAngle] for gaussian distributions of phi
%and theta. sigma for both thta and phi are here assumed to be equal and 
% meanTheta=pi/2-meanSigma. phi and theta defined according to Backer and Moerner 2014.
%numSamples is the number of samples to generate with this motecarlo metho

meanAngle=abs(x(1)); % if negative angles are given, use the magnitude of them insted (phi and thet are defined positive
sigmaAngle=abs(x(2));
GsquaredY=GsquaredX;
g1=GsquaredX(1,1);
g2=GsquaredX(2,2);
GsquaredY(1,1)=g2;
GsquaredY(2,2)=g1;


a=GsquaredX(1,1)-GsquaredX(2,2);
b=GsquaredX(1,1)+GsquaredX(2,2);
c=2*GsquaredX(3,3);

B45=[1/sqrt(2) 1/sqrt(2) 0;1/sqrt(2) -1/sqrt(2) 0;0 0 1]; %Matrix to change base to the +/- 45 base 

meanPhi=meanAngle;
meanTheta=pi/2-meanAngle;
numRotationSteps=1000;
rotationAngles=linspace(0,2*pi,numRotationSteps);
% Ihs=zeros(1,numSamples);
% Ivs=zeros(1,numSamples);
% I45pluses=zeros(1,numSamples);
% I45minuses=zeros(1,numSamples);
diffHV=0;
sumHV=0;
diff45=0;
sum45=0;
phiTot=zeros(1,numRotationSteps*numSamples);
for i=1:numSamples
    phi=rand()*2*pi;%normrnd(meanPhi,sigmaAngle);
    phi=mod(phi,2*pi);
    theta=acos(rand()*2-1);%normrnd(pi/2,0);%
    theta=mod(theta,pi);
    st=sin(theta);
    sp=sin(phi);
    
    ct=cos(theta);
    cp=cos(phi);
    
    my =st*[st*cp;st*sp;ct];
    my45=st*[(st*cp+st*sp)/sqrt(2);(st*cp-st*sp)/sqrt(2);ct];
        
    %Calcualte intensities from rotation around x axis analytically
    diffHV=diffHV+a*((2*pi*my(1)^2)-pi*(my(2)^2+my(3)^2));
    sumHV=sumHV+b*2*pi*my(1)^2+pi*(b+c)*(my(2)^2+my(3)^2);
    
    diff45=diff45+a*((2*pi*my45(1)^2)-pi*(my45(2)^2+my45(3)^2));
    sum45=sum45+b*2*pi*my45(1)^2+pi*(b+c)*(my45(2)^2+my45(3)^2);
    
    phiTot(1+numRotationSteps*(i-1):numRotationSteps*i)=atan2((my(2)*cos(rotationAngles)-my(3)*sin(rotationAngles)),my(1));
    
end

P=zeros(2,1);

% Ihtot=sum(Ihs);
% Ivtot=sum(Ivs);
% I45plustot=sum(I45pluses);
% I45minustot=sum(I45minuses);

P(1)=diffHV/sumHV;
P(2)=diff45/sum45;

end

