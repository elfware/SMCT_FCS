function [P,phiTot] = gaussianPolarizationEqualRotationX(x,GsquaredX,numSamples,numRotationSteps)
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

B45=[1/sqrt(2) 1/sqrt(2) 0;1/sqrt(2) -1/sqrt(2) 0;0 0 1]; %Matrix to change base to the +/- 45 base 

meanPhi=meanAngle;
meanTheta=pi/2-meanAngle;
rotationAngles=linspace(0,2*pi,numRotationSteps);
Ihs=zeros(1,numSamples);
Ivs=zeros(1,numSamples);
I45pluses=zeros(1,numSamples);
I45minuses=zeros(1,numSamples);
phiTot=zeros(1,numRotationSteps*numSamples);
for i=1:numSamples
    phi=normrnd(meanPhi,sigmaAngle);%rand()*2*pi
    theta=pi/2;%normrnd(meanTheta,sigmaAngle);%acos(rand()*2-1)
    my=myEqualExciation(phi,theta);
    phiTot(1+numRotationSteps*(i-1):numRotationSteps*i)=atan2((my(2)*cos(rotationAngles)-my(3)*sin(rotationAngles)),my(1));
    for j=1:numRotationSteps
        rotAngle=rotationAngles(j);
        cR=cos(rotAngle);
        sR=sin(rotAngle);
        Rx=[1 0 0;0 cR -sR;0 sR cR];
        Rxmy=Rx*my;
        %Correct for lower excitation in z direction
        thetaAfterRot=acos(Rxmy(3));
        Rxmy=Rxmy*sin(thetaAfterRot);
       %Change of base for +/+ 45 base
        BRxmy=B45*Rxmy;
        
        Ih=Rxmy'*GsquaredX*Rxmy;
        Iv=Rxmy'*GsquaredY*Rxmy;
        I45plus=BRxmy'*GsquaredX*BRxmy;
        I45minus=BRxmy'*GsquaredY*BRxmy;
        Ihs(i)=Ihs(i)+Ih;
        Ivs(i)=Ivs(i)+Iv;
        I45pluses(i)=I45pluses(i)+I45plus;
        I45minuses(i)=I45minuses(i)+I45minus;
    end  
end

P=zeros(2,1);

Ihtot=sum(Ihs);
Ivtot=sum(Ivs);
I45plustot=sum(I45pluses);
I45minustot=sum(I45minuses);

P(1)=(Ihtot-Ivtot)/(Ihtot+Ivtot);
P(2)=(I45plustot-I45minustot)/(I45plustot+I45minustot);

end

