function [P,allMy,allTheta] = uniformPatchPolarizationRotationX(phiCenter,thetaCenter,Fp,GsquaredX,numSamples,numRotationSteps)
%this function returns the polarasiation signals, P=[Phv;P45] defined by x=
%[meanAngle,sigmaAngle] for gaussian distributions of phi
%and theta. sigma for both thta and phi are here assumed to be equal and 
% meanTheta=pi/2-meanSigma. phi and theta defined according to Backer and Moerner 2014.
%numSamples is the number of samples to generate with this motecarlo metho


center=[cos(phiCenter)*sin(thetaCenter);sin(phiCenter)*sin(thetaCenter);cos(thetaCenter)];
cosNegPhiCenter=cos(-phiCenter);
sinNegPhiCenter=sin(-phiCenter);
GsquaredY=GsquaredX;
g1=GsquaredX(1,1);
g2=GsquaredX(2,2);
GsquaredY(1,1)=g2;
GsquaredY(2,2)=g1;

B45=[1/sqrt(2) 1/sqrt(2) 0;1/sqrt(2) -1/sqrt(2) 0;0 0 1]; %Matrix to change base to the +/- 45 base 

rotationAngles=linspace(0,2*pi,numRotationSteps);
Ihs=zeros(1,numSamples);
Ivs=zeros(1,numSamples);
I45pluses=zeros(1,numSamples);
I45minuses=zeros(1,numSamples);
allMy=zeros(3,numSamples);
allTheta=zeros(1,numSamples);
for i=1:numSamples
    u=rand(); %random number between 0 and 1
    v=(1-Fp)+rand()*Fp; % rnadom numberbetween 1-Fp and 1
    
    phi_patch=2*pi*u;%phi_patch between [0 2pi] uniform over patch area
    %2*v-1 should go from [1-2*Fp to 1] since the surfface area of sampled area of the sampled patch is is 4*pi*r^2*Fp which gives theta_P_max = acos(1-2*FP) and theta_P_min=0=acos(1)
    theta_patch=acos(2*v-1);  %thetapatch between [0 acos(1-2*Fp] over the patch area. this Interval is consistent with that the patch area = Fp/total surface area of spere
    allTheta(i)=theta_patch;
    %coordinates in the patch coordiante system
    x_patch=cos(phi_patch)*sin(theta_patch);
    y_patch=sin(phi_patch)*sin(theta_patch);
    z_patch=cos(theta_patch);
    
    baseVec3=center;
    
    %Rotate down center to phi=0 (around z-axis) to get consistent
    %definition of baseVec1 and baseVec2.
    RotZ=[cosNegPhiCenter -sinNegPhiCenter 0;sinNegPhiCenter cosNegPhiCenter 0;0 0 1];
    baseVec3rot=RotZ*baseVec3;
    %After the rotation the y axis is always orthogonal to baseVec3rot
    baseVec2rot=[0;1;0];
    %The third base vector is the cross product of the two already found
    baseVec1rot=cross(baseVec3rot,baseVec2rot);
    %Rotate back the base vectors
    RotZback=[cosNegPhiCenter sinNegPhiCenter 0;-sinNegPhiCenter cosNegPhiCenter 0;0 0 1];
    baseVec1=RotZback*baseVec1rot;
    baseVec2=RotZback*baseVec2rot;
    %get cartesian coordinates in normal "sphere" coordiante system
    my = baseVec1*x_patch + baseVec2*y_patch + baseVec3*z_patch;
    allMy(:,i)=my;
%     x=coords(1);
%     y=coords(2);
%     z=coords(3);
%     %get phi_P and theta_P in this normal coordiante ssytem fro mteh caresian
%     %coordiantes
%     phi_P = atan(y/x);
%     theta_P = acos(z);
%     phi=normrnd(meanPhi,sigmaAngle);%rand()*2*pi
%     theta=pi/2;%normrnd(meanTheta,sigmaAngle);%acos(rand()*2-1)
%     my=myEqualExciation(phi,theta);
%     phiTot(1+numRotationSteps*(i-1):numRotationSteps*i)=atan2((my(2)*cos(rotationAngles)-my(3)*sin(rotationAngles)),my(1));
    for j=1:numRotationSteps
        if numRotationSteps>1
            rotAngle=j*(2*pi/numRotationSteps);
        else
            rotAngle=0;
        end
        cR=cos(rotAngle);
        sR=sin(rotAngle);
        Rx=[1 0 0;0 cR -sR;0 sR cR];
        Rxmy=Rx*my;
        %Correct for lower excitation in z direction
        %thetaAfterRot=acos(Rxmy(3));
        %Rxmy=Rxmy*sin(thetaAfterRot);
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

