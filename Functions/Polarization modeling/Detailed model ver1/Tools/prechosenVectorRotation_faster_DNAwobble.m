function  [IhExact,IvExact,I45plusExact,I45minusExact,IhDiscrete,IvDiscrete,I45plusDiscrete,I45minusDiscrete] = prechosenVectorRotation_faster_DNAwobble(phiCenter,thetaCenter,Param)
%this function returns the polarasiation signals, P=[Phv;P45] defined by x=
%[meanAngle,sigmaAngle] for gaussian distributions of phi
%and theta. sigma for both thta and phi are here assumed to be equal and 
% meanTheta=pi/2-meanSigma. phi and theta defined according to Backer and Moerner 2014.
%numSamples is the number of samples to generate with this motecarlo metho
%
%All trajectories does the same rotation trajectory in this version of the
%simulation algorithm

StringCom=Param.StringCom;
GsquaredX=Param.GsquaredX;


rotationAxes =Param.rotationAxes;% [1;0;0]; % Every column is an axis around wich rotation is performed 
Dr = Param.Dr;%[5000]; %One dimensiaonal rotatioal diffusion constants for the different axes. In radians^2/s
dt=Param.dt;  %1e-6; %in seconds
dim=Param.dim;%1;
ExScale=Param.ExScale;%scaling coordinate for to include the exitation laser behavior  
averagePhotons=Param.averagePhotonsinStep;%40000;%Average number of photons per second. Number of photons in a timestep is
%assumed to be  poission distribtued.
binning=Param.binning ; %Number of timesteps to bin into one bin (sum intesieties), done after the generatio of the intensieties traces
numTrajs=Param.numTrajs;%20; % number of trajectoreis to simulate, much larger in real data
trajLength=Param.trajLength;%0.2; %Length of trajectories in seconds
burninTime=Param.burninTime;
%startPhi=Param.startPhi;
%startTheta=Param.startTheta;



timeSteps=0:dt:trajLength;
burninIndex=floor(numel(timeSteps)*(burninTime/trajLength)+1);
if isnan(burninIndex)
    burninIndex=1;
end
averagePhotonsinStep=averagePhotons*dt;



GsquaredY=GsquaredX;
g1=GsquaredX(1,1);
g2=GsquaredX(2,2);
GsquaredY(1,1)=g2;
GsquaredY(2,2)=g1;

B45=[1/sqrt(2) 1/sqrt(2) 0;1/sqrt(2) -1/sqrt(2) 0;0 0 1]; %Matrix to change base to the +/- 45 base 

GsquaredX45=B45'*GsquaredX*B45;
GsquaredY45=B45'*GsquaredY*B45;


IhDiscrete=zeros(numTrajs,numel(timeSteps));
IvDiscrete=zeros(numTrajs,numel(timeSteps));    
IhExact=zeros(numTrajs,numel(timeSteps));
IvExact=zeros(numTrajs,numel(timeSteps));
PhvExact=zeros(numTrajs,numel(timeSteps));
startMy=zeros(3,numTrajs);
for i=1:numTrajs
    
%     dnaPhi=randn()*DNAwobble;
%     phiCenterCurr=phiCenter+dnaPhi;
%     rotationAxes=[cos(dnaPhi);sin(dnaPhi);0];
%     
%     center=[cos(phiCenterCurr)*sin(thetaCenter);sin(phiCenterCurr)*sin(thetaCenter);cos(thetaCenter)];
%     cosNegPhiCenter=cos(-phiCenterCurr);
%     sinNegPhiCenter=sin(-phiCenterCurr);
%     
%     u=rand(); %random number between 0 and 1
%     v=(1-Fp)+rand()*Fp; % rnadom numberbetween 1-Fp and 1
%     
%     phi_patch=2*pi*u;%phi_patch between [0 2pi] uniform over patch area
%     %2*v-1 should go from [1-2*Fp to 1] since the surfface area of sampled area of the sampled patch is is 4*pi*r^2*Fp which gives theta_P_max = acos(1-2*FP) and theta_P_min=0=acos(1)
%     theta_patch=acos(2*v-1);  %thetapatch between [0 acos(1-2*Fp] over the patch area. this Interval is consistent with that the patch area = Fp/total surface area of spere
%     allTheta(i)=theta_patch;
%     %coordinates in the patch coordiante system
%     x_patch=cos(phi_patch)*sin(theta_patch);
%     y_patch=sin(phi_patch)*sin(theta_patch);
%     z_patch=cos(theta_patch);
%     
%     baseVec3=center;
%     
%     %Rotate down center to phi=0 (around z-axis) to get consistent
%     %definition of baseVec1 and baseVec2.
%     RotZ=[cosNegPhiCenter -sinNegPhiCenter 0;sinNegPhiCenter cosNegPhiCenter 0;0 0 1];
%     baseVec3rot=RotZ*baseVec3;
%     %After the rotation the y axis is always orthogonal to baseVec3rot
%     baseVec2rot=[0;1;0];
%     %The third base vector is the cross product of the two already found
%     baseVec1rot=cross(baseVec3rot,baseVec2rot);
%     %Rotate back the base vectors
%     RotZback=[cosNegPhiCenter sinNegPhiCenter 0;-sinNegPhiCenter cosNegPhiCenter 0;0 0 1];
%     baseVec1=RotZback*baseVec1rot;
%     baseVec2=RotZback*baseVec2rot;
%     %get cartesian coordinates in normal "sphere" coordiante system
%     startMy (:,i)= baseVec1*x_patch + baseVec2*y_patch + baseVec3*z_patch;
%     
    startMy(1,i)=cos(phiCenter(i))*sin(thetaCenter(i));
    startMy(2,i)=sin(phiCenter(i))*sin(thetaCenter(i));
    startMy(3,i)=cos(thetaCenter(i));
end

coordinates = startMy;
          %Rotate random directed angle around this axis
          % See math.kennesaw.ed/~plaval/math4490/rotgen.pdf for derivation of forrmula for rotation around

%ChooseVector=rand(numTrajs,numel(timeSteps))<0.5;
ChooseVectortemp=rand(1,numel(timeSteps))<0.5;
ChooseVector=repmat(ChooseVectortemp,[numTrajs 1]);%%All trajectories does the same rotation trajectory in this version of the simulation algorithm
thetaStep=sqrt(2*Dr.*dt.*dim); 			%DImensionality comes in in the measn squared angular displacement
RotPath=thetaStep*ChooseVector-thetaStep*not(ChooseVector);
thetaVector=cumsum(RotPath,2);
S=sin(thetaVector);
C=cos(thetaVector);
t=1-C;



coordinates= repmat(reshape(C,1,numTrajs,size(C,2)),3,1).*repmat(coordinates,1,1,size(C,2)) +...                                                 % C*[1 0 0; 0 1 0; 0 0 1] 
             repmat(reshape(t,1,numTrajs,size(t,2)),3,1).*repmat(rotationAxes*rotationAxes'*coordinates,1,1,size(C,2)) +...                    % t*[xA yA zA]'* [xA yA zA]         
             repmat(reshape(S,1,numTrajs,size(S,2)),3,1).*repmat([0 -rotationAxes(3) rotationAxes(2); rotationAxes(3) 0 -rotationAxes(1); -rotationAxes(2) rotationAxes(1) 0]*coordinates,1,1,size(C,2)); %  S*[0  -zA  yA; zA  0 -xA; -yA  xA  0];

rho=sqrt(coordinates(1,:,:).^2+coordinates(2,:,:).^2);

my=repmat(ExScale,1,numTrajs,size(C,2)).*coordinates;
%my(3,:,:)=0;
dipoleMoment=my;%.*repmat(rho,3,1,1); % In camera experiment equal excitation in z irection, therofr no scaling of dipolemoment

Ih=conj(dipoleMoment(1,:,:)).*sum(repmat(transpose(GsquaredX(1,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(2,:,:)).*sum(repmat(transpose(GsquaredX(2,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(3,:,:)).*sum(repmat(transpose(GsquaredX(3,:)),1,numTrajs,size(C,2)).*dipoleMoment,1);
IhT(:,:)=Ih(1,:,:);
Ih=IhT;

Iv=conj(dipoleMoment(1,:,:)).*sum(repmat(transpose(GsquaredY(1,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(2,:,:)).*sum(repmat(transpose(GsquaredY(2,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(3,:,:)).*sum(repmat(transpose(GsquaredY(3,:)),1,numTrajs,size(C,2)).*dipoleMoment,1);
IvT(:,:)=Iv(1,:,:);
Iv=IvT;

I45plus=conj(dipoleMoment(1,:,:)).*sum(repmat(transpose(GsquaredX45(1,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(2,:,:)).*sum(repmat(transpose(GsquaredX45(2,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(3,:,:)).*sum(repmat(transpose(GsquaredX45(3,:)),1,numTrajs,size(C,2)).*dipoleMoment,1);
I45plusT(:,:)=I45plus(1,:,:);
I45plus=I45plusT;


I45minus=conj(dipoleMoment(1,:,:)).*sum(repmat(transpose(GsquaredY45(1,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(2,:,:)).*sum(repmat(transpose(GsquaredY45(2,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(3,:,:)).*sum(repmat(transpose(GsquaredY45(3,:)),1,numTrajs,size(C,2)).*dipoleMoment,1);
I45minusT(:,:)=I45minus(1,:,:);
I45minus=I45minusT;

% if numTrajs==1
%     Ih=Ih';
%     Iv=Iv';
%     I45plus=I45plus';
%     I45minus=I45minus';
%     
% end


Isum=Ih+Iv;
hProb=Ih./Isum;
plus45Prob=I45plus./Isum;
%vProb=Iv./Isum;
%P=(Ih-Iv)./Isum;

% remove values that are part of the burnin time

Ih=Ih(:,burninIndex:end);
Iv=Iv(:,burninIndex:end);

I45plus=I45plus(:,burninIndex:end);
I45minus=I45minus(:,burninIndex:end);


%Do binning
if binning >1
    Ihtemp=reshape(Ih(:,1:end-1),[numTrajs,binning,round(length(Ih)/binning)]);
    Ivtemp=reshape(Iv(:,1:end-1),[numTrajs,binning,round(length(Iv)/binning)]);
    
    Ihtemp2=sum(Ihtemp,2);
    Ivtemp2=sum(Ivtemp,2);
    
    
    IhExact=reshape(Ihtemp2,[numTrajs,round(length(Ih)/binning)]);
    IvExact=reshape(Ivtemp2,[numTrajs,round(length(Iv)/binning)]);
    
    I45plustemp=reshape(I45plus(:,1:end-1),[numTrajs,binning,round(length(I45plus)/binning)]);
    I45minustemp=reshape(I45minus(:,1:end-1),[numTrajs,binning,round(length(I45minus)/binning)]);
    
    I45plustemp2=sum(I45plustemp,2);
    I45minustemp2=sum(I45minustemp,2);
    
    
    I45plusExact=reshape(I45plustemp2,[numTrajs,round(length(I45plus)/binning)]);
    I45minusExact=reshape(I45minustemp2,[numTrajs,round(length(I45minus)/binning)]);
    
    
else
    IhExact=Ih;
    IvExact=Iv;
    
    I45plusExact=I45plus;
    I45minusExact=I45minus;
end
%PhvExact=P;




if strcmp(StringCom,'AddShutnoise')
    avgPhotons=sum(conj(dipoleMoment).*dipoleMoment,1)*averagePhotonsinStep;
    avgPhotonsT(:,:)=avgPhotons(1,:,:);
%     if numTrajs==1
%         avgPhotonsT=avgPhotonsT';
%     end
    numPhotons=poissrnd(avgPhotonsT); %The number of photons in this fixed time step is pooison distributed
    
    
    IhD=binornd(numPhotons,hProb);% The number of photons on the horizotal channel in this time  step is binomialy distributed with number of trials numPhotons and probability of sucess hProb
    IvD=numPhotons - IhD;
    
    I45plusD=binornd(numPhotons,plus45Prob);
    I45minusD=numPhotons-I45plusD;
    
    IhD=IhD(:,burninIndex:end);
    IvD=IvD(:,burninIndex:end);
    
    I45plusD=I45plusD(:,burninIndex:end);
    I45minusD=I45minusD(:,burninIndex:end);
    
    if binning>1
        Ihtemp=reshape(IhD(:,1:end-1),[numTrajs,binning,round(length(IhD)/binning)]);
        Ivtemp=reshape(IvD(:,1:end-1),[numTrajs,binning,round(length(IvD)/binning)]);
        
        Ihtemp2=sum(Ihtemp,2);
        Ivtemp2=sum(Ivtemp,2);
        IhDiscrete=reshape(Ihtemp2,[numTrajs,round(length(Ih)/binning)]);
        IvDiscrete=reshape(Ivtemp2,[numTrajs,round(length(Iv)/binning)]);
        
        I45plustemp=reshape(I45plusD(:,1:end-1),[numTrajs,binning,round(length(I45plusD)/binning)]);
        I45minustemp=reshape(I45minusD(:,1:end-1),[numTrajs,binning,round(length(I45minusD)/binning)]);
        
        I45plustemp2=sum(I45plustemp,2);
        I45minustemp2=sum(I45minustemp,2);
        
        
        I45plusExact=reshape(I45plustemp2,[numTrajs,round(length(I45plus)/binning)]);
        I45minusExact=reshape(I45minustemp2,[numTrajs,round(length(I45minus)/binning)]);
    
        
    else
        IhDiscrete=IhD;
        IvDiscrete=IvD;
        
        I45plusDiscrete=I45plusD;
        I45minusDiscrete=I45minusD;
    end
    
    
else
    IhDiscrete=[];
    IvDiscrete=[];
    
    I45plusDiscrete=[];
    I45minusDiscrete=[];
    
end

%%%%%
%PHVExact=(IhExact-IvExact)./(IhExact+IvExact);
%PHVDiscrete=(IhDiscrete-IvDiscrete)./(IhDiscrete+IvDiscrete);

%P45Exact=(I45plusExact-I45minusExact)./(I45plusExact+I45minusExact);
%P45Discrete=(I45plusDiscrete-I45minusDiscrete)./(I45plusDiscrete-I45minusDiscrete);

end

