function [IhExact IvExact IhDiscrete IvDiscrete]=PhotonTraceRotSIm_binning_newinformat(Param,Dr,dt,trajLength)

%  Param.rotationAxes=[1;0;0];
%  Param.Dr=[5000];
%  Param.dim=1;
%  Param.dt=1e-6;
%  Param.numTrajs=20;
%  Param.trajLength=0.2;
%  Param.averagePhotonsinStep=40000;
%  Param.ExScale=[1;1;0];
%  
%  Param.GsquaredX=GsquaredX;
%  Param.GsquaredY=GsquaredY;
%  
%  Param.StringCom='NoShutnoise';%'AddShutnoise';

%% Simulate polarization intensities, signals and autocrorreltation for rotational diffusion with photon shot noise
StringCom=Param.StringCom;
GsquaredX=Param.GsquaredX;
GsquaredY=Param.GsquaredY;

rotationAxes =Param.rotationAxes;% [1;0;0]; % Every column is an axis around wich rotation is performed 
%Dr = Param.Dr;%[5000]; %One dimensiaonal rotatioal diffusion constants for the different axes. In radians^2/s
%dt=Param.dt;  %1e-6; %in seconds
dim=Param.dim;%1;
ExScale=Param.ExScale;%scaling coordinate for to include the exitation laser behavior  
averagePhotons=Param.averagePhotonsinStep;%40000;%Average number of photons per second. Number of photons in a timestep is
%assumed to be  poission distribtued.
binning=Param.binning ; %Number of timesteps to bin into one bin (sum intesieties), done after the generatio of the intensieties traces
numTrajs=Param.numTrajs;%20; % number of trajectoreis to simulate, much larger in real data
%trajLength=Param.trajLength;%0.2; %Length of trajectories in seconds
burninTime=Param.burninTime;
startPhi=Param.startPhi;
startTheta=Param.startTheta;
burninIndex=trajLength/burninTime;

rotationForces = [0]; % harmonic potential force constant *k), to get acting force 2*k(x-x_0)
numAxes=numel(rotationForces);

timeSteps=0:dt:trajLength;

averagePhotonsinStep=averagePhotons*dt;

IhDiscrete=[];%cell(1,numTrajs);
IvDiscrete=[];%cell(1,numTrajs);

IhExact=[];%cell(1,numTrajs);
IvExact=[];%cell(1,numTrajs);
PhvExact=[];%cell(1,numTrajs);

  
IhExact=zeros(numTrajs,numel(timeSteps));
IvExact=zeros(numTrajs,numel(timeSteps));
PhvExact=zeros(numTrajs,numel(timeSteps));
phi = ones(1,numTrajs)*startPhi;
theta = ones(1,numTrajs)*startTheta;

coordinates = myEqualExciation(phi,theta);
          %Rotate random directed angle around this axis
          % See math.kennesaw.ed/~plaval/math4490/rotgen.pdf for derivation of forrmula for rotation around

ChooseVector=rand(numTrajs,numel(timeSteps))<0.5;
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
dipoleMoment=my.*repmat(rho,3,1,1); % Because of no Z exctation dipolemoment is my scaled with rho = sin(theta)

Ih=conj(dipoleMoment(1,:,:)).*sum(repmat(transpose(GsquaredX(1,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(2,:,:)).*sum(repmat(transpose(GsquaredX(2,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(3,:,:)).*sum(repmat(transpose(GsquaredX(3,:)),1,numTrajs,size(C,2)).*dipoleMoment,1);
IhT(:,:)=Ih(1,:,:); Ih=IhT;

Iv=conj(dipoleMoment(1,:,:)).*sum(repmat(transpose(GsquaredY(1,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(2,:,:)).*sum(repmat(transpose(GsquaredY(2,:)),1,numTrajs,size(C,2)).*dipoleMoment,1) +...
   conj(dipoleMoment(3,:,:)).*sum(repmat(transpose(GsquaredY(3,:)),1,numTrajs,size(C,2)).*dipoleMoment,1);
IvT(:,:)=Iv(1,:,:); Iv=IvT;

Isum=Ih+Iv;
hProb=Ih./Isum;
%vProb=Iv./Isum;
%P=(Ih-Iv)./Isum;

% remove values that are part of the burnin time

Ih=Ih(:,burninIndex:end);
Iv=Iv(:,burninIndex:end);

%Do binning
if binning >1
    Ihtemp=reshape(Ih(:,1:end-1),[numTrajs,binning,round(length(Ih)/binning)]);
    Ivtemp=reshape(Iv(:,1:end-1),[numTrajs,binning,round(length(Iv)/binning)]);
    
    Ihtemp2=sum(Ihtemp,2);
    Ivtemp2=sum(Ivtemp,2);
    
    
    IhExact=reshape(Ihtemp2,[numTrajs,round(length(Ih)/binning)]);
    IvExact=reshape(Ivtemp2,[numTrajs,round(length(Iv)/binning)]);
else
    IhExact=Ih;
    IvExact=Iv;
end
%PhvExact=P;




if strcmp(StringCom,'AddShutnoise')
    IhDiscrete=zeros(numTrajs,numel(timeSteps));
    IvDiscrete=zeros(numTrajs,numel(timeSteps));  
    avgPhotons=sum(conj(dipoleMoment).*dipoleMoment,1)*averagePhotonsinStep;
    avgPhotonsT(:,:)=avgPhotons(1,:,:);
    numPhotons=poissrnd(avgPhotonsT); %The number of photons in this fixed time step is pooison distributed
    
    
    IhD=binornd(numPhotons,hProb);% The number of photons on the horizotal channel in this time  step is binomialy distributed with number of trials numPhotons and probability of sucess hProb
    IvD=numPhotons - IhD;
    
    IhD=IhD(:,burninIndex:end);
    IvD=IvD(:,burninIndex:end);
    if binning>1
        Ihtemp=reshape(IhD(:,1:end-1),[numTrajs,binning,round(length(IhD)/binning)]);
        Ivtemp=reshape(IvD(:,1:end-1),[numTrajs,binning,round(length(IvD)/binning)]);
        
        Ihtemp2=sum(Ihtemp,2);
        Ivtemp2=sum(Ivtemp,2);
        IhDiscrete=reshape(Ihtemp2,[numTrajs,round(length(Ih)/binning)]);
        IvDiscrete=reshape(Ivtemp2,[numTrajs,round(length(Iv)/binning)]);
    else
        IhDiscrete=IhD;
        IvDiscrete=IvD;
    end
    
    
else
    IhDiscrete=[];
    IvDiscrete=[];
end


