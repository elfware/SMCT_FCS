Param=[1.45 1.51 2000 pi/2 0];

numericalAperture=Param(1);%1.45;%1.45
refractiveIndex=Param(2);%1.51;%1.4811;%of oil
thetaMax=asin(numericalAperture/refractiveIndex);

xSteps=Param(3);%2000; %number of steps in x in the backfocal plane (in the integration)
ySteps=xSteps;

thextemp=linspace(-1,1,xSteps);
theytemp=linspace(-1,1,ySteps);
thex=thextemp(1:end-1)+(thextemp(2)-thextemp(1))/2; %To get integration according to trapzoid method
they=theytemp(1:end-1)+(theytemp(2)-theytemp(1))/2;
dxdy=(thex(2)-thex(1))*(they(2)-they(1));
numx=numel(they);
numy=numel(thex);

Gtot=zeros(3);
EX=zeros(numy,numx); 
EY=zeros(numy,numx);
Etot=zeros(numy,numx);

GsquaredX=zeros(3);%Integrated G^T*G components for the x component in the intensity in the image plane
GsquaredY=zeros(3);
Gtot=zeros(3);
GsquaredXys=cell(1,numy); 
GsquaredYys=cell(1,numy);
GtotTys=cell(1,numy);

x = thex;
y = they;
[xx,yy]=meshgrid(x,y);
xxL=transpose(xx(:));
yyL=transpose(yy(:));
rho=sqrt(xxL.^2+yyL.^2);  %% this has values larger then 1  figure();imagesc(reshape(rho,size(xx)));colorbar
%rho(rho>sin(thetaMax))=inf; %% Set all values larger then aloud by max Na to 0
theta=asin(rho);          %this is imagenary if to larg values are used in the above figure();imagesc(reshape(theta,size(xx)));colorbar
phi=atan2(yyL,xxL);       %% figure();imagesc(reshape(phi,size(xx)));colorbar
Gbfp=greenTensorBFP_v2(phi,theta);
%GbfpTest=reshape(Gbfp(1,1,:),xSteps-1,ySteps-1,1).*(rho>sin(thetaMax));
%imagesc(Gbfp)

%imagesc(reshape(rho,1999,1999));     
%imagesc(reshape(phi,1999,1999));

%Create Greentensor part for each polarization       
GbfpXT=[];
GbfpX=Gbfp(1,:,:);  
GbfpXT(:,1,:)=conj(GbfpX);

factor=repmat(reshape(rho<sin(thetaMax),[1 1 numel(rho)]),3);%% Set all values larger then aloud by max Na to 0
GsquaredX=sum(repmat(GbfpXT,1,3).*repmat(GbfpX,3,1)*dxdy.*factor,3); %This takes Gbfp(1,:)'*Gbfp(1,:) and sums over all theta and phi

GbfpYT=[];
GbfpY=Gbfp(2,:,:);  
GbfpYT(:,1,:)=conj(GbfpY);
GsquaredY=sum(repmat(GbfpYT,1,3).*repmat(GbfpY,3,1)*dxdy.*factor,3);


%%Theory
NA=Param(1);
n_oil=Param(2);
rhoM=NA/n_oil;

h11=(9*pi/8)*asin(rhoM) + (3*pi/8)*rhoM.*sqrt(1-rhoM.^2) + (pi/2)*rhoM;
h22=(3*pi/8)*asin(rhoM) + (pi/8)*rhoM.*sqrt(1-rhoM.^2)   - (pi/2)*rhoM;
h33=(pi/2)*asin(rhoM)   - (pi/2)*rhoM.*sqrt(1-rhoM.^2);

GsquaredXTheory=[h11 0  0;
                  0 h22 0;
                  0  0  h33];
                  

[diag(GsquaredXTheory) diag(GsquaredX) diag(GsquaredY)]
[diag(GsquaredXTheory)./diag(GsquaredX)]