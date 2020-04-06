function [C, LagTime]=PointerBasedCorrx_V5(TimeDataX,TimeDataY,B,ncas,extraTLags,NormFlag)
%Calculates autocorrelation using photon arival times. 
%The algorithme is based on: 
%Wahl Michael, Gregor Ingo, Patting Matthias, Enderlein J??rg
%2003 Vol. 11, No. 26 OPTICS EXPRESS 3583

%TimeData: Vector with photon arival times in units of the sampleling clock 
%B: Amount of lag times in a linjear section 
%ncas: Amount of lag timme sections with increases by a factor 2
%NormFlag can be: None, Biased, UnBiased, Coff 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Testing
%clear all
%close all
%Test Data Pointer Corr Data
%TimeDataX=sort(round(10^9*rand(1,10000)));%In units of 5ns
%TimeDataY=sort(round(10^9*rand(1,20000)));%In units of 5ns
%Duplicates=find(diff(TimeData)==0);
%TimeDataX=sort([TimeDataX 4+TimeDataX(1:2000:end) 8+TimeDataX(1:2200:end)]);
%TimeDataY=TimeDataX;

%B=2;
%ncas=10;
Offset=2^(0);
Scale=1;
%NormFlag='None';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NTot=B*ncas+1;
j=[1:NTot-1];
Rescale=Offset*[2.^(Scale*floor((j-1)/B))];
[~,ind1]=min(abs(Rescale-extraTLags(1)));
[~,ind2]=min(abs(Rescale-extraTLags(end)));
Rescale=[Rescale(1:ind1-1) diff(extraTLags) Rescale(ind2+1:end)];
Rescale=[1 Rescale]; %added by Emil 2018_02_26
NTot=numel(Rescale)+1; %added by Emil 2018_02_26
LagTime=[0 cumsum(Rescale)];%lagtime in units of the data
Rescale=[1 Rescale];

LX=TimeDataX(end);
SX=size(TimeDataX,2);
TimeDataRepX=repmat(TimeDataX,NTot,1);
RescaleRepX=repmat(Rescale',1,SX);
TimeDataRepX=sort(floor(TimeDataRepX./RescaleRepX),2); %Rescale data

LY=TimeDataY(end);
SY=size(TimeDataY,2);
TimeDataWLagY=repmat(TimeDataY,NTot,1) + repmat(LagTime',1,SY);
RescaleRepY=repmat(Rescale',1,SY);
TimeDataWLagY=sort(floor(TimeDataWLagY./RescaleRepY),2); %Rescale data

N=max(LX,LY);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiple rows
%a=[floor([1 1   2   2   3   5   5   5   7   8   8 20 30 35 35 35 35 50 50]/1);
%   floor([1 1   2   2   3   5   5   5   7   8   8 20 30 35 35 35 35 50 50]/4);
%   floor([1 1   2   2   3   5   5   5   7   8   8 20 30 35 35 35 35 50 50]/8)];
%
%b=a+[1 ;4; 8];%k=1 4 8
X=TimeDataRepX;
Y=TimeDataWLagY;

Z=zeros(NTot,1);
O=ones(NTot,1);

DiffX=[(diff(X,1,2)>0) O];
WX=[Z diff(sort(DiffX.*repmat([1:SX],NTot,1),2),1,2)];
YX=sort(DiffX.*X,2);
YX(WX==0)=-1;

DiffY=[(diff(Y,1,2)>0) O];
WY=[Z diff(sort(DiffY.*repmat([1:SY],NTot,1),2),1,2)];
YY=sort(DiffY.*Y,2);
YY(WY==0)=-1;

[A B]=sort([YX YY],2);
w=[WX WY];
IndexRep=repmat([1:NTot]',1,SX+SY)';
BT=B';
W=reshape(w(sub2ind([NTot SX+SY],IndexRep(:), BT(:))),SX+SY,NTot)';

%W1=[diff(A,1,2)==0].*W(:,1:end-1);
%W2=[diff(A,1,2)==0].*W(:,2:end);
%WP=W1.*W2;
WP=[diff(A,1,2)==0].*W(:,1:end-1).*W(:,2:end);
Corr=sum(WP,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CorrData=Corr'./Rescale;
LagTime;

switch NormFlag
case {'None'}
  C=CorrData;
case {'Biased'}
  C=CorrData/N;
case {'UnBiased'}
  C=CorrData./(N-LagTime);
 case {'Coff'}
  C=CorrData./CorrData(1);
end


%%%%
%%
% close all
% Res=5e-9;
% figure()
% semilogx(LagTime(2:end)*Res*1e6,C(2:end))
% xlabel('us')
% 
% figure()
% semilogx(LagTime(2:end)*Res*1e6,Rescale(2:end)*Res*1e6,'r')
% xlabel('us')
% ylabel('us')
% ['Largest time:' num2str(LagTime(end)*Res)]