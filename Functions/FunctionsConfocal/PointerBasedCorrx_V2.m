function [C LagTime]=PointerBasedCorrx_V1(TimeData,B,ncas,NormFlag)
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
test=0;% Use for testing 
if test==1
%Test Data Pointer Corr Data
TimeData=sort(round(10^9*rand(1,10000)));%In units of 5ns
Duplicates=find(diff(TimeData)==0);

B=50;
ncas=20;
NormFlag='None';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=1:B*ncas;
Rescale=[2.^floor((j-1)/B)]; %sort(repmat([2.^(0:ncas-1)],1,B));
LagTime=cumsum(Rescale);%lagtime in units of the data

N=TimeData(end);
TimeDataRep=repmat(TimeData,size(LagTime,2),1);
TimeDataWLag=TimeDataRep + repmat(LagTime',1,size(TimeData,2));
RescaleRep=repmat(Rescale',1,size(TimeData,2));

%Rescale data
TimeDataRep=sort(floor(TimeDataRep./RescaleRep),2);
TimeDataWLag=sort(floor(TimeDataWLag./RescaleRep),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiple rows
%a=[floor([1 1   2   2   3   5   5   5   7   8   8 20 30 35 35 35 35 50 50]/1);
%   floor([1 1   2   2   3   5   5   5   7   8   8 20 30 35 35 35 35 50 50]/4);
%   floor([1 1   2   2   3   5   5   5   7   8   8 20 30 35 35 35 35 50 50]/8)];
%
%b=a+[1 ;4; 8];%k=1 4 8
a=TimeDataRep;
b=TimeDataWLag;

Z=zeros(size(a,1),1);
O=1+Z;

Wa=[Z diff(sort([(diff(a,1,2)>0) O].*repmat([1:size(a,2)],size(a,1),1),2),1,2)];
Ya=sort([(diff(a,1,2)>0) O].*a,2);
Ya(Wa==0)=-1;

Wb=[Z diff(sort([(diff(b,1,2)>0) O].*repmat([1:size(b,2)],size(b,1),1),2),1,2)];
Yb=sort([(diff(b,1,2)>0) O].*b,2);
Yb(Wb==0)=-1;

[A B]=sort([Ya Yb],2);
w=[Wa Wb];
%W=reshape(w(sub2ind(size([Ya Yb]),[repmat([1:size(B,1)]',1,size(B,2))'(:)], B'(:))),size([Ya Yb],2),size([Ya Yb],1))';
SYab=size([Ya Yb]);
SB=size(B);
SBRep=repmat([1:SB(1)]',1,SB(2))';
BT=B';
W=reshape(w(sub2ind(SYab,SBRep(:), BT(:))),SYab(2),SYab(1))';

%W1=[diff(A,1,2)==0].*W(:,1:end-1);
%W2=[diff(A,1,2)==0].*W(:,2:end);
%WP=W1.*W2;
WP=[diff(A,1,2)==0].*W(:,1:end-1).*W(:,2:end);
Corr=sum(WP,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CorrData=[size(TimeData,2) Corr'./Rescale];
LagTime=[0 LagTime];

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
if test==1
Res=5e-9;
figure()
semilogx(LagTime(2:end)*Res,C(2:end))
end
%
%%%
%figure()
%plot(LagTime(2:end)*Res,C(2:end))
%axis([1e-5 1e-2 0 1])



