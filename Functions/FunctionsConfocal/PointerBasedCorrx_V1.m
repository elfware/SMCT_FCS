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
%Test Data Pointer Corr Data
%TimeData=sort(round(10^9*rand(1,10000)));%In units of 5ns
%Duplicates=find(diff(TimeData)==0);

%B=100;
%ncas=20;
%NormFlag='UnBiased';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=1:B*ncas;
Rescale=[2.^floor((j-1)/B)]; %sort(repmat([2.^(0:ncas-1)],1,B));
LagTime=cumsum(Rescale);%lagtime in units of the data
%NRescale=[2.^floor(j/B)];

N=TimeData(end);
TimeDataRep=repmat(TimeData,size(LagTime,2),1);
TimeDataWLag=TimeDataRep + repmat(LagTime',1,size(TimeData,2));
RescaleRep=repmat(Rescale',1,size(TimeData,2));

%Rescale data
TimeDataRep=sort(floor(TimeDataRep./RescaleRep),2);
TimeDataWLag=sort(floor(TimeDataWLag./RescaleRep),2);

CorrDataSet=sort([TimeDataRep TimeDataWLag],2);
EqualElement=diff(CorrDataSet,1,2)==0;
CorrData=[size(TimeData,2) sum(EqualElement,2)'./Rescale];
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







