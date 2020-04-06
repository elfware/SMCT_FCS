function [RCh1 RCh2 k]=BinaryAoutoCorr_V5(DataSelectionRules,TTTRCh1Time,TTTRCh2Time,Res,points)
%  clear all
%  close all
%  TTTRPos=sort(round(rand(1e5,1)*1e6))*5e-9;
%  points=1000;
%  m=1;
%Res=4e-3;

MinTotTime=DataSelectionRules.Trajectory.MinTrajectoryTimeLength;
A=min(TTTRCh1Time(1),TTTRCh2Time(1));
B=max(TTTRCh1Time(end),TTTRCh2Time(end));

Minpoints=floor((B-A)/Res);       %Use Actual trajectory length
%Minpoints=floor(MinTotTime/Res); Use chortest trajectory length 
pointsTraj=Minpoints;

LinearTime=A+Res*[0:pointsTraj];

TTTRCh1LinEvents=hist(TTTRCh1Time,LinearTime);
TTTRCh2LinEvents=hist(TTTRCh2Time,LinearTime);

TTTRCh1LinEvents=TTTRCh1LinEvents(1:end-1); %%Remove last element where all higher counts are placed
TTTRCh2LinEvents=TTTRCh2LinEvents(1:end-1); %%Remove last element where all higher counts are placed

%%Autocorrelation 
maxpoints=points;%;points-1;
[Rxcorr Lag]=xcov(TTTRCh1LinEvents,maxpoints,'none');%unbiased
%RCh1=(Rxcorr(maxpoints+1:end)/Rxcorr(maxpoints+1))';
RCh1=(Rxcorr(maxpoints+1:end))';

[Rxcorr Lag]=xcov(TTTRCh2LinEvents,maxpoints,'none');
%RCh2=(Rxcorr(maxpoints+1:end)/Rxcorr(maxpoints+1))';
RCh2=(Rxcorr(maxpoints+1:end))';

k=(Res*Lag(maxpoints+1:end))';
%%%%%%%%%%%%%%%%%%%%%
