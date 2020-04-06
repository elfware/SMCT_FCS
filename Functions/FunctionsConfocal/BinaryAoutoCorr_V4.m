function [R k]=BinaryAoutoCorr_V4(TTTRPos,points,m)
%  clear all
%  close all
%  TTTRPos=sort(round(rand(1e5,1)*1e6))*5e-9;
%  points=1000;
%  m=1;

%Binary_pointer=round((TTTRPos-TTTRPos(1))/5e-9)+1;
%%xcorr test
%tic

%photonDens=1/diff((Binary_pointer-1)*5e-9);
%photonDensTime=(Binary_pointer(1:end-1)-1)*5e-9;

photonDens=1/diff(TTTRPos-TTTRPos(1));
photonDensTime=TTTRPos(1:end-1)-TTTRPos(1);

maxpoints=min(points, size(photonDensTime,1)-1);
[Rxcorr Lag]=xcorr(photonDens,maxpoints,'biased');
k=(photonDensTime(Lag(maxpoints+1:end)+1))';
R=(Rxcorr(maxpoints+1:end)/Rxcorr(maxpoints+1))';
%toc
%figure()
%plot(photonDensTime(Lag(points+1:end)+1),Rxcorr(points+1:end)/Rxcorr(points+1))
%%
