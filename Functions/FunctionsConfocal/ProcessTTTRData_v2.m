function Trajobject=ProcessTTTRData_v2(Scanobjects,Trajobject,Ind, AutoCorrelationParam,TimeDiffParam)


RCh1=[];         RCh2=[];
ACfftsignalCh1=[]; ACfftsignalCh2=[];
%  fftsignalCh1=[]; fftsignalCh2=[];
%  Ch1Bind=[];      Ch2Bind=[];
Ch1Pos=[];       Ch2Pos=[];

for TrajToRetrive=Ind;

[A B]=ind2sub(size(Trajobject.INFO.Traj'),strfind(Trajobject.INFO.Traj','Scan file index')-2);
FileIndexOfInterest=str2num(Trajobject.INFO.Traj'(1:A,B)');

FileOfInterest=Trajobject.Traj(TrajToRetrive,FileIndexOfInterest);
Scanobject=Scanobjects{FileOfInterest}.Scanobject;
DataSelectionRules=Scanobjects{FileOfInterest}.DataSelectionRules;


%%Extract Scanning Start stop for trajectory
StartTraj=Trajobject.Traj(TrajToRetrive,2);
StopTraj=Trajobject.Traj(TrajToRetrive,3);

TTTRStart=Scanobject.ScanData.TTTREventStructure.StartEvents([StartTraj:StopTraj]);
TTTRStop=Scanobject.ScanData.TTTREventStructure.StopEvents([StartTraj:StopTraj]);

Ch1Pos=Scanobject.ScanData.TTTREventStructure.APDEvents.Ch1Pos;%H
Ch2Pos=Scanobject.ScanData.TTTREventStructure.APDEvents.Ch2Pos;%V
Time=Scanobject.ScanData.TTTREventStructure.TimeEvents.Sync;

Ch1Pos=Ch1Pos(find( and(Ch1Pos >= TTTRStart(1), Ch1Pos <= TTTRStop(end)) ));
Ch2Pos=Ch2Pos(find( and(Ch2Pos >= TTTRStart(1), Ch2Pos <= TTTRStop(end)) ));

%TTTRCh1Pos=Time(Ch1Pos);
TTTRCh1Time{TrajToRetrive}=Time(Ch1Pos);
%TTTRCh2Pos=Time(Ch2Pos);
TTTRCh2Time{TrajToRetrive}=Time(Ch2Pos);

%  TTTRStartTime{TrajToRetrive}=Time(TTTRStart);
%  TTTRStopTime{TrajToRetrive}=Time(TTTRStop);

Ch1PosCounts(TrajToRetrive)=size(Ch1Pos,1);
Ch2PosCounts(TrajToRetrive)=size(Ch2Pos,1);

%TTTRTrajTimeLength=diff(sort([TTTRCh1Pos;TTTRCh2Pos])([1 end]));

%  %%Binned Trajectory building 
%  MinTotTime=DataSelectionRules.Trajectory.MinTrajectoryTimeLength;
%  Res=AutoCorrelationParam(1);
%  Minpoints=floor(MinTotTime/Res);
%  points=Minpoints;
%  
%  A=min(TTTRCh1Time{TrajToRetrive}(1),TTTRCh2Time{TrajToRetrive}(1));
%  LinearTime(TrajToRetrive,:)=A+Res*[0:points-1];
%  
%  TTTRCh1LinEvents(TrajToRetrive,:)=hist(TTTRCh1Time{TrajToRetrive},LinearTime(TrajToRetrive,:));
%  TTTRCh2LinEvents(TrajToRetrive,:)=hist(TTTRCh2Time{TrajToRetrive},LinearTime(TrajToRetrive,:));
%  
%  %%FFT of binned trajectory 
%  ACfftBinnedCh1(TrajToRetrive,:)=fftshift(fft(TTTRCh1LinEvents(TrajToRetrive,:)));
%  ACfftBinnedCh2(TrajToRetrive,:)=fftshift(fft(TTTRCh2LinEvents(TrajToRetrive,:)));
%  
%  
%  MinTotTime=DataSelectionRules.Trajectory.MinTrajectoryTimeLength;
%  N=floor(MinTotTime/Res);
%  
%  f_sBinned=1/Res;
%  ACFBinned(TrajToRetrive,:)=f_sBinned/N*[-N/2:N/2-1];

%  %% Time diff analysis
%  Bins=linspace(0,TimeDiffParam(1),TimeDiffParam(2));
%  Ch1Bind(TrajToRetrive,:)=hist(diff(TTTRCh1Pos),Bins);
%  Ch2Bind(TrajToRetrive,:)=hist(diff(TTTRCh2Pos),Bins);
%  
%  % FFT of time diff analysis
%  fftsignalCh1(TrajToRetrive,:)=fftshift(fft(Ch1Bind(TrajToRetrive,:)));
%  fftsignalCh2(TrajToRetrive,:)=fftshift(fft(Ch2Bind(TrajToRetrive,:)));
%  
%  f_s=1/mean(diff(Bins));
%  N=size(Bins,2);
%  F=f_s/N*[-N/2:N/2-1];

%%% autocorrelation
%Binning=AutoCorrelationParam(1);
Res=AutoCorrelationParam(1);
points=AutoCorrelationParam(2);
[RCh1(TrajToRetrive,:) RCh2(TrajToRetrive,:) k(TrajToRetrive,:)]=BinaryAoutoCorr_V5(DataSelectionRules,TTTRCh1Time{TrajToRetrive},TTTRCh2Time{TrajToRetrive},Res,points);

% FFT of autocorrelation
ACfftsignalCh1(TrajToRetrive,:)=fftshift(fft(RCh1(TrajToRetrive,:)));
ACfftsignalCh2(TrajToRetrive,:)=fftshift(fft(RCh2(TrajToRetrive,:)));


MinTotTime=DataSelectionRules.Trajectory.MinTrajectoryTimeLength;
N=min(points,floor(MinTotTime/Res));

f_s=1/Res;
ACF(TrajToRetrive,:)=f_s/N*[-N/2:N/2];

%%Save data in the Trajobject 

%  %TTTR photon pointer to time events 
%  %Trajobject.TTTTRTrajs.TimePointers.Ch1Pos{TrajToRetrive}=Ch1Pos;
%  %Trajobject.TTTTRTrajs.TimePointers.Ch2Pos{TrajToRetrive}=Ch2Pos;
%  Trajobject.TTTTRTrajs.TimePointers.TTTRCh1Time=TTTRCh1Time;
%  Trajobject.TTTTRTrajs.TimePointers.TTTRCh2Time=TTTRCh2Time;
%  Trajobject.TTTTRTrajs.TimePointers.TTTRStartTime=TTTRStartTime;
%  Trajobject.TTTTRTrajs.TimePointers.TTTRStopTime=TTTRStopTime;

% Basic stat
Trajobject.TTTTRTrajs.BasicStat.Ch1PosCounts=Ch1PosCounts;
Trajobject.TTTTRTrajs.BasicStat.Ch2PosCounts=Ch2PosCounts;
%
%  %Binned trajectory
%  Trajobject.TTTTRTrajs.BinnedTraj.LinearTime=LinearTime;
%  Trajobject.TTTTRTrajs.BinnedTraj.TTTRCh1LinEvents=TTTRCh1LinEvents;
%  Trajobject.TTTTRTrajs.BinnedTraj.TTTRCh2LinEvents=TTTRCh2LinEvents;
%  
%  Trajobject.TTTTRTrajs.BinnedTraj.ACFBinned=ACFBinned;
%  Trajobject.TTTTRTrajs.BinnedTraj.ACfftBinnedCh1=ACfftBinnedCh1;
%  Trajobject.TTTTRTrajs.BinnedTraj.ACfftBinnedCh2=ACfftBinnedCh2;

%  %TimeDiff
%  Trajobject.TTTTRTrajs.TimeDiff.Bins=Bins;
%  Trajobject.TTTTRTrajs.TimeDiff.Ch1Bind=Ch1Bind;
%  Trajobject.TTTTRTrajs.TimeDiff.Ch2Bind=Ch2Bind;
%  
%  Trajobject.TTTTRTrajs.TimeDiff.F=F;
%  Trajobject.TTTTRTrajs.TimeDiff.fftsignalCh1=fftsignalCh1;
%  Trajobject.TTTTRTrajs.TimeDiff.fftsignalCh2=fftsignalCh2;

%autocorrelation
Timechift=k;
Trajobject.TTTTRTrajs.Autocorrelation.Timechift=Timechift;
Trajobject.TTTTRTrajs.Autocorrelation.RCh1=RCh1;
Trajobject.TTTTRTrajs.Autocorrelation.RCh2=RCh2;

Trajobject.TTTTRTrajs.Autocorrelation.ACF=ACF;
Trajobject.TTTTRTrajs.Autocorrelation.ACfftsignalCh1=ACfftsignalCh1;
Trajobject.TTTTRTrajs.Autocorrelation.ACfftsignalCh2=ACfftsignalCh2;

end

%  figure()
%  plot(Bins*1e6,Ch1Bind,'-b')
%  hold on
%  plot(Bins*1e6,Ch2Bind,'-r')
%  xlabel('DTime [us]')
%  ylabel('#N')
%  axis([0 500 0 max(Ch1Bind(1),Ch1Bind(2)*1.2)])
%  
%  figure()
%  plot(Bins,(Ch1Bind-Ch2Bind)./(Ch1Bind+Ch1Bind),'-b')


%  %%%
%  fftsignal=fftshift(abs(fft(Ch1Bind)));
%  f_s=1/mean(diff(Bins));
%  N=size(Bins,2);
%  F=f_s/N*[-N/2:N/2-1];


%  Index=1;
%  TrajToRetrive=slidingInd(Index);
%  DataSelectionRules.Plot.Mode='XYTimePol';%% XYTimePol, XPYTimePol, XYPTimePol, XPYPTimePol, XYPolPol, XYTimeInt
%  DataSelectionRules.Trajectory.Mode='FixWindowAv';%Use: FixWindowAv, RollingAv or DirectEstimation3Points
%  ShowDataPlots(Scanobjects,Trajobject,DataSelectionRules,TrajToRetrive,'sShowRawData')
%  %%%%%%%%
%  
%  figure()
%  plot(F,fftsignalCh1(Index,:),'r')
%  hold on
%  plot(F,fftsignalCh2(Index,:),'b')
%  grid on
%  
%  figure()
%  plot([1:points]*Res,RCh1(Index,:),'r')
%  hold on
%  plot([1:points]*Res,RCh2(Index,:),'b')
%  grid on



%  size(TTTRCh1Pos,1);
%  size(TTTRCh2Pos,1);
%  TotAPD1Scan=Trajobject.Traj(TrajToRetrive,6);
%  TotAPD2Scan=Trajobject.Traj(TrajToRetrive,7);
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Intervall=213641+[-100:100];
%  TTTRStart=Scanobject.ScanData.TTTREventStructure.StartEvents(Intervall);
%  TTTRStop=Scanobject.ScanData.TTTREventStructure.StopEvents(Intervall);
%  MTime=(Time(TTTRStop)-Time(TTTRStart))*1e6;
%  
%  
