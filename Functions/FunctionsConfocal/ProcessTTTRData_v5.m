function Trajobject=ProcessTTTRData_v5(Scanobjects,Trajobject,Ind, AutoCorrelationParam,TimeDiffParam,p)

RCh1=[];         RCh2=[];
ACfftsignalCh1=[]; ACfftsignalCh2=[];
Ch1Pos=[];       Ch2Pos=[];


for TrajToRetrive=Ind;
for Index=[1:size(Trajobject.INFO.Traj,1)]     
    if not(isempty(strfind(Trajobject.INFO.Traj{Index},'Scan file index')))
       FileIndexOfInterest=str2double(Trajobject.INFO.Traj{Index}(1));
    end
end
%[A B]=ind2sub(size(Trajobject.INFO.Traj'),strfind(Trajobject.INFO.Traj','Scan file index')-2);
%FileIndexOfInterest=str2num(Trajobject.INFO.Traj'(1:A,B)');

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

TTTRCh1Time{TrajToRetrive}=Time(Ch1Pos);
TTTRCh2Time{TrajToRetrive}=Time(Ch2Pos);

TTTRStartTime{TrajToRetrive}=Time(TTTRStart);
TTTRStopTime{TrajToRetrive}=Time(TTTRStop);

Ch1PosCounts(TrajToRetrive)=size(Ch1Pos,1);
Ch2PosCounts(TrajToRetrive)=size(Ch2Pos,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Remove time between measurements
StartstopTimeDiff=[0; TTTRStartTime{TrajToRetrive}(2:end)-TTTRStopTime{TrajToRetrive}(1:end-1)];

CTime=zeros(size(Time));
CTime(TTTRStart)=StartstopTimeDiff;
CompTime=cumsum(CTime);
TimeNoInbetweens=Time-CompTime;

if p==1
TTTRCh1Time{TrajToRetrive}=TimeNoInbetweens(Ch1Pos);
TTTRCh2Time{TrajToRetrive}=TimeNoInbetweens(Ch2Pos);
end

%%% autocorrelation
%Binning=AutoCorrelationParam(1);
Res=AutoCorrelationParam(1);
points=AutoCorrelationParam(2);
%[RCh1(TrajToRetrive,:) RCh2(TrajToRetrive,:) k(TrajToRetrive,:)]=BinaryAoutoCorr_V5(DataSelectionRules,TTTRCh1Time{TrajToRetrive},TTTRCh2Time{TrajToRetrive},Res,points);
[RCh1(TrajToRetrive,:) RCh2(TrajToRetrive,:) RChSum(TrajToRetrive,:) RCrosscorr(TrajToRetrive,:) k(TrajToRetrive,:)]=BinaryAoutoCorr_V6(DataSelectionRules,TTTRCh1Time{TrajToRetrive},TTTRCh2Time{TrajToRetrive},Res,points);

%System inpulse 
Res=AutoCorrelationParam(1);
points=AutoCorrelationParam(2);
MinTotTime=DataSelectionRules.Trajectory.MinTrajectoryTimeLength;
A=min(TTTRCh1Time{TrajToRetrive}(1),TTTRCh2Time{TrajToRetrive}(1));
B=max(TTTRCh1Time{TrajToRetrive}(end),TTTRCh2Time{TrajToRetrive}(end));
Minpoints=floor((B-A)/Res);       %Use Actual trajectory length
pointsTraj=Minpoints;

LinearTime=A+Res*[0:pointsTraj];
omega=2*pi/3e-3;
TTTRSystemTimeSin=(1+sin(LinearTime*omega))/2;

[Rxcorr Lag]=xcov(TTTRSystemTimeSin,points,'none');
RxcorrSin(TrajToRetrive,:)=(Rxcorr(points+1:end))';
kSin(TrajToRetrive,:)=(Res*Lag(points+1:end))';

% FFT of autocorrelation
ACfftsignalCh1(TrajToRetrive,:)=fftshift(fft(RCh1(TrajToRetrive,:)));
ACfftsignalCh2(TrajToRetrive,:)=fftshift(fft(RCh2(TrajToRetrive,:)));


MinTotTime=DataSelectionRules.Trajectory.MinTrajectoryTimeLength;
N=min(points,floor(MinTotTime/Res));

f_s=1/Res;
ACF(TrajToRetrive,:)=f_s/N*[-N/2:N/2];

%%Save data in the Trajobject 

%  %TTTR photon pointer to time events 
%  Trajobject.TTTTRTrajs.TimePointers.Ch1Pos{TrajToRetrive}=Ch1Pos;
%  Trajobject.TTTTRTrajs.TimePointers.Ch2Pos{TrajToRetrive}=Ch2Pos;
Trajobject.TTTTRTrajs.TimePointers.TTTRCh1Time=TTTRCh1Time;
Trajobject.TTTTRTrajs.TimePointers.TTTRCh2Time=TTTRCh2Time;
%  Trajobject.TTTTRTrajs.TimePointers.TTTRStartTime=TTTRStartTime;
%  Trajobject.TTTTRTrajs.TimePointers.TTTRStopTime=TTTRStopTime;

% Basic stat
Trajobject.TTTTRTrajs.BasicStat.Ch1PosCounts=Ch1PosCounts;
Trajobject.TTTTRTrajs.BasicStat.Ch2PosCounts=Ch2PosCounts;

%autocorrelation
Timechift=k;
Trajobject.TTTTRTrajs.Autocorrelation.Timechift=Timechift;
Trajobject.TTTTRTrajs.Autocorrelation.RCh1=RCh1;
Trajobject.TTTTRTrajs.Autocorrelation.RCh2=RCh2;
Trajobject.TTTTRTrajs.Autocorrelation.RChSum=RChSum;
Trajobject.TTTTRTrajs.Autocorrelation.RCrosscorr=RCrosscorr;

Trajobject.TTTTRTrajs.Autocorrelation.ACF=ACF;
Trajobject.TTTTRTrajs.Autocorrelation.ACfftsignalCh1=ACfftsignalCh1;
Trajobject.TTTTRTrajs.Autocorrelation.ACfftsignalCh2=ACfftsignalCh2;

%System autocorrelation 
Timechift=kSin;
Trajobject.TTTTRTrajs.Autocorrelation.System.Timechift=Timechift;
Trajobject.TTTTRTrajs.Autocorrelation.System.RxcorrSin=RxcorrSin;

end
%  
