function Trajobject=ProcessTTTRData_v9(Scanobjects,Trajobject,Ind, AutoCorrelationParam,p)
%This is a new function to use a logarithmic scale for calculating
%autocorrelation. 
%First Modified 2017-06-02

RCh1=[];         RCh2=[];
ACfftsignalCh1=[]; ACfftsignalCh2=[];
Ch1Pos=[];       Ch2Pos=[];
RChPol=[];

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
TTTRRes=Scanobject.ScanSpeeds.TTTRResolution;

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
B=AutoCorrelationParam(1);
ncas=AutoCorrelationParam(2);
[RCh1(TrajToRetrive,:) RCh2(TrajToRetrive,:) RChSum(TrajToRetrive,:) RChCross(TrajToRetrive,:) k(TrajToRetrive,:)]=BinaryAoutoCorr_V10(DataSelectionRules,TTTRCh1Time{TrajToRetrive},TTTRCh2Time{TrajToRetrive},B,ncas,TTTRRes);

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

%Autocorrelation
Timechift=k;
Trajobject.TTTTRTrajs.Autocorrelation.Timechift=Timechift;
Trajobject.TTTTRTrajs.Autocorrelation.RCh1=RCh1;
Trajobject.TTTTRTrajs.Autocorrelation.RCh2=RCh2;
Trajobject.TTTTRTrajs.Autocorrelation.RChSum=RChSum;
Trajobject.TTTTRTrajs.Autocorrelation.RChCross=RChCross;

end
%  
