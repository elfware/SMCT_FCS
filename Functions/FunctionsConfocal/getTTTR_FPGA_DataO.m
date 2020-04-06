function [ EventStructure ] = getTTTR_FPGA_DataO(Scanobject, ShowFlag, V)

%SEE_TTTR_FPGA_Data Summary of this function goes here
%   Takes a .bin TTTR file and extracts all events 

SeeData=1;%For debugging

%%%%%%%%% Defenisions 
pathname=Scanobject.DataPath.pathname;
filename=Scanobject.DataPath.filenameTTTRData;

BetweenMTimeScan=(1/Scanobject.ScanSpeeds.ScanSpeed)-Scanobject.ScanSpeeds.MeasurementOnTime;  %time between measurement when scaning
BetweenMTimeTrack=(1/Scanobject.ScanSpeeds.ScanSpeedTracking)-Scanobject.ScanSpeeds.MeasurementOnTime;  %time between measurement when scaning 

WRAPAROUND=2^16;
Resolution=5e-9;%secondes per tick (200Mhz clock)

%%% Load data 
fid=fopen([pathname filename]);
T3Record = fread(fid, 'uint32', 'ieee-be');
fclose(fid);

if V
  disp('File loded.')
  disp('  ')
end


if isempty(T3Record)
   EventStructure.noevent=1;
   if V
      disp('File contains no events.')
      disp('  ')
   end

else
    EventStructure.noevent=0;
    %+-------------------------------+  +-------------------------------+ 
    %|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
    %+-------------------------------+  +-------------------------------+ 
    nsync = bitand(T3Record,2^16-1);       % the lowest 16 bits:  
    %+-------------------------------+  +-------------------------------+ 
    %| | | | | | | | | | | | | | | | |  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
    %+-------------------------------+  +-------------------------------+    
    chan = bitand(bitshift(T3Record,-16),2^5-1);   % the upper 5 bits:
    %+-------------------------------+  +-------------------------------+ 
    %| | | | | | | | | | | |x|x|x|x|x|  | | | | | | | | | | | | | | | | |
    %+-------------------------------+  +-------------------------------+ 

    Ev1=find(bitand(chan,1));%event 1 occured APD 1
    Ev2=find(bitand(chan,2));%event 2 occured APD 2
    Ev4=find(bitand(chan,4));%event 4 occured Start marker
    Ev8=find(bitand(chan,8));%event 8 occured Rollover marker
    Ev16=find(bitand(chan,16));%event 16 occured Stop marker

   if V
    disp('Events has been unwraped.')
    disp('  ')
   end
end

%%% Test Data for inconsistencis
%% See if we are missing start stop values !!!
StartStopError=size(Ev4,1)~=size(Ev16,1);
if V
   if StartStopError
    disp(['Start and/or Stop events are missing/to many/notaligned, total missmatch is:' num2str(size(Ev4,1)-size(Ev16,1))])
    disp('  ')
  end
end

%% See if the starat and stop values starts in the wright order 
StartStopError=(Ev16(1)-Ev4(1))<0;
if V
  if StartStopError
    disp('Events starts with an stop marker. ')
    disp('  ')
  end
end

for LoopIndex=[1 2]

%%% Scandata:
if LoopIndex==1
  ScanCount=sum(Scanobject.ScanData.ScanMatrix(:,4:5),2);
  ReleasThresholdOK=Scanobject.ScanData.ScanMatrix(:,7);
  TrackingOnSwitch=[0;diff(ReleasThresholdOK)];
  TrackingOn=Scanobject.ScanData.ScanMatrix(:,6);
end
%%%

%%% Evaluate counts per measurement
chanNorolloverPointer=find( ((bitand(chan,8)==8).*chan)~=8 );%Points to all elemets that are not a rollover
chanNorollover=chan(chanNorolloverPointer);%all events which are not only rollover
Ev4Norollover=find(bitand(chanNorollover,4));%event 4 occured Start marker, pointer that is pointing to chanNorollover
Ev16Norollover=find(bitand(chanNorollover,16));%event 16 occured Stop marker, pointer that is pointing to chanNorollover

%[find(bitand(chanNorollover(Ev4Norollover),4)~=4) chanNorollover(Ev4Norollover)(find(bitand(chanNorollover(Ev4Norollover),4)~=4))]
%[find(bitand(chanNorollover(Ev16Norollover),16)~=16) chanNorollover(Ev16Norollover)(find(bitand(chanNorollover(Ev16Norollover),16)~=16))]


if V
if LoopIndex==1
  disp(['LoopIndex: ' num2str(LoopIndex) ' Compensating missmatch'])
  disp('  ')
elseif LoopIndex==2
  disp(['LoopIndex: ' num2str(LoopIndex) ' Rerun with new aligment'])
  disp('  ')
end
end

if and(V,SeeData)
  disp('First 10 scan and fast photon count values before doing anything:')
  [ScanCount(1:10) Ev16Norollover(1:10)-Ev4Norollover(1:10)-1]
  
  disp('')
  disp('See first tracker on switch event: ')
  TrakingPointer=find(abs(TrackingOnSwitch)==1);
  FirstOff=TrakingPointer(1);
  IntervalToSee=FirstOff+[-10:10];
  [IntervalToSee' Ev4Norollover(IntervalToSee) Ev16Norollover(IntervalToSee) ScanCount(IntervalToSee) Ev16Norollover(IntervalToSee)-Ev4Norollover(IntervalToSee)-1 TrackingOnSwitch(IntervalToSee)]
%  [3163043+[-10:20]' chanNorollover(3163043+[-10:20])]  1005115   3163043
end

%%Compensate for events that are causing missaligment 
%% The pettern is:
%%  each time the ReleasThresholdOK switches from low to high or high to low one extra start event occure
%%  each time the ReleasThresholdOK switches from  high to low on extra stop event occure
%% Thus each time ReleasThresholdOK switches from low to high a single start event needs to be removed
%% and each time ReleasThresholdOK switches from high to low a start and stop event needs to be removed
if LoopIndex==1
  %Remove first
  FindIndex=find(abs(TrackingOnSwitch)==1) + [0:sum(abs(TrackingOnSwitch)==1)-1]';
  VComp=ones(size(Ev4Norollover));
  VComp(FindIndex)=0;
  Ev4Norollover=Ev4Norollover(find(VComp));

  FindIndex=find(TrackingOnSwitch==1) + [0:sum(TrackingOnSwitch==1)-1]';
  VComp=ones(size(Ev16Norollover));
  VComp(FindIndex)=0;
  Ev16Norollover=Ev16Norollover(find(VComp));
  %%%%%%%%%%%%%%
  
  %% Find error in the start stop sequence due to removal of the start at trackingSwitchOff
  FindIndex=find(TrackingOnSwitch==-1);
  VcompStart=zeros(size(chanNorollover));
  VcompStart(Ev4Norollover(FindIndex))=1;
 
  VcompStop=zeros(size(chanNorollover));
  VcompStop(Ev16Norollover(FindIndex))=1;
  Vlogvector=cumsum(VcompStart)-cumsum(VcompStop);
  
  Vcomp=[diff(Vlogvector);0]==1 & Vlogvector<0 & bitand(chanNorollover,16)==16 ;
  
%    VlogvectorPoint=find(Vlogvector<0);
%    Vcomp=zeros(size(chanNorollover));
%    Vcomp(Ev16Norollover)=1;
%    Vcomp(VlogvectorPoint)=Vcomp(VlogvectorPoint)+1;
%    
%  find(Vcomp(Ev16Norollover)==2)(2)
%  Vcomp(Ev16Norollover)(find(Vcomp(Ev16Norollover)==2)(2)+[0:10])

%  FindIndex=find(Vcomp(Ev16Norollover)>1)(1);
  FindIndex=find(Vcomp(Ev16Norollover));
  VComp=ones(size(Ev16Norollover));
  VComp(FindIndex)=0;
  Ev16Norollover=Ev16Norollover(find(VComp));
  
end 

if and(V,SeeData)
  disp('')
  disp('After compensation:')
    
  disp('See first tracker on switch event: ')
  TrakingPointer=find(abs(TrackingOnSwitch)==1);
  FirstOff=TrakingPointer(1);
  IntervalToSee=FirstOff+[-10:10]; %[1005110:1005120];
  disp(['sum of |diff|: ' num2str(sum(abs(ScanCount(IntervalToSee)-(Ev16Norollover(IntervalToSee)-Ev4Norollover(IntervalToSee)-1))))])
  disp(['Amount of diff events: ' num2str(sum((ScanCount(IntervalToSee)-(Ev16Norollover(IntervalToSee)-Ev4Norollover(IntervalToSee)-1))~=0))])
  
  SwitchOneventsINinterval=TrakingPointer(find(and( TrakingPointer>=IntervalToSee(1),TrakingPointer<=IntervalToSee(end))));
  size(SwitchOneventsINinterval,1)
  %disp(['Amount of all events and SwitchOnevents: ' num2str(sum( abs(ScanCount(SwitchOneventsINinterval)-(Ev16Norollover(SwitchOneventsINinterval)-Ev4Norollover(SwitchOneventsINinterval)-1)) ))])
  disp(['Sum of all |diff| events cond with SwitchOnevents: ' num2str(sum( abs(ScanCount(SwitchOneventsINinterval)-(Ev16Norollover(SwitchOneventsINinterval)-Ev4Norollover(SwitchOneventsINinterval)-1)) ))])
  
  [IntervalToSee' Ev4Norollover(IntervalToSee) Ev16Norollover(IntervalToSee) ScanCount(IntervalToSee) Ev16Norollover(IntervalToSee)-Ev4Norollover(IntervalToSee)-1 TrackingOnSwitch(IntervalToSee) ReleasThresholdOK(IntervalToSee) TrackingOn(IntervalToSee)]
  %[Ev4Norollover(FirstOff)+[-10:25]' chanNorollover(Ev4Norollover(FirstOff)+[-10:25]) VcompStart(Ev4Norollover(FirstOff)+[-10:25])  VcompStop(Ev4Norollover(FirstOff)+[-10:25]) Vlogvector(Ev4Norollover(FirstOff)+[-10:25])]
  
  
  disp('See tracker on switch event: ')
  FirstOff=TrakingPointer(1:10);
  IntervalToSee=FirstOff;
  [IntervalToSee Ev4Norollover(IntervalToSee) Ev16Norollover(IntervalToSee) ScanCount(IntervalToSee) Ev16Norollover(IntervalToSee)-Ev4Norollover(IntervalToSee)-1 TrackingOnSwitch(IntervalToSee)]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract photon counts which are on top of a start or stop event
%% For some reson the start events only makes the missmatch worse
AtLeast4events=((bitand(chanNorollover(Ev4Norollover),4)==4).*chanNorollover(Ev4Norollover));
Ev4p1NorolloverPointer=find( AtLeast4events==4+1); 
Ev4p2NorolloverPointer=find( AtLeast4events==4+2); 
Ev4p1Norollover=Ev4Norollover(Ev4p1NorolloverPointer);
Ev4p2Norollover=Ev4Norollover(Ev4p2NorolloverPointer);

AtLeast16events=((bitand(chanNorollover(Ev16Norollover),16)==16).*chanNorollover(Ev16Norollover));
Ev16p1NorolloverPointer=find( AtLeast16events==16+1); 
Ev16p2NorolloverPointer=find( AtLeast16events==16+2); 
Ev16p1Norollover=Ev16Norollover(Ev16p1NorolloverPointer);
Ev16p2Norollover=Ev16Norollover(Ev16p2NorolloverPointer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=size(Ev4Norollover,1)-size(Ev16Norollover,1);
if D~=0
if V
  disp('Missmatch in size between chanels, will trim down to the shortest')
  disp(['Missmatch is: ' num2str(D)])
  disp('________________________________________________')
end
  M=min(size(Ev4Norollover,1),size(Ev16Norollover,1));
  IntervallN=[1:M];%[1:numElements];
  Ev4Norollover=Ev4Norollover(IntervallN);
  Ev16Norollover=Ev16Norollover(IntervallN); 
end
TTRAPDCount=Ev16Norollover-Ev4Norollover-1; 

Countatend=zeros(size(chanNorollover));
Countatend(Ev16Norollover)=TTRAPDCount;
Countatend(Ev16p1Norollover)=Countatend(Ev16p1Norollover)+1;
Countatend(Ev16p2Norollover)=Countatend(Ev16p2Norollover)+1;

TTRAPDCountEnd=Countatend(Ev16Norollover);

CountatStart=zeros(size(chanNorollover));
CountatStart(Ev4Norollover)=TTRAPDCount;
CountatStart(Ev4p1Norollover)=CountatStart(Ev4p1Norollover)+1;
CountatStart(Ev4p2Norollover)=CountatStart(Ev4p2Norollover)+1;

TTRAPDCountStart=CountatStart(Ev4Norollover); %For some reson this do not give a better aligment  

D=size(TTRAPDCount,1)-size(ScanCount,1);
if D~=0
if V
  disp('Missmatch in size between Scan and Fast chanels, will trim down to the shortest')
  disp(['Missmatch is: ' num2str(D)])
  disp('________________________________________________')
end
  M=min(size(TTRAPDCount,1),size(ScanCount,1));
  IntervallN=[1:M];
  TTRAPDCount=TTRAPDCount(IntervallN);
  Ev4Norollover=Ev4Norollover(IntervallN);
  Ev16Norollover=Ev16Norollover(IntervallN);
  
  ScanCount=ScanCount(IntervallN);
  ReleasThresholdOK=ReleasThresholdOK(IntervallN);
  TrackingOnSwitch=TrackingOnSwitch(IntervallN);
  TrackingOn=TrackingOn(IntervallN);
end
Diffcount=abs(TTRAPDCount-ScanCount);

if V
 D=find(Diffcount~=0);  %All diff
 D1=find(and(Diffcount~=0,abs(TrackingOnSwitch)~=1)); % diff not with TrackingOnSwitch 
 D2=find(and(Diffcount~=0,abs(TrackingOnSwitch)==1)); % diff with TrackingOnSwitch
 DTemp=(abs(TrackingOnSwitch)==1);
 D3=find(and(Diffcount~=0,[DTemp(2:end); 0]  ) ); % diff 1 element before TrackingOnSwitch
 D4=find(and(Diffcount~=0,[0; DTemp(1:end-1)]  ) );% diff 1 element after TrackingOnSwitch
 D5=find( Diffcount~=0 & not([0; DTemp(1:end-1)]) & not([DTemp(2:end); 0]) & abs(TrackingOnSwitch)~=1  ); % diff not with TrackingOnSwitch offsets -1 0 or 1
 
 disp('________________________________________________')
 if not(isempty(D))
  disp(['Amount of points differing between Scan and fast channles are: ' num2str(size(D,1)) ' out of ' num2str(size(Diffcount,1))])
  disp(['Largest missmatch is:' num2str(max(Diffcount(D)))])
 end
 
 if not(isempty(D1))
  disp(['Amount of points not asosieted to a change in the TrackingOnSwitch are: ' num2str(size(D1,1))])
  disp(['Largest missmatch is:' num2str(max(Diffcount(D1)))])
 end
 
 if not(isempty(D2))
  disp(['Amount of points asosieted to a change in the TrackingOnSwitch are: ' num2str(size(D2,1))])
  disp(['Largest missmatch is:' num2str(max(Diffcount(D2)))])
   
 end

 if not(isempty(D3))
  disp(['Amount of points asosieted to a change in the TrackingOnSwitch-1 are: ' num2str(size(D3,1))])
  disp(['Largest missmatch is:' num2str(max(Diffcount(D3)))])
 end
 
 if not(isempty(D4))
  disp(['Amount of points asosieted to a change in the TrackingOnSwitch+1 are: ' num2str(size(D4,1))])
  disp(['Largest missmatch is:' num2str(max(Diffcount(D4)))])
 end

 if not(isempty(D5))
  disp(['Amount of points not asosieted to a change in the TrackingOnSwitch [-1 0 1] are: ' num2str(size(D5,1))])
  disp(['Largest missmatch is:' num2str(max(Diffcount(D5)))])
  disp('________________________________________________') 
 end
 
 
 SEE=1;
 IntervalToSee=D(SEE)+[-10:10]';
 [IntervalToSee Ev4Norollover(IntervalToSee) Ev16Norollover(IntervalToSee) ScanCount(IntervalToSee) Ev16Norollover(IntervalToSee)-Ev4Norollover(IntervalToSee)-1 TTRAPDCountStart(IntervalToSee) TTRAPDCountEnd(IntervalToSee) TrackingOnSwitch(IntervalToSee)+10*TrackingOnSwitch(IntervalToSee+1)+100*TrackingOnSwitch(IntervalToSee-1)]
 %[Ev4Norollover(D(SEE))+[-10:20]' chanNorollover(Ev4Norollover(D(SEE))+[-10:20]') Vlogvector(Ev4Norollover(D(SEE))+[-10:20]') Vcomp(Ev4Norollover(D(SEE))+[-10:20]') VcompStart(Ev4Norollover(D(SEE))+[-10:20]') VcompStop(Ev4Norollover(D(SEE))+[-10:20]')] 
 %Interval=find(Vcomp)(3)+[-10:20]';
 %[Interval chanNorollover(Interval) Vlogvector(Interval) VcompT(Interval) VcompStart(Interval) VcompStop(Interval)] 

 
end

if V
%Use to display values and tracking errors
      %IntervallN=D';%+[-10:100];
      %IntervallN=[1:100];
      %[IntervallN' Ev4Norollover(IntervallN) Ev16Norollover(IntervallN) Diffcount(IntervallN) TTRAPDCount(IntervallN) ScanCount(IntervallN) TrackingOn(IntervallN) ReleasThresholdOK(IntervallN) TrackingOnSwitch(IntervallN)]
     %[IntervallN' Diffcount(IntervallN) Ev4Norollover(IntervallN) Ev16Norollover(IntervallN) TTRAPDCount(IntervallN) ScanCount(IntervallN) TrackingOnSwitch(IntervallN)]
end

%%% Go back to the original data 
if LoopIndex==1
  chanNew=chan;
  Ev4New=chanNorolloverPointer(Ev4Norollover);
  Ev16New=chanNorolloverPointer(Ev16Norollover);
  
  chanNew(Ev4)=0;
  chanNew(Ev4New)=chan(Ev4New);
  chanNew(Ev16)=0;
  chanNew(Ev16New)=chan(Ev16New);
  chan=chanNew;
  
  dummy=zeros(size(chan));
  dummy(Ev4New)=1;
  start=dummy;

  dummy=zeros(size(chan));
  dummy(Ev16New)=1;
  stop=dummy;

  Vlogvector=cumsum(start)-cumsum(stop);
  Vlogvector(Ev16)=1;
  Vpointer=find(Vlogvector);%
  
  %StartStop=diff(Vlogvector);
  %find(bitand(chan(find(StartStop==1)+1),4)~=4)
  %find(bitand(chan(find(stop)),16)~=16)
  %[find(bitand(chanNorollover(Ev4Norollover),4)~=4) chanNorollover(Ev4Norollover)(find(bitand(chanNorollover(Ev4Norollover),4)~=4))]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  chan=chan(Vpointer);
  nsync=nsync(Vpointer);

  Ev1=find(bitand(chan,1));%event 1 occured APD 1
  Ev2=find(bitand(chan,2));%event 2 occured APD 2
  Ev4=find(bitand(chan,4));%event 4 occured Start marker
  Ev8=find(bitand(chan,8));%event 8 occured Rollover marker
  Ev16=find(bitand(chan,16));%event 16 occured Stop marker
  %%%%%
end
end

%%Stor old data
nsyncOld=nsync;
chanOld=chan;


%%% Set all start events to 0
% This was wrong, the FPGA sets all starts to 0, but the out put form the FPGA is the value before 0
%  OffsetCorrection_vectorStart=zeros(size(nsync));
%  OffsetCorrection_vectorStart(Ev4)=nsync(Ev4);
%  
%  OffsetCorrection_vectorStop=zeros(size(nsync));
%  OffsetCorrection_vectorStop(Ev16)=nsync(Ev4);
%  
%  OffsetCorrection_vector=zeros(size(nsync));
%  OffsetCorrection_vector=cumsum(OffsetCorrection_vectorStart)-cumsum(OffsetCorrection_vectorStop);
%  OffsetCorrection_vector(Ev16)=nsync(Ev4);
%  
%  nsync=nsync-OffsetCorrection_vector;

TimeBetween=nsync(Ev4);
nsync(Ev4)=0;

%%% WRAPAROUND Correction
%%% Extract markers that contian a channel 8 rollover event (event of the from ?1???)
WRAPAROUND_vector=zeros(size(nsync));
WRAPAROUND_vector(Ev8)=WRAPAROUND;
WRAPAROUND_vector=cumsum(WRAPAROUND_vector);
nsync= nsync + WRAPAROUND_vector;

WRAPAROUND_vectorStart=zeros(size(nsync));
WRAPAROUND_vectorStart(Ev4)=nsync(Ev4);

WRAPAROUND_vectorStop=zeros(size(nsync));
WRAPAROUND_vectorStop(Ev16)=nsync(Ev4);

WRAPAROUND_vector=zeros(size(nsync));
WRAPAROUND_vector=cumsum(WRAPAROUND_vectorStart)-cumsum(WRAPAROUND_vectorStop);
WRAPAROUND_vector(Ev16)=nsync(Ev4);
WRAPAROUND_vector(Ev8)=WRAPAROUND_vector(Ev8)+WRAPAROUND;

nsync= nsync - WRAPAROUND_vector;

%%% WRAPAROUND to Stop time Correction for next start event
StartStopCorrection_vector=zeros(size(nsync));
StartStopCorrection_vector(Ev4(2:end))=nsync(Ev16(1:end-1))-nsync(Ev4(2:end));
StartStopCorrection_vector=cumsum(StartStopCorrection_vector);

truensync=nsync+StartStopCorrection_vector;

%%% Finding the true time for each channel 
nsync=nsync*Resolution;%Each start starts at 0
truensync=truensync*Resolution;%Each start is alignied to the previos stop
TimeBetween=TimeBetween*Resolution;%Time between stop and start modulo 2^16

%%% Correct for Time between measurements
TrackingFlag=or(ReleasThresholdOK,TrackingOn);

if and(BetweenMTimeScan<2^16*Resolution,BetweenMTimeTrack<2^16*Resolution)
BetweenMTimeVector=zeros(size(TrackingFlag));
BetweenMTimeVector=TimeBetween;
BetweenMTimeVector(1)=0; %set the first ellement to 0 to enforce starting at 0

BETWEENTIME_vector=zeros(size(nsync));
BETWEENTIME_vector(Ev4)=BetweenMTimeVector;
BETWEENTIME_vector=cumsum(BETWEENTIME_vector);
else
BetweenMTimeVector=zeros(size(TrackingFlag));
BetweenMTimeVector(find(TrackingFlag==0))=BetweenMTimeScan;
BetweenMTimeVector(find(TrackingFlag==1))=BetweenMTimeTrack;

BETWEENTIME_vector=zeros(size(nsync));
BETWEENTIME_vector(Ev4(2:end))=BetweenMTimeVector(1:end-1);
BETWEENTIME_vector=cumsum(BETWEENTIME_vector);
end 


truensync=truensync+BETWEENTIME_vector;

%[nsync(1:10) chan(1:10) truensync(1:10)*1e3 WRAPAROUND_vector(1:10)*1e3 StartStopCorrection_vector(1:10)*1e3 BETWEENTIME_vector(1:10)*1e3]
%[chan(find(TrackingFlag)(1)-10:find(TrackingFlag)(1)+100) nsyncOld(find(TrackingFlag)(1)-10:find(TrackingFlag)(1)+100)]

%%
if ShowFlag==1
figure(1)
plot(nsync)
hold on
plot(Ev1,nsync(Ev1),'.r');%APD 1
plot(Ev2,nsync(Ev2),'.b');%APD 2
plot(Ev4,nsync(Ev4),'.g');%Start
plot(Ev8,nsync(Ev8),'.k');%Rollover
plot(Ev16,nsync(Ev16),'.m');%Stop
xlabel('Event index')
ylabel('Folded Time [s]')
legend('All events','APD1', 'APD2', 'Start', 'Rollover', 'Stop')

figure(2)
plot(truensync)
hold on
plot(Ev1,truensync(Ev1),'.r')
plot(Ev2,truensync(Ev2),'.b')
plot(Ev4,truensync(Ev4),'.g')
plot(Ev8,truensync(Ev8),'.k')
plot(Ev16,truensync(Ev16),'.m')
xlabel('Event index')
ylabel('Time [s]')
legend('All events','APD1', 'APD2', 'Start', 'Rollover', 'Stop')
end


%% Events related to APD channels 
EventStructure.TimeEvents.Sync=truensync;%[S] All events
EventStructure.APDEvents.Ch1Pos=Ev1;%All positions in truesync related to ch1
EventStructure.APDEvents.Ch2Pos=Ev2;%All positions in truesync related to ch2
EventStructure.StartEvents=Ev4;%All positions in truesync related Start 
EventStructure.StopEvents=Ev16;%All positions in truesync related Stop
EventStructure.Old.nsync=nsyncOld;
EventStructure.Old.chan=chanOld;

