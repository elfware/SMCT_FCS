%clear all
close all
%V10: Throw away beginning of all trajectories (tracking artifact when
%capturing a dot) specfied by parameter trajStarTime
%% APD_1=H and APD_2=V.
%addpath(genpath('/home/samba/elflab/Projects/polarization/Pol_MatLabScript/Functions'))
%  javaaddpath ("/path/to/xerces-2_11_0/xercesImpl.jar");
%  javaaddpath ("/path/to/xerces-2_11_0/xml-apis.jar");


pathnamesold={'/Calibration_deg0_HRtoL/ScanData/'; ...
    '/Pos1M1_deg0_HRtoL_Lockin8_release0,4/ScanData/'; ...
    '/Pos1M2_deg0_HRtoL_Lockin8_release0,4/ScanData/'; ...
    '/Pos2M1_deg0_HRtoL_Lockin8_release0,4/ScanData/'; ...
    '/Pos2M2_deg0_HRtoL_Lockin8_release0,4/ScanData/'; ...
    '/Pos3M1_deg0_HRtoL_Lockin8_release0,4/ScanData/'}; %Något fel i denna fil

pathnames=cell(numel(pathnamesold),1);
for i=1:numel(pathnamesold)
    pathnames{i}=[currDataPath pathnamesold{i}];
end


for Scanfiles=[2:4];%[1:size(pathnames,1)]
disp('________________________________________________________')
disp(pathnames{Scanfiles})
%Definitions
DirCont=dir(pathnames{Scanfiles});

Scanobject.DataPath.pathname=pathnames{Scanfiles,:};
Scanobject.DataPath.filenameScanData=DirCont(4).name;
Scanobject.DataPath.filenameTTTRData=DirCont(3).name;
Scanobject.DataPath.filenameXMLFile=DirCont(5).name;

%XMLData=xmlread([pathnames(Scanfiles,:) Scanobject.DataPath.filenameXMLFile]);

Scanobject.PosCalibration.Max=[21129 29550];
Scanobject.PosCalibration.Offset=[998 1003];
Scanobject.PosCalibration.ScaleFactor=[1 1];%[1.2287 1.7425];
Scanobject.PosCalibration.PositionLength=20;%[um]

%APDShift  %1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
ShiftValue=[2 1 1 1 2 1];
Scanobject.PosCalibration.APDShift=ShiftValue(Scanfiles);%[um]

Scanobject.ScanSpeeds.TTTRResolution=5e-9;%[s] 200Mhz clock
%            1
ScanSpeed=[1900];
Scanobject.ScanSpeeds.ScanSpeed=ScanSpeed(1);%sampels/second 
Scanobject.ScanSpeeds.ScanSpeedTracking=1/(33*20e-6);%sampels/second

Scanobject.ScanSpeeds.LaserOnTime=500e-6;%[s]
Scanobject.ScanSpeeds.MeasurementOnTime=500e-6;%[s]


%%Parameters for Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Data Selection 

DataSelectionRules.APDThreshold=[0 0]; %First threshold for each APD, 

DataSelectionRules.PhotonCounts.CountRange=[0 inf];% Second threshold for the sum APD_1+APD_2
DataSelectionRules.PhotonCounts.ErrorBound=[0 inf];%Select data within the error bound  (LowerErrorBound < UpperErrorBound)
DataSelectionRules.PositionTime.xBounds=[-1 21]; %[um] 
DataSelectionRules.PositionTime.yBounds=[-1 21]; %[um] 
DataSelectionRules.PositionTime.TBounds=[0 inf];  %[s]

%% Trajectory Building 
DataSelectionRules.Trajectory.PosStepDiviation=inf;%[um] 
DataSelectionRules.Trajectory.CutDicplacment=[5 10 100];%max um/s (first) over max amount of steps (second)- And last is absolut max um/s for any step
DataSelectionRules.Trajectory.MinTrajectoryTimeLength=0.2;%[S]
DataSelectionRules.Trajectory.BackgroundThrechold=12;%[counts] is the minimum photon count what we like to interpret as a signal, this is used for calculating a real background
DataSelectionRules.Trajectory.CutTraj=[100 1 1]; %Cut when counts rate [1/s] in window is less then [first value ] over an interval of [second value] amount of frames defined by wndw 
					           %if flag 1 is set the count rate will be set by the function 'BasicPhotonStat' which will calculate a bakgrond and add [first value ]
DataSelectionRules.Trajectory.CutTrajWhenLaserHigh=[1000 1 1 300];%Cut when counts rate [1/s] in window is less then [first value ] over an interval of [second value] amount of frames defined by wndw 
					           %if flag 1 is set the count rate will be set by the function 'BasicPhotonStat'
					           %as the value for low laser + [forth value]
DataSelectionRules.Trajectory.Mode='FixWindowAv';%Use: FixWindowAv, RollingAv or Directestimation3Points DirectEstimation3Points
DataSelectionRules.Trajectory.wndw=12;%Amount of points to use in avreage for FixWindowAv and RollingAv 
DataSelectionRules.Trajectory.Sigma=1.1*(0.21*540e-9/1.45)*1e6;%[um] Used if Directestimation3Points is used.
DataSelectionRules.verbos=0;
DataSelectionRules.Trajectory.TrimStartAndEndOfTraj=[0 0];%[s] Trim off time from beiging and end of found trajectories
% paramters for classifciation of sliding
classparams.varwindow=5;%10each
classparams.PrincRelVar=[0.8 1];%0.8 ~ovalnes 
classparams.PrincAbsVar=[0 400];%2   ~sliding size
classparams.AngleInterval=[0-pi/9 0+pi/9];%[min max] and use the interval -pi/2 to pi/2
classparams.minNotStuck=0.01;%0.5
classparams.minMinWindowVar=0;%0.0005;
classparams.maxMeanDiff=inf;%0.3; %Maximu average step distance in the trajecotry (in abosulte number, + or - doesn't matter) this is set to get trajectories that are sliding and not just moving with direction of flow
%%%%%%%%%%

%%Plot scan data
DataSelectionRules.Plot.PlotMode={'3dScatter' 'jet'};% use 3dScatter, 2dScatter, 2dPix or 2dPixTrajectoryMovie and colormap: gray, jet, ...
DataSelectionRules.Plot.Mode='XYTimePol';% use [Int Pol Angle TimeInt TimePol] for intensity [counts], polarization [p], angle [deg], scan Time [ms]
%p=1;%plot image and trajectories
%Print=0;
%SeeTrajectory=1; %use if PlotMode= 2dPixTrajectoryMovie, it specifies which trajectory to image

%MAVPlot=[1 0];%adds rolling avreage data to the plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Scanobject=getScanData(Scanobject);
[Scanobject DataSelectionRules]=BasicPhotonStat(Scanobject,DataSelectionRules,1);

Trajectories=getTrajectories(Scanobject,DataSelectionRules);
Scanobject=getScanDataSelectionRules(Scanobject,DataSelectionRules);

%% Evaluate other traj estimators
TrajectoryToUse=[1:size(Trajectories.Traj,1)];
DataSelectionRules.Trajectory.Mode='FixWindowAv';
Trajectories=getTrajPosEstimation(Scanobject,Trajectories,DataSelectionRules,TrajectoryToUse);

%  DataSelectionRules.Trajectory.Mode='RollingAv';
%  Trajectories=getTrajPosEstimation(Scanobject,Trajectories,DataSelectionRules,TrajectoryToUse);

%  DataSelectionRules.Trajectory.Mode='DirectEstimation3Points';
%  Trajectories=getTrajPosEstimation(Scanobject,Trajectories,DataSelectionRules,TrajectoryToUse);
  
%%%%

getTrajMeanPol(Scanobject,Trajectories,0,0.3);
disp(['Number of trajectories found without selection rules: ' num2str(size(Trajectories.Traj,1))]);

%DataSelectionRules.Plot.Mode='XYTimePol';
%ShowTraj(Scanobject,Trajectories,DataSelectionRules, 2,'noShowRawData');

%% Select Traj according to selection rules 	
DataSelectionRules.Trajectory.Mode='FixWindowAv';
Trajectories=getTrajSelectionRules_V3(Scanobject,Trajectories,DataSelectionRules);
disp(['Number of trajectories found with selection rules: ' num2str(size(Trajectories.Traj,1))]);

%DataSelectionRules.APDThreshold=[0 0]; %First threshold for each APD, 
%DataSelectionRules.Plot.Mode='XYTimeInt';
%ShowRawData(Scanobject,DataSelectionRules,0)
%DataSelectionRules.Plot.Mode='XYTimePol';
%ShowTraj(Scanobject,Trajectories,DataSelectionRules, 80,'noShowRawData');


if not(isempty(Trajectories.Traj))
  Trajectories=getTrajMeanPol(Scanobject,Trajectories,0,0.1);
  
  Scanobjects{Scanfiles}.Scanobject=Scanobject;
  Scanobjects{Scanfiles}.DataSelectionRules=DataSelectionRules;
  Scanobjects{Scanfiles}.classparams=classparams;
  Scanobjects{Scanfiles}.Trajectories=Trajectories;

  %%Create one object that includes all trajectories 
  Trajectories.Traj(:,end+1)=Scanfiles*ones(size(Trajectories.Traj,1),1);
  if exist('Trajobject')
    Trajobject.Traj=[Trajobject.Traj; Trajectories.Traj];
    %Trajobject.TrajMeanPol=[Trajobject.TrajMeanPol; Trajectories.TrajMeanPol];
    
    for fieldIndex=1:size(fieldnames(Trajectories.TrajPosEstimation),1)
     FieldNames= fieldnames(Trajectories.TrajPosEstimation);
     eval(['Trajobject.TrajPosEstimation.' FieldNames{fieldIndex} '=[Trajobject.TrajPosEstimation.' FieldNames{fieldIndex} ' Trajectories.TrajPosEstimation.' FieldNames{fieldIndex} '];']);
    end
  else
    Trajobject.Traj=Trajectories.Traj;
    %Trajobject.TrajMeanPol=Trajectories.TrajMeanPol;
    Trajobject.INFO=Trajectories.INFO;
    Trajobject.INFO.Traj=[Trajectories.INFO.Traj;...
			  num2str(size(Trajectories.INFO.Traj,1)) '.Scan file index'];
  
    for fieldIndex=1:size(fieldnames(Trajectories.TrajPosEstimation),1)
       FieldNames= fieldnames(Trajectories.TrajPosEstimation);
     eval(['Trajobject.TrajPosEstimation.' FieldNames{fieldIndex} '=Trajectories.TrajPosEstimation.' FieldNames{fieldIndex} ';']);
    end
  end
end

%%clean up
  clear Scanobject Trajectories %DataSelectionRules classparams
  disp(' ')
end
%disp('________________________________________________________')
disp(['Total number of trajectories: ' num2str(size(Trajobject.Traj,1))]);

%%Apply frequency filter cut of for position data: Hopfully removes only high frequency instrument noise 
%DataSelectionRules.Trajectory.Mode='DirectEstimation3Points';%Use: FixWindowAv, RollingAv or DirectEstimation3Points
%Trajobject=PosFrequencyFilter(Trajobject,DataSelectionRules,50);%cuts att 50Mhz=> 20nhz period 

%Trajobject=
%getTrajPosEstimationPolyfit(Trajobject,DataSelectionRules,TrajectoryToUse);


%  %Cal diffusion coeff
DataSelectionRules.Trajectory.Mode='FixWindowAv';%Use: FixWindowAv, RollingAv or DirectEstimation3Points
Trajobject=getDiffCVE(Trajobject,DataSelectionRules,1);%cal DiffCVE, if last arg is 0 no plot, 1 plot x/y diff coff, 2 plot abs x/y diff coff

figure()
a=Trajobject.Traj(:,[4 10]);
plot(a(:,1),a(:,2),'r*')
hold on 
a=Trajobject.Traj(:,[4 9]);
plot(a(:,1),a(:,2),'g*')
grid on 

%% Classification
DataSelectionRules.Trajectory.Mode='FixWindowAv';%Use: FixWindowAv, RollingAv or DirectEstimation3Points
[slidingVector,notstuck,angles,rellat,trajPrinclatent] = slidingClassifyPCAver_V2(Trajobject,DataSelectionRules,classparams,1);
disp(['Number of trajectories passing classification: ' num2str(sum(slidingVector))]);
slidingInd=find(slidingVector);

%%% Load TTTR Data
V=1;
Scanobjects=GetTTTRData(Scanobjects,Trajobject,V);

%%
%% Process TTTR Data
TimeDiffParam=[1e-3 5000];% 0 to first element, second element is the amount of elements 
AutoCorrelationParam=[1 22];% B ncas
Ind=[1:  size(Trajobject.Traj,1)];%slidingInd;
Trajobject=ProcessTTTRData_v9(Scanobjects,Trajobject,Ind, AutoCorrelationParam,1);

%save('ScanTrajData_deg_0_B1ncas22_Main_v10_BackgroundThrechold_12.mat','-v7.3','Scanobjects', 'Trajobject','TimeDiffParam','AutoCorrelationParam','Ind','slidingVector','notstuck','angles','rellat','trajPrinclatent','slidingInd')

% %%
% %AutoCorrelationParam=[20 20];% B ncas
% %j=1:AutoCorrelationParam(1)*AutoCorrelationParam(2);
% %R=2.^floor((j-1)/AutoCorrelationParam(1));
% %T=cumsum(R)*5e-9;
% %plot(T,'-*')

% %%
% close all
% for J=1:1
% x=Trajobject.TTTTRTrajs.Autocorrelation.Timechift(slidingInd(J),1:end);
% ySum=Trajobject.TTTTRTrajs.Autocorrelation.RChSum(slidingInd(J),1:end);
% yCross=Trajobject.TTTTRTrajs.Autocorrelation.RChCross(slidingInd(J),1:end);
% semilogx(x,ySum,'-*b')
% hold on
% semilogx(x,yCross,'-*r')
% grid on
% end
%%
% %%
% close all
% x=Trajobject.TTTTRTrajs.Autocorrelation.Timechift(slidingInd(1),1:end);
% ySum=mean(Trajobject.TTTTRTrajs.Autocorrelation.RChSum(slidingInd,1:end));
% yCross=mean(Trajobject.TTTTRTrajs.Autocorrelation.RChCross(slidingInd,1:end));
% 
% semilogx(x,ySum,'-*b')
% hold on
% semilogx(x,yCross,'-*r')
% grid on
% %% 
 
% %%
% % for J=1:size(Ind,2);
% % x=Trajobject_2.TTTTRTrajs.Autocorrelation.Timechift(Ind(J),2:end);
% % y=Trajobject_2.TTTTRTrajs.Autocorrelation.RChSum(Ind(J),2:end);
% % 
% % plot(x,y)
% % hold on 
% % end
% % axis([1e-6 1e-4 0 5])
% 
% %  figure()
% %  [A B]=hist(Trajobject.Traj(:,4),100);
% %  hist(Trajobject.Traj(:,4),100)
% %  hold on
% %  hist(Trajobject.Traj(slidingInd,4),B,'r')
% %  title(['Mean polarisation:' num2str(mean(Trajobject.Traj(:,4))) ' mean pol over sliding:' num2str(num2str(mean(Trajobject.Traj(slidingInd,4))))])
% 
% %hist(Trajobject.Traj(slidingInd([16 17 19 22 24]),4),B,'r')
% %  [A B]=sort(Trajobject.Traj(slidingInd,4),'descend');  
% 
%% Plot trajecotries
   close all
   SEE=2;
   ViewTraj=slidingInd(SEE);%slidingInd(SEE);
   DataSelectionRules.Plot.Mode='XYTimePol';%% XYTimePol, XPYTimePol, XYPTimePol, XPYPTimePol, XYPolPol, XYTimeInt
   DataSelectionRules.Trajectory.Mode='FixWindowAv';%Use: FixWindowAv, RollingAv or DirectEstimation3Points
   ShowDataPlots(Scanobjects,Trajobject,DataSelectionRules,ViewTraj,'sShowRawData')
   Trajobject.Traj(ViewTraj,4)

 
 
  
