function Trajectories=getTrajSelectionRules_V2(Scanobject,Trajectories,DataSelectionRules)
%DataSelectionRules.Trajectory.Mode
DicplacmentThreshold=DataSelectionRules.Trajectory.CutDicplacment(1);
DicplacmentCutAfter=DataSelectionRules.Trajectory.CutDicplacment(2);
DicplacmentThresholdmax=DataSelectionRules.Trajectory.CutDicplacment(3);

APDThreshold=DataSelectionRules.Trajectory.CutTraj(1);
CutAfter=DataSelectionRules.Trajectory.CutTraj(2); %Cuting tolerance 


APDThresholdLaserHigh=DataSelectionRules.Trajectory.CutTrajWhenLaserHigh(1);
CutAfterLaserHigh=DataSelectionRules.Trajectory.CutTrajWhenLaserHigh(2);

nn=1;
for n=1:size(Trajectories.Traj,1)

  eval(['TrajData=Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '{n};']);
  APDRatesignal=sum(TrajData(:,6:7),2)./TrajData(:,10);
  
  CutingVector=zeros(size(TrajData,1),1);
  
  APDRatesignalLow=APDRatesignal(find(not(TrajData(:,8))));
  APDRatesignalHigh=APDRatesignal(find(TrajData(:,8)));
  
  if not(isempty(APDRatesignalLow))
  CutingVectorLow=conv(ones(CutAfter,1),double(APDRatesignalLow>APDThreshold));
  CutingVectorLow=CutingVectorLow(CutAfter:end+1-CutAfter)>0;  
  else
      CutingVectorLow=[];
  end
  
  if not(isempty(APDRatesignalHigh))
  CutingVectorHigh=conv(ones(CutAfterLaserHigh,1),double(APDRatesignalHigh>APDThresholdLaserHigh));
  CutingVectorHigh=CutingVectorHigh(CutAfterLaserHigh:end+1-CutAfterLaserHigh)>0;  
  else
      CutingVectorHigh=[];
  end
  
  CutingVectorAPD=[CutingVectorLow; ones(CutAfter-1,1); CutingVectorHigh ;ones(CutAfterLaserHigh-1,1)];
  
  %Cut according to max Displacment
  DisplacmentSpeed=[sqrt(sum((diff(TrajData(:,2:3),1,1).^2),2))./diff(TrajData(:,1)); 0];
  
  CutingVectorDicplacment=conv(ones(DicplacmentCutAfter,1),double(DisplacmentSpeed<DicplacmentThreshold));
  CutingVectorDicplacment=CutingVectorDicplacment(DicplacmentCutAfter:end+1-DicplacmentCutAfter)>0;
  
  CutingVectorDicplacment=[CutingVectorDicplacment; ones(DicplacmentCutAfter-1,1)];
  CutingVectorDicplacmentMax=DisplacmentSpeed<DicplacmentThresholdmax;
  
  %Accumileted cuting
  %[n size(CutingVectorAPD) size(CutingVectorDicplacment) size(CutingVectorDicplacmentMax)]
  CutingVector=CutingVectorAPD.*CutingVectorDicplacment.*CutingVectorDicplacmentMax;
  %CutingVector=CutingVectorAPD.*CutingVectorDicplacment;
    
  
  if CutingVector(1)==0
    CutingVectorDiff=[0; diff(CutingVector)];
  elseif CutingVector(1)==1
    CutingVectorDiff=[1; diff(CutingVector)];
  end
  
  if CutingVector(end-1)==1
    CutingVectorDiff(end)=-1;
  elseif CutingVector(end-1)==0
    CutingVectorDiff(end)=0;
  end
  
  %Find if some parts are passing the time lenngth requierment
  CutingVectorDiffStart=find(CutingVectorDiff==1);
  CutingVectorDiffStop=find(CutingVectorDiff==-1);
  
  TimeDiff=TrajData(CutingVectorDiffStop,1)-TrajData(CutingVectorDiffStart,1);
  TimeDiffpointerPass=find(TimeDiff>DataSelectionRules.Trajectory.MinTrajectoryTimeLength);
 
  CutingVectorStarPointer=CutingVectorDiffStart(TimeDiffpointerPass);
  CutingVectorStopPointer=CutingVectorDiffStop(TimeDiffpointerPass);
  
  if not(isempty(TimeDiffpointerPass))
  for TrajSectionPased=[1:size(TimeDiffpointerPass,1)]
    
    %%Reselect trajectory based on minimum time length requirement 
    %%And remove elements in Traj
    %Traj(n,:)=Trajectories.Traj(n,:);
    Traj=[];
    if strcmp(DataSelectionRules.Trajectory.Mode,'RollingAv') 
      Traj(n,2:3)=[Trajectories.Traj(n,2)+(0-1)+(CutingVectorStarPointer(TrajSectionPased)-1) Trajectories.Traj(n,2)+(DataSelectionRules.Trajectory.wndw-1)+(CutingVectorStopPointer(TrajSectionPased)-1)];
      Traj(n,1)=Scanobject.ScanData.ScanMatrix(Traj(n,3),1)-Scanobject.ScanData.ScanMatrix(Traj(n,2),1);
    elseif strcmp(DataSelectionRules.Trajectory.Mode,'FixWindowAv') 
      Traj(n,2:3)=[Trajectories.Traj(n,2)+(DataSelectionRules.Trajectory.wndw*CutingVectorStarPointer(TrajSectionPased))-1 Trajectories.Traj(n,2)+(DataSelectionRules.Trajectory.wndw*CutingVectorStopPointer(TrajSectionPased))-1];
      Traj(n,1)=Scanobject.ScanData.ScanMatrix(Traj(n,3),1)-Scanobject.ScanData.ScanMatrix(Traj(n,2),1);    
    else 
    disp('Wrong Mode. Plese use RollingAv or FixWindowAv')
    end
    
    NewTrajectoryTimeLength=Traj(n,1);%diff(Trajectories.TrajAvriage{n}([1 CutHere],1));
    if NewTrajectoryTimeLength> DataSelectionRules.Trajectory.MinTrajectoryTimeLength %Double check that the time requierment is passed
      TrajDataNew.TrajAvriage{nn}(:,:)=TrajData(CutingVectorStarPointer(TrajSectionPased):CutingVectorStopPointer(TrajSectionPased),:);
      %TrajectoriesNew.TrajPosEstimation.TrajAvriage{nn}(:,:)=TrajData(1:CutHere,:);
      TrajDataNew.Traj(nn,1:3)=Traj(n,:);
      nn=nn+1;
    end
  end
  end
end

if exist('TrajDataNew')
  eval(['Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '=TrajDataNew.TrajAvriage;']);
  Trajectories.Traj=TrajDataNew.Traj;
else
  eval(['Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '=[];']);
  Trajectories.Traj=[];
end

%  %Select Trajectories where at least one point is passing the selection criterias 
%  TrajectoriesT=[];
%  for n=[1:size(Trajectories,1)]
%    T= sum((Trajectories(n,1)<=Index) & ( Trajectories(n,2)>=Index) )~=0;
%    CutingVector=conv(ones(CutTrajAfter,1),sum(TotData(Trajectories(n,1):Trajectories(n,2),4:5),2)>CutTraj)(1:size(Trajectories(n,1):Trajectories(n,2),2))>0;  
%    pasCondition=T & sum(CutingVector)>0;
%    if pasCondition
%      CutingStartPositions=Trajectories(n,1)+find(diff(CutingVector)==1)-1;
%      if isempty(CutingStartPositions)
%      CutingStartPositions=Trajectories(n,1);
%      end
%      CutingEndPositions=Trajectories(n,1)+find(diff(CutingVector)==-1)-1;
%      if isempty(CutingEndPositions)
%      CutingEndPositions=Trajectories(n,2);
%      end
%      
%      if CutingStartPositions(1)>CutingEndPositions(1)
%      CutingStartPositions=[Trajectories(n,1); CutingStartPositions];
%      end
%      if CutingStartPositions(end)>CutingEndPositions(end)
%      CutingEndPositions=[CutingEndPositions; Trajectories(n,2)];
%      end
%      
%      if size(CutingStartPositions,1)~=size(CutingEndPositions,1)
%      disp('Trajectory cutting vectors are not of same length')
%      end
%      for m=[1:size(CutingEndPositions,1)]
%      	    TrajectoriesT=[TrajectoriesT; [CutingStartPositions(m) CutingEndPositions(m)]];
%      end
%    end
%  end