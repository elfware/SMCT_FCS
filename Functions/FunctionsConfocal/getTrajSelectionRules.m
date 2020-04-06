function Trajectories=getTrajSelectionRules(Scanobject,Trajectories,DataSelectionRules)
%DataSelectionRules.Trajectory.Mode
%DicplacmentThreshold=DataSelectionRules.Trajectory.CutDicplacment(1);
%DicplacmentCutAfter=DataSelectionRules.Trajectory.CutDicplacment(2);

APDThreshold=DataSelectionRules.Trajectory.CutTraj(1);
CutAfter=DataSelectionRules.Trajectory.CutTraj(2);

APDThresholdLaserHigh=DataSelectionRules.Trajectory.CutTrajWhenLaserHigh(1);
CutAfterLaserHigh=DataSelectionRules.Trajectory.CutTrajWhenLaserHigh(2);

nn=1;
for n=1:size(Trajectories.Traj,1)
  eval(['TrajData=Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '{n};']);
  APDRatesignal=sum(TrajData(:,6:7),2)./TrajData(:,10);

  %CutingVector=conv(ones(CutAfter,1),APDRatesignal>APDThreshold)(CutAfter:end+1-CutAfter)>0;  
  APDRatesignalLow=APDRatesignal(find(not(TrajData(:,8))));
  APDRatesignalHigh=APDRatesignal(find(TrajData(:,8)));
     
  CutingVectorLow=conv(ones(CutAfter,1),APDRatesignalLow>APDThreshold)(CutAfter:end+1-CutAfter)>0;  
  CutingVectorHigh=conv(ones(CutAfterLaserHigh,1),APDRatesignalHigh>APDThresholdLaserHigh)(CutAfterLaserHigh:end+1-CutAfterLaserHigh)>0;  

  CutingVector=[CutingVectorLow; ones(CutAfter,1); CutingVectorHigh];
  
  if CutingVector(1)==1
    CutHereVector=find(diff(CutingVector)==-1);
    if isempty(CutHereVector)==0
      CutHere=CutHereVector(1);
    else
      CutHere=size(APDRatesignal,1);
    end
    
    %Cut according to max displacment 
    %Displacment=sqrt(sum((diff(TrajData(:,2:3),1,1).^2),2));
    %CutingVectorDicplacment=conv(ones(DicplacmentCutAfter,1),Displacment<DicplacmentThreshold)(DicplacmentCutAfter:end+1-DicplacmentCutAfter)>0;  
    
    
    %%Reselect trajectory based on minimum time length requirement 
    %%And remove elements in Traj
    Traj=[];
    if strcmp(DataSelectionRules.Trajectory.Mode,'RollingAv') 
      Traj(n,2:3)=[Trajectories.Traj(n,2) Trajectories.Traj(n,2)+(DataSelectionRules.Trajectory.wndw-1)+(CutHere-1)];
      Traj(n,1)=Scanobject.ScanData.ScanMatrix(Traj(n,3),1)-Scanobject.ScanData.ScanMatrix(Traj(n,2),1);
    elseif strcmp(DataSelectionRules.Trajectory.Mode,'FixWindowAv') 
      Traj(n,2:3)=[Trajectories.Traj(n,2) Trajectories.Traj(n,2)+(DataSelectionRules.Trajectory.wndw*CutHere)-1];
      Traj(n,1)=Scanobject.ScanData.ScanMatrix(Traj(n,3),1)-Scanobject.ScanData.ScanMatrix(Traj(n,2),1);    
    else 
    disp('Wrong Mode. Plese use RollingAv or FixWindowAv')
    end
    
    NewTrajectoryTimeLength=Traj(n,1);%diff(Trajectories.TrajAvriage{n}([1 CutHere],1));
    if NewTrajectoryTimeLength> DataSelectionRules.Trajectory.MinTrajectoryTimeLength
      TrajDataNew.TrajAvriage{nn}(:,:)=TrajData(1:CutHere,:);
      %TrajectoriesNew.TrajPosEstimation.TrajAvriage{nn}(:,:)=TrajData(1:CutHere,:);
      TrajDataNew.Traj(nn,1:3)=Traj(n,:);
      nn=nn+1;
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