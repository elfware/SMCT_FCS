function Trajectories=getTrajectories(Scanobject,DataSelectionRules)
verbos=DataSelectionRules.verbos;

%FindTrajectorys
TrajectoryThreshold=Scanobject.ScanData.ScanMatrix(:,6) & Scanobject.ScanData.ScanMatrix(:,7);
  
StartTrajectory=find((diff(TrajectoryThreshold)==1));
EndTrajectory=find(diff(TrajectoryThreshold)==-1);
  
StartTrajectory=StartTrajectory(1:size(EndTrajectory,1)); %Trimdown to equal length

TrajectoryTimeLength=Scanobject.ScanData.ScanMatrix(EndTrajectory,1)-Scanobject.ScanData.ScanMatrix(StartTrajectory,1);

%% Remove Trajectories shorter than MinTrajectoryTimeLength
SelectTrajectoryLength=find(TrajectoryTimeLength> DataSelectionRules.Trajectory.MinTrajectoryTimeLength);
  
  
StartTrajectoryOfIntrest=StartTrajectory(SelectTrajectoryLength);
EndTrajectoryOfIntrest=EndTrajectory(SelectTrajectoryLength);
TrajectoryTimeLengthOfIntrest=TrajectoryTimeLength(SelectTrajectoryLength);

if verbos
  disp(['Found ' num2str(size(StartTrajectoryOfIntrest,1)) ' trajectories'])
end

Trajectories.INFO.Traj={['Found ' num2str(size(StartTrajectoryOfIntrest,1)) ' trajectories']; '1.Trajectory time length [s]'; '2.Start index'; '3.End index'};
Trajectories.Traj=[TrajectoryTimeLengthOfIntrest StartTrajectoryOfIntrest EndTrajectoryOfIntrest];
  
  
