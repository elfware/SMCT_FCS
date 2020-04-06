function Trajectories=getTrajPosEstimation(Scanobject,Trajectories,DataSelectionRules,nn)


if or(strcmp(DataSelectionRules.Trajectory.Mode,'FixWindowAv'),strcmp(DataSelectionRules.Trajectory.Mode,'RollingAv')) 
  Trajectories=getTrajAvraging(Scanobject,Trajectories,DataSelectionRules,nn);
elseif strcmp(DataSelectionRules.Trajectory.Mode,'DirectEstimation3Points')
  Trajectories=getTrajDirectEstimation3Points(Scanobject,Trajectories,DataSelectionRules,nn);
else
disp('Wrong Mode value use: FixWindowAv, RollingAv or DirectEstimation3Points')
end
