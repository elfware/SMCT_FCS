function [Scanobject DataSelectionRules]=BasicPhotonStat(Scanobject,DataSelectionRules,p)

AllMean=mean(sum(Scanobject.ScanData.ScanMatrix(:,4:5),2));
%AllPhotonRate=AllMean/mean(diff(Scanobject.ScanData.ScanMatrix(:,1)));
AllPhotonRate=mean(sum(Scanobject.ScanData.ScanMatrix(:,4:5),2)/Scanobject.ScanSpeeds.MeasurementOnTime);

OfIntrest=find(sum(Scanobject.ScanData.ScanMatrix(:,4:5),2)<DataSelectionRules.Trajectory.BackgroundThrechold);
%ThrecholdPhotonRate=mean(sum(Scanobject.ScanData.ScanMatrix(OfIntrest,4:5),2))/mean(diff(Scanobject.ScanData.ScanMatrix(OfIntrest,1)));
ThrecholdPhotonRate=mean(sum(Scanobject.ScanData.ScanMatrix(OfIntrest,4:5),2)/Scanobject.ScanSpeeds.MeasurementOnTime);

%Mean photon count over all data 
MenaPolAll=(sum(Scanobject.ScanData.ScanMatrix(:,4),1)-sum(Scanobject.ScanData.ScanMatrix(:,5),1))./sum(sum(Scanobject.ScanData.ScanMatrix(:,4:5),1));

%Mean photon count over data with more then the background threchold
OfIntrest=find(sum(Scanobject.ScanData.ScanMatrix(:,4:5),2)>=DataSelectionRules.Trajectory.BackgroundThrechold);
MenaPolover=(sum(Scanobject.ScanData.ScanMatrix(OfIntrest,4),1)-sum(Scanobject.ScanData.ScanMatrix(OfIntrest,5),1))./sum(sum(Scanobject.ScanData.ScanMatrix(OfIntrest,4:5),1));


if p
  disp(['Mean photon count over all data: ' num2str(AllMean)])
  disp(['Mean photon rate over all data: ' num2str(AllPhotonRate) ' Hz'])
  disp(['Mean photon rate over data with less then ' num2str(DataSelectionRules.Trajectory.BackgroundThrechold) ' counts: ' num2str(ThrecholdPhotonRate) ' Hz'])
  disp(['Mena polarization over all data:' num2str(MenaPolAll)])
  disp(['Mena polarization over data with more then ' num2str(DataSelectionRules.Trajectory.BackgroundThrechold) ' counts: ' num2str(MenaPolover)])
end 

Scanobject.ScanData.BasicStat.AllMean=AllMean;
Scanobject.ScanData.BasicStat.AllPhotonRate=AllPhotonRate;
Scanobject.ScanData.BasicStat.ThrecholdPhotonRate=ThrecholdPhotonRate;
Scanobject.ScanData.BasicStat.AllMeanPol=MenaPolAll;
Scanobject.ScanData.BasicStat.ThrecholdMeanPol=MenaPolover;

if DataSelectionRules.Trajectory.CutTraj(3)==1
DataSelectionRules.Trajectory.CutTraj(1)=DataSelectionRules.Trajectory.CutTraj(1)+ThrecholdPhotonRate;
DataSelectionRules.Trajectory.CutTrajWhenLaserHigh(1)=ThrecholdPhotonRate+DataSelectionRules.Trajectory.CutTrajWhenLaserHigh(4);
end