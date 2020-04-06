function Trajectories=getTrajAvraging(Scanobject,Trajectories,DataSelectionRules,nn)

wndw=DataSelectionRules.Trajectory.wndw;
if strcmp(DataSelectionRules.Trajectory.Mode,'FixWindowAv') 
  wndwStep=wndw;
elseif strcmp(DataSelectionRules.Trajectory.Mode,'RollingAv')
  wndwStep=1;
else
disp('Wrong Mode value used in getTrajAvraging, use: FixWindowAv or RollingAv')
end

for n=nn; 
Interval=Trajectories.Traj(n,2):Trajectories.Traj(n,3);

  X=Scanobject.ScanData.ScanMatrix(Interval,2);
  Y=Scanobject.ScanData.ScanMatrix(Interval,3);
  Z=Scanobject.ScanData.ScanMatrix(Interval,1);
  TotCounts=Scanobject.ScanData.ScanMatrix(Interval,4)+Scanobject.ScanData.ScanMatrix(Interval,5);
  LaserFlag=Scanobject.ScanData.ScanMatrix(Interval,9);
  
  %%% Timing 
  TrajTime=Z(end)-Z(1);

  if wndw<size(Z,1)
    WindowTime=Z(wndw)-Z(1);
  else
    disp('Fewer elemets compared to requested time window')
    WindowTime=Z(end)-Z(1);
  end

  %% Calculating avriage waighted by the photon counts
  %% if wndwStep=1 then we have a rolling mean, 
  %% if wndwStep=wndw then we have the mean ofver fix steps 
  XXW=X.*TotCounts;
  YYW=Y.*TotCounts;
  
  SAPD=conv(ones(wndw,1), TotCounts);
  SAPD=SAPD(wndw:wndwStep:end+1-wndw);

  XXRM=conv(ones(wndw,1), XXW);
  XXRM=XXRM(wndw:wndwStep:end+1-wndw)./SAPD;
  YYRM=conv(ones(wndw,1), YYW);
  YYRM=YYRM(wndw:wndwStep:end+1-wndw)./SAPD;
  TTRM=Z(wndw:wndwStep:end);
  %%%%%%

  %%%% Polarization estimation 
  RMCAPD1=conv(Scanobject.ScanData.ScanMatrix(Interval,4),ones(wndw,1));
  RMCAPD1=RMCAPD1(wndw:wndwStep:end+1-wndw);
  
  RMCAPD2=conv(Scanobject.ScanData.ScanMatrix(Interval,5),ones(wndw,1));
  RMCAPD2=RMCAPD2(wndw:wndwStep:end+1-wndw);
  
  WindowTime=conv(Z,[1; zeros(wndw-2,1); -1]);
  WindowTime=WindowTime(wndw:wndwStep:end+1-wndw);
  AcumilatedMeasurementTime=conv(ones(size(Z)),ones(wndw,1));
  AcumilatedMeasurementTime=AcumilatedMeasurementTime(wndw:wndwStep:end+1-wndw)*Scanobject.ScanSpeeds.MeasurementOnTime;
  
  Pol=(RMCAPD1-RMCAPD2)./(RMCAPD1+RMCAPD2);
  
  %%Error assuming poisson statistics 
  DAPD_1=sqrt(RMCAPD1);
  DAPD_1(find(DAPD_1==0))=1;
  DAPD_2=sqrt(RMCAPD2);
  DAPD_2(find(DAPD_2==0))=1;
  
  DPDAPD_1=2*RMCAPD2./((RMCAPD1+RMCAPD2).^2);
  DPDAPD_2=-2*RMCAPD1./((RMCAPD1+RMCAPD2).^2);
  
  DPol=sqrt((DAPD_1.*DPDAPD_1).^2+(DAPD_2.*DPDAPD_2).^2);

  %%Laser LaserFlag
  LaserFlag=LaserFlag(wndw:wndwStep:end);
  
  %Adding all information to the Trajectori object
 
  INFOText={[DataSelectionRules.Trajectory.Mode ' mean with a window of ' num2str(wndw) ' elements']; ...
            ['1.' DataSelectionRules.Trajectory.Mode ' mean time [s]']; ...
            ['2.' DataSelectionRules.Trajectory.Mode ' mean x ']; ...
            ['3.' DataSelectionRules.Trajectory.Mode ' mean y']; ...
            ['4.' DataSelectionRules.Trajectory.Mode ' mean polarization']; ...
            ['5.' DataSelectionRules.Trajectory.Mode ' mean poisson polarization error']; ...
            ['6.APD1 counts over the' DataSelectionRules.Trajectory.Mode '  window']; ...
            ['7.APD2 counts over the' DataSelectionRules.Trajectory.Mode '  window']; ...
            ['8.Laser high Flag']; ...
            ['9.' DataSelectionRules.Trajectory.Mode ' window time [s]'];...
            ['10.' DataSelectionRules.Trajectory.Mode ' acumilated measurement time [s]']};
  
  eval(['Trajectories.INFO.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '=INFOText;']);
  eval(['Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '{n}=[TTRM XXRM YYRM Pol DPol RMCAPD1 RMCAPD2 LaserFlag WindowTime AcumilatedMeasurementTime];']);
  eval(['Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '{n}=Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '{n}(find(prod(not(isnan(Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '{n})),2)),:);']);
end

