function A=ShowRawData(Scanobject,DataSelectionRules,p)


%Use Data selection rules
DataOfInterest=Scanobject.ScanData.ScanMatrix(:,4)>=DataSelectionRules.APDThreshold(1) & Scanobject.ScanData.ScanMatrix(:,5)>=DataSelectionRules.APDThreshold(2);
DataOfInterest=DataOfInterest & sum(Scanobject.ScanData.ScanMatrix(:,4:5),2)>=DataSelectionRules.PhotonCounts.CountRange(1) & sum(Scanobject.ScanData.ScanMatrix(:,4:5),2)<=DataSelectionRules.PhotonCounts.CountRange(2);
DataOfInterest=DataOfInterest & Scanobject.ScanData.ScanMatrix(:,2)>=DataSelectionRules.PositionTime.xBounds(1) & Scanobject.ScanData.ScanMatrix(:,2)<=DataSelectionRules.PositionTime.xBounds(2);
DataOfInterest=DataOfInterest & Scanobject.ScanData.ScanMatrix(:,3)>=DataSelectionRules.PositionTime.yBounds(1) & Scanobject.ScanData.ScanMatrix(:,3)<=DataSelectionRules.PositionTime.yBounds(2);
DataOfInterest=DataOfInterest & Scanobject.ScanData.ScanMatrix(:,1)>=DataSelectionRules.PositionTime.TBounds(1) & Scanobject.ScanData.ScanMatrix(:,1)<=DataSelectionRules.PositionTime.TBounds(2);


switch DataSelectionRules.Plot.Mode
  case {'Pol'}
%  Z=TotDataT(:,6);
%  S=TotDataT(:,7)*50;
%  C=Z;
  case {'Int'}
%  Z=TotDataT(:,4)+TotDataT(:,5);
%  S=3*ones(size(TotDataT(:,6)));
%  C=Z;
  case {'Angle'}
%  Z=(180/pi)*(acos(TotDataT(:,6)/sqrt(2))-(pi/4));
%  S=(180/pi)*abs((sqrt(2)./(2-TotDataT(:,6).^2)).*TotDataT(:,7))*2;
%  C=Z;
  case{'XYTimeInt'}
  Z=Scanobject.ScanData.ScanMatrix(DataOfInterest,1);
  X=Scanobject.ScanData.ScanMatrix(DataOfInterest,2);
  Y=Scanobject.ScanData.ScanMatrix(DataOfInterest,3);
  S=1e-10 + 10*sum(Scanobject.ScanData.ScanMatrix(DataOfInterest,4:5),2)/max(sum(Scanobject.ScanData.ScanMatrix(DataOfInterest,4:5),2));
  C=sum(Scanobject.ScanData.ScanMatrix(DataOfInterest,4:5),2);%TotDataT(:,6);%Z;
  case{'TimePol'}
  %Z=TotDataT(:,1);%Time
  %S=12*(sum(TotDataT(:,4:5),2)/max(sum(TotDataT(:,4:5),2))).^1.5;%Tot APD counts
  %C=TotDataT(:,6);% Pol
end

  
switch DataSelectionRules.Plot.PlotMode{1}
  case {'3dScatter'}
  h1=figure;
  scatter3(X,Y,Z,S,C,'filled')
  hold on
  if p
  plot3(X,Y,Z,'-*')
  end
  grid on
  xlabel('x')
  ylabel('y')
  zlabel(DataSelectionRules.Plot.Mode)

 % title(['File:' pathname([9:end-1])])
  colormap(DataSelectionRules.Plot.PlotMode{2})
  colorbar
  axis([DataSelectionRules.PositionTime.xBounds(1) DataSelectionRules.PositionTime.xBounds(2) DataSelectionRules.PositionTime.yBounds(1) DataSelectionRules.PositionTime.yBounds(2) min(Z) max(Z)])
  view(2)
end
