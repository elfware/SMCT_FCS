function A=ShowTraj(Scanobject,Trajectories,DataSelectionRules, nn,PlotArg)


h1=figure;
h2=figure;

for n=nn 
eval(['TrajData=Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode ';']);

%%Restrict datat to ROI
DataOfInterest=TrajData{n}(:,2)>=DataSelectionRules.PositionTime.xBounds(1) & TrajData{n}(:,2)<=DataSelectionRules.PositionTime.xBounds(2);
DataOfInterest=DataOfInterest & TrajData{n}(:,3)>=DataSelectionRules.PositionTime.yBounds(1) & TrajData{n}(:,3)<=DataSelectionRules.PositionTime.yBounds(2);
DataOfInterest=DataOfInterest & TrajData{n}(:,1)>=DataSelectionRules.PositionTime.TBounds(1) & TrajData{n}(:,1)<=DataSelectionRules.PositionTime.TBounds(2);


X=TrajData{n}(DataOfInterest,2); % mean x
Y=TrajData{n}(DataOfInterest,3); % mean y
Z=TrajData{n}(DataOfInterest,1); % mean time
C=TrajData{n}(DataOfInterest,4); % mean polarization
S=TrajData{n}(DataOfInterest,5); % mean polarization poissonn error
RMAPD1=TrajData{n}(DataOfInterest,6); % window counts 
RMAPD2=TrajData{n}(DataOfInterest,7); % window counts 
LaserFlag=TrajData{n}(DataOfInterest,8); %Laser Flag
WindowTime=TrajData{n}(DataOfInterest,10); %measurement time over rolling window

if strcmp(PlotArg,'ShowRawData')
  %Use Data selection rules
  DataOfInterestRaw=Scanobject.ScanData.ScanMatrix(:,4)>=DataSelectionRules.APDThreshold(1) & Scanobject.ScanData.ScanMatrix(:,5)>=DataSelectionRules.APDThreshold(2);
  DataOfInterestRaw=DataOfInterestRaw & sum(Scanobject.ScanData.ScanMatrix(:,4:5),2)>=DataSelectionRules.PhotonCounts.CountRange(1) & sum(Scanobject.ScanData.ScanMatrix(:,4:5),2)<=DataSelectionRules.PhotonCounts.CountRange(2);
  DataOfInterestRaw=DataOfInterestRaw & Scanobject.ScanData.ScanMatrix(:,2)>=DataSelectionRules.PositionTime.xBounds(1) & Scanobject.ScanData.ScanMatrix(:,2)<=DataSelectionRules.PositionTime.xBounds(2);
  DataOfInterestRaw=DataOfInterestRaw & Scanobject.ScanData.ScanMatrix(:,3)>=DataSelectionRules.PositionTime.yBounds(1) & Scanobject.ScanData.ScanMatrix(:,3)<=DataSelectionRules.PositionTime.yBounds(2);
  DataOfInterestRaw=DataOfInterestRaw & Scanobject.ScanData.ScanMatrix(:,1)>=DataSelectionRules.PositionTime.TBounds(1) & Scanobject.ScanData.ScanMatrix(:,1)<=DataSelectionRules.PositionTime.TBounds(2);
  DataOfInterestRaw=DataOfInterestRaw & Scanobject.ScanData.ScanMatrix(:,1)>=Z(1) & Scanobject.ScanData.ScanMatrix(:,1)<=Z(end);
  
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
      ZRaw=Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,1);
      XRaw=Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,2);
      YRaw=Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,3);
      SRaw=1e-10 + 10*sum(Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,4:5),2)/max(sum(Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,4:5),2));
      CRaw=sum(Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,4:5),2);%TotDataT(:,6);%Z;
    case{'XYTimePol'}
      ZRaw=Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,1);
      XRaw=Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,2);
      YRaw=Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,3);
      SRaw=1e-10 + 10*sum(Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,4:5),2)/max(sum(Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,4:5),2));
      CRaw=sum(Scanobject.ScanData.ScanMatrix(DataOfInterestRaw,4:5),2);%TotDataT(:,6);%Z;

      %Z=TotDataT(:,1);%Time
      %S=12*(sum(TotDataT(:,4:5),2)/max(sum(TotDataT(:,4:5),2))).^1.5;%Tot APD counts
      %C=TotDataT(:,6);% Pol
  end
end

%% Find if there is a laser high flag
HighFlagIndex=find(diff(LaserFlag));

if size(HighFlagIndex,1)==0
  dispHighflag=0;
elseif size(HighFlagIndex,1)>1
  disp('More then one Laser High flag!!')
  disp('Only the first Laser High flag is shown')
  HighFlagIndex=HighFlagIndex(1);
  dispHighflag=1;
else
  dispHighflag=1;
end
  


%%%%%%%

  figure(h1)
  subplot(3,1,1)
  plot(Z,X-mean(X),'b')
  hold on
  plot(Z,Y-mean(Y),'r')
  if strcmp(PlotArg,'ShowRawData')
    scatter(ZRaw,XRaw-mean(X),SRaw,CRaw,'filled')
    scatter(ZRaw,YRaw-mean(Y),SRaw,CRaw,'filled')
  end
  if dispHighflag
    plot(Z(HighFlagIndex)*[1 1], [min([X-mean(X); Y-mean(Y)]) max([X-mean(X); Y-mean(Y)])],'-m')
  end
  hold off
  xlabel('Time [s]')
  legend('x-mean(x)','y-mean(y)','location','northeastoutside')
  ylabel('position [um]')
  title(['File:' Scanobject.DataPath.pathname([9:end-1]) ' Trajectory:' num2str(n)])
  grid on 
  axis([Z(1) Z(end) min([X-mean(X); Y-mean(Y)]) max([X-mean(X); Y-mean(Y)])])

  subplot(3,1,2)
  plot(Z,C,'b')	
  hold on 
  plot(Z,C+S,'g')
  plot(Z,C-S,'g')	
  if dispHighflag
    plot(Z(HighFlagIndex)*[1 1], [-1 1],'-m')
  end
  hold off
  ylabel('polarization')
  xlabel('Time [s]')
  legend('Pol','Error','location','northeastoutside')
  grid on 
  axis([Z(1) Z(end) -1 1])

  subplot(3,1,3)
  plot(Z,RMAPD1./WindowTime,'b')
  hold on
  plot(Z,RMAPD2./WindowTime,'r')
  plot(Z,(RMAPD1+RMAPD2)./WindowTime,'m')
  plot(Z,DataSelectionRules.Trajectory.CutTraj(1)*Z./Z,'g')
  if dispHighflag
    plot(Z(HighFlagIndex)*[1 1], [min([RMAPD1./WindowTime; RMAPD2./WindowTime]) max([RMAPD1./WindowTime; RMAPD2./WindowTime])],'-m')
  end
  hold off
  xlabel('Time [s]')
  legend('APD1','APD2','Tot','Background','location','northeastoutside')
  ylabel('Count rate [1/s]')
  grid on 
  axis([Z(1) Z(end) min(([RMAPD1./WindowTime; RMAPD2./WindowTime])) max((RMAPD1+RMAPD2)./WindowTime)])



%    figure(h1)
%    subplot(4,1,1)
%    plot(Z,X-mean(X),'b')
%    hold on
%    plot(Z,Y-mean(Y),'r')
%    if strcmp(PlotArg,'ShowRawData')
%      scatter(ZRaw,XRaw-mean(X),SRaw,CRaw,'filled')
%      scatter(ZRaw,YRaw-mean(Y),SRaw,CRaw,'filled')
%    end
%    if dispHighflag
%      plot(Z(HighFlagIndex)*[1 1], [min([X-mean(X); Y-mean(Y)]) max([X-mean(X); Y-mean(Y)])],'-m')
%    end
%    hold off
%    xlabel('Time [s]')
%    legend('x-mean(x)','y-mean(y)','location','northeastoutside')
%    ylabel('position [um]')
%    title(['File:' Scanobject.DataPath.pathname([9:end-1]) ' Trajectory:' num2str(n)])
%    grid on 
%    axis([Z(1) Z(end) min([X-mean(X); Y-mean(Y)]) max([X-mean(X); Y-mean(Y)])])
%  
%    subplot(4,1,2)
%    plot(Z,C,'b')	
%    hold on 
%    plot(Z,C+S,'g')
%    plot(Z,C-S,'g')	
%    if dispHighflag
%      plot(Z(HighFlagIndex)*[1 1], [-1 1],'-m')
%    end
%    hold off
%    ylabel('polarization')
%    xlabel('Time [s]')
%    legend('Pol','Error','location','northeastoutside')
%    grid on 
%    axis([Z(1) Z(end) -1 1])
%  
%    subplot(4,1,3)
%    plot(Z,RMAPD1./WindowTime,'b')
%    hold on
%    plot(Z,RMAPD2./WindowTime,'r')
%    plot(Z,(RMAPD1+RMAPD2)./WindowTime,'m')
%    plot(Z,DataSelectionRules.Trajectory.CutTraj(1)*Z./Z,'g')
%    if dispHighflag
%      plot(Z(HighFlagIndex)*[1 1], [min([RMAPD1./WindowTime; RMAPD2./WindowTime]) max([RMAPD1./WindowTime; RMAPD2./WindowTime])],'-m')
%    end
%    hold off
%    xlabel('Time [s]')
%    legend('APD1','APD2','location','northeastoutside')
%    ylabel('Count rate [1/s]')
%    grid on 
%    axis([Z(1) Z(end) min(([RMAPD1./WindowTime; RMAPD2./WindowTime])) max((RMAPD1+RMAPD2)./WindowTime)])
%  
%    subplot(4,1,4)
%    plot(Z,RMAPD1,'b')
%    hold on
%    plot(Z,RMAPD2,'r')
%    plot(Z,RMAPD1+RMAPD2,'m')
%    if dispHighflag
%      plot(Z(HighFlagIndex)*[1 1], [min([RMAPD1; RMAPD2]) max([RMAPD1; RMAPD2])],'-m')
%    end
%    hold off
%    xlabel('Time [s]')
%    legend('APD1','APD2','location','northeastoutside')
%    ylabel('Counts')
%    grid on 
%    axis([Z(1) Z(end) min([RMAPD1; RMAPD2]) max(RMAPD1+RMAPD2)])
  
  
  %% Re set dependent on the plot mode
  switch DataSelectionRules.Plot.Mode
    case {'Pol'}
    %  Z=TotDataT(:,6);
    %  S=TotDataT(:,7)*50;
    %  C=Z;
    case {'XYTimeInt'}
      Z=Z;
      X=X;
      Y=Y;
      C=RMAPD1+RMAPD2;
    case {'Angle'}
    %  Z=(180/pi)*(acos(TotDataT(:,6)/sqrt(2))-(pi/4));
    %  S=(180/pi)*abs((sqrt(2)./(2-TotDataT(:,6).^2)).*TotDataT(:,7))*2;
    %  C=Z;
    case{'XYTimePol'}
      Z=Z;
      X=X;
      Y=Y;
      C=C;
    case{'XYPolPol'}
      Z=C;
      X=X;
      Y=Y;
      C=C;
    case{'XPYPTimePol'}
      Z=Z;
      X=mean(X).*ones(size(X));
      Y=mean(Y).*ones(size(Y));
      C=C;
    case{'XYPTimePol'}
      Z=Z;
      X=X;
      Y=mean(Y).*ones(size(Y));
      C=C;
    case{'XPYTimePol'}
      Z=Z;
      X=mean(X).*ones(size(X));
      Y=Y;
      C=C;      
end
  %%%Create 3D tube 
  figure(h2)
  N=6;
  R=Trajectories.Traj(n,1)*0.002/0.8;
  Rz=3*R;
  
  %%%Find plan orthogonal to path gradient
  Normal=diff([X Y Z]);
  Normal=[Normal; Normal(end,:)]; %Add a last elemet to make it the same length as the data set
  Normal= Normal./repmat(sqrt(sum(Normal.^2,2)),1,3);
  Vector1=[-X -Y ((Normal(:,1).*X+Normal(:,2).*Y+Normal(:,3).*Z)./Normal(:,3))-Z];
  Vector1=Vector1./repmat(sqrt(sum(Vector1.^2,2)),1,3);
  Vector2=cross(Normal,Vector1,2);
  Vector1=Vector2./repmat(sqrt(sum(Vector2.^2,2)),1,3);
  
  %%%Generate points on circles on the plan orthogonal to path gradient
  Ang=linspace(0,2*pi,N+1);
  xcir=R*cos(Ang); ycir=R*sin(Ang);
  
  X0=repmat(X,1,N+1); Y0=repmat(Y,1,N+1); Z0=repmat(Z,1,N+1);
  XX=X0+ repmat(Vector1(:,1),1,N+1).*repmat(xcir,size(X0,1),1) + repmat(Vector2(:,1),1,N+1).*repmat(ycir,size(Y0,1),1); 
  YY=Y0+ repmat(Vector1(:,2),1,N+1).*repmat(xcir,size(X0,1),1) + repmat(Vector2(:,2),1,N+1).*repmat(ycir,size(Y0,1),1);
  ZZ=Z0+ Rz*(repmat(Vector1(:,3),1,N+1).*repmat(xcir,size(X0,1),1) + repmat(Vector2(:,3),1,N+1).*repmat(ycir,size(Y0,1),1));
  
  CC=repmat(C,1,N+1);
  surface(XX,YY,ZZ,CC)
  shading interp;
   
  text(mean(X)+0.5,mean(Y),mean(Z),num2str(n))
  if strcmp(PlotArg,'ShowRawData')
    hold on
    scatter3(XRaw,YRaw,ZRaw,SRaw,CRaw,'filled')
  end
  %plot3(X,Y, Z,'-r','LineWidth',0.5);
  grid on
  xlabel('x [um]')
  ylabel('y [um]')
  zlabel('Time [s]')
  title(['File:' Scanobject.DataPath.pathname([9:end-1]) ' Trajectory:' num2str(n) ' time length: ' num2str(Trajectories.Traj(n,1)) ' [s]'])
  colormap('jet')
  colorbar
  view(3)
  if strcmp(PlotArg,'ShowRawData')
    axis([min(XRaw) max(XRaw) min(YRaw) max(YRaw) Z(1) Z(end)])  
  end 
end
hold off
