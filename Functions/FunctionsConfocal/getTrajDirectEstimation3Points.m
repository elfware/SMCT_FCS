function Trajectories=getTrajDirectEstimation3Points(Scanobject,Trajectories,DataSelectionRules,nn)
%Select traj
for n=nn
  %RawData 1 2 3 4 5 
  %        t x y A1A2
  RawData=Scanobject.ScanData.ScanMatrix(Trajectories.Traj(n,2):Trajectories.Traj(n,3),:);
  sigma=DataSelectionRules.Trajectory.Sigma;%[um]
  %% Position estimations
  %% Direct estimation of position
  %Para1=getTrajDirectMeasCalc(RawData,sigma,0);%[T xp yp]
  %Para2=getTrajDirectMeasCalc(RawData(IntervaltoUse,:),sigma,1);%[T sigma xp yp]

  RMCAPD1f=RawData(:,4);
  RMCAPD2f=RawData(:,5);
  Cf=sum(RawData(:,4:5),2);
  xf=RawData(:,2);
  yf=RawData(:,3);
  Tf=RawData(:,1);
  LaserFlagf=RawData(:,9);	
  Para=[];
  
  index=1;
  for offset=[1:6:size(xf,1)-6];
    x=xf(offset+[0 2 4]);%[um]
    y=yf(offset+[0 2 4]);%[um]
    C=Cf(offset+[0 2 4]);
    T=Tf(offset+[0 2 4]);
    PT1=DistansDirectCal(x,y,C,sigma^2,0);%[xp yp]
    TT1=mean(T);
    
    
    x=xf(offset+[0 2 4]+1);
    y=yf(offset+[0 2 4]+1);
    C=Cf(offset+[0 2 4]+1);
    T=Tf(offset+[0 2 4]+1);
    PT2=DistansDirectCal(x,y,C,sigma^2,0);%[xp yp]
    TT2=mean(T);
    
    %%%% Polarization estimation 
    RMCAPD1=RMCAPD1f(offset+[0 1 2 3 4 5]);
    RMCAPD2=RMCAPD2f(offset+[0 1 2 3 4 5]);

    RMCAPD1=sum(RMCAPD1);
    RMCAPD2=sum(RMCAPD2);
    
    if (RMCAPD1+RMCAPD2)>0
    Pol=(RMCAPD1-RMCAPD2)./(RMCAPD1+RMCAPD2);
      
    %%Error assuming poisson statistics 
    DAPD_1=sqrt(RMCAPD1);
    DAPD_1(find(DAPD_1==0))=1;
    DAPD_2=sqrt(RMCAPD2);
    DAPD_2(find(DAPD_2==0))=1;
  
    DPDAPD_1=2*RMCAPD2./((RMCAPD1+RMCAPD2).^2);
    DPDAPD_2=-2*RMCAPD1./((RMCAPD1+RMCAPD2).^2);
  
    DPol=sqrt((DAPD_1.*DPDAPD_1).^2+(DAPD_2.*DPDAPD_2).^2);
    else 
    Pol=nan;
    DPol=nan;
    end
    %%Laser LaserFlag
    LaserFlag=sum(LaserFlagf(offset+[0 1 2 3 4 5]))>0;

    if or(abs(sum(PT1))<inf , abs(sum(PT2))<inf)
       if and(abs(sum(PT1))<inf , abs(sum(PT2))<inf)
	  P=(PT1+PT2)/2;%[xp yp]
	  TimeStamp=(TT1+TT2)/2;
	  %%Window time 
          WindowTime=diff(Tf(offset+[0 5]));
          AcumilatedMeasurementTime= 6*Scanobject.ScanSpeeds.MeasurementOnTime;
       elseif and(abs(sum(PT1))<inf , not(abs(sum(PT2))<inf))
	  P=PT1;%[xp yp]
	  TimeStamp=TT1;
	  %%Window time 
          WindowTime=diff(Tf(offset+[0 4]));
          AcumilatedMeasurementTime= 3*Scanobject.ScanSpeeds.MeasurementOnTime;
       elseif and(not(abs(sum(PT1))<inf) , abs(sum(PT2))<inf)
	  P=PT2;%[xp yp]
	  TimeStamp=TT2;
	  %%Window time 
          WindowTime=diff(Tf(offset+[1 5]));
          AcumilatedMeasurementTime= 3*Scanobject.ScanSpeeds.MeasurementOnTime;
       else
	  disp('Something wrong')	  
       end
    else
       P=[nan; nan];%[xp yp]
       TimeStamp=(TT1+TT2)/2;
       %%Window time 
       WindowTime=diff(Tf(offset+[0 5]));
       AcumilatedMeasurementTime=0;
    end
    
    Para(index,:)=[TimeStamp P' Pol DPol RMCAPD1 RMCAPD2 LaserFlag WindowTime AcumilatedMeasurementTime];
    index=index+1;
  end 
  
  %Adding all information to the Trajectori object
  Trajectories.INFO.TrajPosEstimation.DirectEstimation3Points={[DataSelectionRules.Trajectory.Mode ' mean with a window of 6 elements']; ...
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
                                 
  Trajectories.TrajPosEstimation.DirectEstimation3Points{n}=[Para];
  Trajectories.TrajPosEstimation.DirectEstimation3Points{n}=Trajectories.TrajPosEstimation.DirectEstimation3Points{n}(find(prod(not(isnan(Trajectories.TrajPosEstimation.DirectEstimation3Points{n})),2)),:);% Remove rows with NaN 

end



%% plot part
%[Para1 Para2 sqrt((Para1(:,1)-Para2(:,2)).^2+(Para1(:,2)-Para2(:,3)).^2)]
%  close all
%  figure()
%  scatter3(RawData(:,2),RawData(:,3),RawData(:,1),10*sum(RawData(:,4:5),2)/max(sum(RawData(:,4:5),2)),sum(RawData(:,4:5),2)/max(sum(RawData(:,4:5),2)),'filled')
%  title(['Data overview, number of points:' num2str(size(RawData(:,1),1))])
%  
%  
%  figure()
%  scatter3(RawData(IntervaltoUse,2),RawData(IntervaltoUse,3),RawData(IntervaltoUse,1),10*sum(RawData(IntervaltoUse,4:5),2)/max(sum(RawData(IntervaltoUse,4:5),2)),sum(RawData(IntervaltoUse,4:5),2)/max(sum(RawData(IntervaltoUse,4:5),2)),'filled')
%  hold on 
%  %plot3(RawData(IntervaltoUse,2),RawData(IntervaltoUse,3),RawData(IntervaltoUse,1),'-')
%  
%  plot3(Para1(:,2),Para1(:,3),Para1(:,1),'-g')
%  %plot3(Para2(:,3),Para2(:,4),Para2(:,1),'-r*')
%  
%  plot3(cc1(:,2),cc1(:,3),cc1(:,1),'-k')
%  %plot3(cc2(:,2),cc2(:,3),cc2(:,1),'-m')
%  %plot3(cc3(:,1),cc3(:,2),RawData(3:size(Para1,1)+2,1),'-r*')
%  colorbar
%  
%  figure()
%  scatter(RawData(IntervaltoUse,1),RawData(IntervaltoUse,2),10*sum(RawData(IntervaltoUse,4:5),2)/max(sum(RawData(IntervaltoUse,4:5),2)),sum(RawData(IntervaltoUse,4:5),2)/max(sum(RawData(IntervaltoUse,4:5),2)),'filled')
%  hold on
%  plot(Para1(:,1),Para1(:,2),'-b')
%  %plot(Para2(:,1),Para2(:,3),'-r')
%  plot(cc1(:,1),cc1(:,2),'-r')
%  %plot(cc2(:,1),cc2(:,2),'-m')
%  
%  figure()
%  scatter(RawData(IntervaltoUse,1),RawData(IntervaltoUse,3),10*sum(RawData(IntervaltoUse,4:5),2)/max(sum(RawData(IntervaltoUse,4:5),2)),sum(RawData(IntervaltoUse,4:5),2)/max(sum(RawData(IntervaltoUse,4:5),2)),'filled')
%  hold on 
%  plot(Para1(:,1),Para1(:,3),'-b')
%  %plot(Para2(:,1),Para2(:,4),'-r')
%  plot(cc1(:,1),cc1(:,3),'-r')
%plot(cc2(:,1),cc2(:,3),'-m')