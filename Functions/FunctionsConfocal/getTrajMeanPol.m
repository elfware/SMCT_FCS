function Trajectories=getTrajMeanPol(Scanobject,Trajectories,p,PolPlotIndexThrechold)

m=size(Trajectories.Traj,2);
nn=1:size(Trajectories.Traj,1);
for n=nn; 
Interval=Trajectories.Traj(n,2):Trajectories.Traj(n,3);

  APD1=sum(Scanobject.ScanData.ScanMatrix(Interval,4));
  APD2=sum(Scanobject.ScanData.ScanMatrix(Interval,5));

  %%%%Polarization estimation 

  Pol=(APD1-APD2)./(APD1+APD2);
  
  %%Error assuming poisson statistics 
  DAPD_1=sqrt(APD1);
  DAPD_1(find(DAPD_1==0))=1;
  DAPD_2=sqrt(APD2);
  DAPD_2(find(DAPD_2==0))=1;
  
  DPDAPD_1=2*APD2./((APD1+APD2).^2);
  DPDAPD_2=-2*APD1./((APD1+APD2).^2);
  
  DPol=sqrt((DAPD_1.*DPDAPD_1).^2+(DAPD_2.*DPDAPD_2).^2);

  %Adding all information to the Trajectori object                                   
  Trajectories.Traj(n,m+[1:4])=[Pol DPol APD1 APD2];                                 
  %Trajectories.TrajMeanPol(n,:)=[Pol DPol APD1 APD2];
end
  %Adding INfO 
  Trajectories.INFO.Traj=[Trajectories.INFO.Traj;...
			  num2str(m+1) '.Mean polarization'; ...
                          num2str(m+2) '.Poisson polarization error'; ...
                          num2str(m+3) '.APD1 counts'; ...
                          num2str(m+4) '.APD2 counts'];


if p
Index=find(abs(Trajectories.TrajMeanPol(:,1))>PolPlotIndexThrechold);

figure()
hist(Trajectories.TrajMeanPol(:,1),100)
xlabel('Polarization')
ylabel('Occurence')
hold on 
h=1;
for I=1:size(Index,1) 
text(Trajectories.TrajMeanPol(Index(I),1),h,num2str(Index(I)))
h=h+0.1;
end
hold off
end

