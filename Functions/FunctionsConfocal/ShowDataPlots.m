function A=ShowDataPlots(Scanobjects,Trajobject,DataSelectionRulesNew,n,PlotString)
  
for Index=[1:size(Trajobject.INFO.Traj,1)]     
    if not(isempty(strfind(Trajobject.INFO.Traj{Index},'Scan file index')))
       dataplace=str2double(Trajobject.INFO.Traj{Index}(1));
    end
end

  %dataplace=str2num(Trajobject.INFO.Traj'(strfind(Trajobject.INFO.Traj','Scan file index')-2));
  Scanobject=Scanobjects{Trajobject.Traj(n,dataplace)}.Scanobject;
  DataSelectionRules=Scanobjects{Trajobject.Traj(n,dataplace)}.DataSelectionRules;
  DataSelectionRulesNew;
  %DataSelectionRules need to be updated with the enteries of DataSelectionRulesNew
  
%    ShowTraj(Scanobject,Trajobject,DataSelectionRules,n);
%    ShowTraj(Scanobject,Trajobject_RollingAv,DataSelectionRules,n);
  
  
  RawData=Scanobject.ScanData.ScanMatrix(Trajobject.Traj(n,2):Trajobject.Traj(n,3),1:5);
  RawDataS=[RawData(:,1:3) sum(RawData(:,4:5),2)];
  
  figure()
  subplot(1,2,1)
  scatter(RawDataS(:,1),RawDataS(:,2),1e-5+50*RawDataS(:,4)/max(RawDataS(:,4)),1e-5+RawDataS(:,4)/max(RawDataS(:,4)),'filled')
  hold on 
  for i=size(fieldnames(Trajobject.TrajPosEstimation),1)
    FieldName=fieldnames(Trajobject.TrajPosEstimation);
    eval(['TrajData=Trajobject.TrajPosEstimation.' FieldName{i} ';']);
    plot(TrajData{n}(:,1),TrajData{n}(:,2),'b')
  end
  colorbar
%  hold Off
  grid on 
  xlabel('time [s]')
  ylabel('x-pos [um]')
  
  subplot(1,2,2)
  scatter(RawDataS(:,1),RawDataS(:,3),1e-5+50*RawDataS(:,4)/max(RawDataS(:,4)),1e-5+RawDataS(:,4)/max(RawDataS(:,4)),'filled')
  hold on 
  for i=size(fieldnames(Trajobject.TrajPosEstimation),1)
      FieldName=fieldnames(Trajobject.TrajPosEstimation);
    eval(['TrajData=Trajobject.TrajPosEstimation.' FieldName{i} ';']);
    plot(TrajData{n}(:,1),TrajData{n}(:,3),'b')
  end
  colorbar
%  hold Off
  grid on 
  title(['Raw data, Max value' num2str(max(RawDataS(:,4)))])
  xlabel('time [s]')
  ylabel('y-pos [um]')
  
  
  
  ShowTraj(Scanobject,Trajobject,DataSelectionRules, n,PlotString) 

for Index=[1:size(Trajobject.INFO.Traj,1)]     
    if not(isempty(strfind(Trajobject.INFO.Traj{Index},'Diffusion coefficent x [um^2/s]')))
       xpos=str2double(Trajobject.INFO.Traj{Index}(1));
    end
end
%[A B]=ind2sub(size(Trajobject.INFO.Traj'),strfind(Trajobject.INFO.Traj','Diffusion coefficent x [um^2/s]')-2);
%xpos=str2num(Trajobject.INFO.Traj'(1:A,B)');

for Index=[1:size(Trajobject.INFO.Traj,1)]     
    if not(isempty(strfind(Trajobject.INFO.Traj{Index},'Diffusion coefficent y [um^2/s]')))
       ypos=str2double(Trajobject.INFO.Traj{Index}(1));
    end
end
%[A B]=ind2sub(size(Trajobject.INFO.Traj'),strfind(Trajobject.INFO.Traj','Diffusion coefficent y [um^2/s]')-2);
%ypos=str2num(Trajobject.INFO.Traj'(1:A,B)');

disp(['Dx ' num2str(Trajobject.Traj(n,xpos)) ' [um^2/s]' ' Dy ' num2str(Trajobject.Traj(n,ypos)) ' [um^2/s]']);


