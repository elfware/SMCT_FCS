function Scanobjects=GetTTTRData(Scanobjects,Trajobject,V)
%%If V=1 the Loding info is shown

for Index=[1:size(Trajobject.INFO.Traj,1)]     
    if not(isempty(strfind(Trajobject.INFO.Traj{Index},'Scan file index')))
       FileIndexOfInterest=str2double(Trajobject.INFO.Traj{Index}(1));
    end
end
FilesInUse=Trajobject.Traj([1 find(diff(sort(Trajobject.Traj(:,FileIndexOfInterest))))'+1],FileIndexOfInterest)';

  
for Index=FilesInUse
  Scanobject=Scanobjects{Index}.Scanobject;
  [ EventStructure ] = getTTTR_FPGA_DataO(Scanobject, 0, V);
  Scanobjects{Index}.Scanobject.ScanData.TTTREventStructure=EventStructure;
end
