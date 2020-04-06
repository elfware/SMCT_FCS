function [Scanobject]=getScanData(Scanobject)

% Retrive data 
StringFind=strfind(Scanobject.DataPath.pathname,'/');
Scanobject.DataPath.pathname=Scanobject.DataPath.pathname(1:StringFind(end));
fid=fopen([Scanobject.DataPath.pathname Scanobject.DataPath.filenameScanData]);
Data = fread(fid, 'uint64=>uint64','ieee-be');
fclose(fid);

%_____________________________________
% Extract data from 64 bit number
% F= Flag 
% |FF|    APD_1    |FF|    APD_2    |      x-data    |    y-data      |
%  xxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxx
ySG=double(bitand(Data,2^16-1)); 
xSG=double(bitand(bitshift(Data,-16),2^16-1));
APD_2=double(bitand(bitshift(Data,-32),2^14-1));     % V polarization  
TrakingFlag(:,1)=bitand(bitshift(Data,-32),2^15)~=0; %Tracking on
TrakingFlag(:,3)=bitand(bitshift(Data,-32),2^14)~=0; %
APD_1=double(bitand(bitshift(Data,-48),2^14-1));     % H polarization  
TrakingFlag(:,2)=bitand(bitshift(Data,-48),2^15)~=0; %Releas Threshold not reached 
TrakingFlag(:,4)=bitand(bitshift(Data,-48),2^14)~=0; %Laser Switch during tracking


x=Scanobject.PosCalibration.ScaleFactor(1)*Scanobject.PosCalibration.PositionLength*(xSG-Scanobject.PosCalibration.Offset(1))/(Scanobject.PosCalibration.Max(1)-Scanobject.PosCalibration.Offset(1));
y=Scanobject.PosCalibration.ScaleFactor(2)*Scanobject.PosCalibration.PositionLength*(ySG-Scanobject.PosCalibration.Offset(2))/(Scanobject.PosCalibration.Max(2)-Scanobject.PosCalibration.Offset(2));

ScanSpeedIndex=(TrakingFlag(:,1)|TrakingFlag(:,2)|TrakingFlag(:,3)|TrakingFlag(:,4))/Scanobject.ScanSpeeds.ScanSpeedTracking + not((TrakingFlag(:,1)|TrakingFlag(:,2)|TrakingFlag(:,3)|TrakingFlag(:,4)))/Scanobject.ScanSpeeds.ScanSpeed;
ScanTime=cumsum(ScanSpeedIndex);

Scanobject.ScanData.INFO={'1.Time[s]'; '2.x[um]'; '3.y[um]'; '4.APD1[counts]'; '5.APD2[counts]'; '6.TrackingON'; '7.ReleasThresholdOK'; '8.NotUsed'; '9.LaserSwitchMode' };
Scanobject.ScanData.ScanMatrix=[ScanTime x y APD_1 APD_2 TrakingFlag(:,1) TrakingFlag(:,2) TrakingFlag(:,3) TrakingFlag(:,4)];


%%% 
%%% Shift and remove exece elements after TrackingOnSwitch switches off
%%% This is done to conforme to the TTTR Data set.  
APDShift=Scanobject.PosCalibration.APDShift;   
Scanobject.ScanData.ScanMatrix=Scanobject.ScanData.ScanMatrix(1+APDShift:end,:);
ReleasThresholdOK=Scanobject.ScanData.ScanMatrix(:,7);
TrackingOnSwitch=[0;diff(ReleasThresholdOK)];

V=ones(size(Scanobject.ScanData.ScanMatrix,1),1);
V([find(TrackingOnSwitch==-1)+1; find(TrackingOnSwitch==-1)+2])=0;
V=repmat(V,1,size(Scanobject.ScanData.ScanMatrix,2));

Scanobject.ScanData.ScanMatrix=reshape(Scanobject.ScanData.ScanMatrix(find(V)),size(find(V),1)/size(Scanobject.ScanData.ScanMatrix,2),size(Scanobject.ScanData.ScanMatrix,2));


%%%
