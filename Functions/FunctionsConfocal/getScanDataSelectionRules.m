function Scanobject=getScanDataSelectionRules(Scanobject,DataSelectionRules);

APDThreshold=DataSelectionRules.APDThreshold;
CountRange=DataSelectionRules.PhotonCounts.CountRange;
ErrorBound=DataSelectionRules.PhotonCounts.ErrorBound;
xBounds=DataSelectionRules.PositionTime.xBounds;
yBounds=DataSelectionRules.PositionTime.yBounds;
TBounds=DataSelectionRules.PositionTime.TBounds;

TotData=Scanobject.ScanData.ScanMatrix(:,1:5);


%Select Interesting values
  % First row:  Select on min value for each APD
  % Second row: Select by minimum total counts and error size
  % Thired row: Select over positions and time 

%  %%%              1              2               3                  4                    5        
%  TotData=[ScanTime(1:size(x,1)) x(1:size(x,1)) y(1:size(x,1)) APD_1(1:size(x,1)) APD_2(1:size(x,1))];
 
Index=find( (TotData(:,4)>=APDThreshold(1)|TotData(:,5)>=APDThreshold(2)) & ... 
((TotData(:,4)+TotData(:,5))>CountRange(1)) & ((TotData(:,4)+TotData(:,5))<CountRange(2)) & ...
((TotData(:,2)>xBounds(1))&(TotData(:,2)<xBounds(2))) & ...
((TotData(:,3)>yBounds(1))&(TotData(:,3)<yBounds(2))) & ...
((TotData(:,1)>TBounds(1))&(TotData(:,1)<TBounds(2))));                


[A B]=size(Scanobject.ScanData.ScanMatrix);
NewColumn=zeros(A,1);
NewColumn(Index)=1;

Scanobject.ScanData.ScanMatrix(:,B+1)=NewColumn;
Scanobject.ScanData.INFO=[Scanobject.ScanData.INFO; num2str(B+1) '.PasedDataSelectionRules'];