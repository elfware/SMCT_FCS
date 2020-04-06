function Trajobject=getTrajPosEstimationPolyfit(Trajobject,DataSelectionRules,nn);


for FrajIndex=nn
%TrajData=Trajobject.TrajPosEstimation.DirectEstimation3Points{FrajIndex};
TrajData=Trajobject.TrajPosEstimation.FixWindowAv{FrajIndex};
Time=TrajData(:,1);
X=TrajData(:,2);
Y=TrajData(:,3);
Pol=TrajData(:,4);
round((Time(end)-Time(1))*8)

ppX = splinefit (Time, X, round((Time(end)-Time(1))*7.6), "order", 6, "beta", 0.25);
ppY = splinefit (Time, Y, round((Time(end)-Time(1))*7.6), "order", 6, "beta", 0.25);

%% Plot
XFit = ppval (ppX, Time);
YFit = ppval (ppY, Time);

figure()
plot(Time, X,'.')
hold on
plot(Time, XFit,'-')

figure()
plot(Time, Y,'.r')
hold on
plot(Time, YFit,'-r')

figure()
plot3(X,Y,Time,'.')
hold on
plot3(XFit, YFit, Time,'-r')

end
