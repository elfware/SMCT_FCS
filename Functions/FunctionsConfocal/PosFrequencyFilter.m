function Trajobject=PosFrequencyFilter(Trajobject,DataSelectionRules,Cutoff)

eval(['TrajData=Trajobject.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode ';'])

for Trajindex=[1:size(Trajobject.Traj,1)]

%%% FFT position filtering 
x=TrajData{Trajindex}(:,2);
y=TrajData{Trajindex}(:,3);
t=TrajData{Trajindex}(:,1);

fftX=fft(x);
fftY=fft(y);

fftXfilterpar=fftshift(fftX);
fftYfilterpar=fftshift(fftY);

f_s=1/mean(diff(t));
N=size(t,1);
F=f_s/N*[-N/2:N/2-1];

SetTo0=abs(F)>Cutoff;
fftXfilterpar(SetTo0)=0;
fftYfilterpar(SetTo0)=0;

Xnew=ifft(fftshift(fftXfilterpar));
Ynew=ifft(fftshift(fftYfilterpar));

TrajData{Trajindex}(:,2)=real(Xnew);
TrajData{Trajindex}(:,3)=real(Ynew);

eval(['Trajobject.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode 'FFTFilterd{Trajindex}=Trajobject.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode '{Trajindex};']);
eval(['Trajobject.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode 'FFTFilterd{Trajindex}=TrajData{Trajindex};']);
end

  %Adding all information to the Trajectori object
 
  INFOText=[DataSelectionRules.Trajectory.Mode 'FFTFilterd is position filtered with a cutoff at ' num2str(Cutoff) ' Mhz'; ...
            '1.' DataSelectionRules.Trajectory.Mode ' mean time [s]'; ...
            '2.' DataSelectionRules.Trajectory.Mode ' mean x '; ...
            '3.' DataSelectionRules.Trajectory.Mode ' mean y'; ...
            '4.' DataSelectionRules.Trajectory.Mode ' mean polarization'; ...
            '5.' DataSelectionRules.Trajectory.Mode ' mean poisson polarization error'; ...
            '6.APD1 counts over the' DataSelectionRules.Trajectory.Mode '  window'; ...
            '7.APD2 counts over the' DataSelectionRules.Trajectory.Mode '  window'; ...
            '8.Laser high Flag'; ...
            '9.' DataSelectionRules.Trajectory.Mode ' window time [s]';...
            '10.' DataSelectionRules.Trajectory.Mode ' acumilated measurement time [s]'];
  
  eval(['Trajobject.INFO.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode 'FFTFilterd=INFOText;']);
