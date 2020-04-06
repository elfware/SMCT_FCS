function [R k]=BinaryAoutoCorr_V2(TTTRPos,points,m)

%m=1; %Binning  amount m=1 is one step or no binning => resolution 5ns, m=2 is two step binning of 2 => resolution 10ns, m*5 ns is the resolution for binningn size m  
Binary_pointer=round((TTTRPos-TTTRPos(1))/5e-9)+1;
N1=size(Binary_pointer,1); %Amount of ones 
N=Binary_pointer(end);%last pointer which is also the total amount of enterise zeros and ones

n=(N/m); %Amount of binns
MeanBinary=N1/n;


if points>n-1
disp(['To large displacment restrict to max points: ' num2str(n)])
points=n;
end

for k=0:points;
  MdiagP=0;
  for OFFset=[-(m-1):(m-1)]
    Binary_pointerNoOffset=Binary_pointer(find((Binary_pointer)<=N-k*m-OFFset));
  
    Binary_pointerOffset=Binary_pointer(find((Binary_pointer-k*m-OFFset)>=1))-k*m-OFFset;
    SMdiagP=sum(mod(find(diff(sort([Binary_pointerNoOffset; Binary_pointerOffset]))==0)-1,m)<=(m-1)-abs(OFFset));
  
    MdiagP=MdiagP+SMdiagP;
  end
  
  %%Calculation of variance for binned data
  if k==0
    VarBinMDdiag=MdiagP;
    VarBinary=abs((VarBinMDdiag/n)-MeanBinary.^2);
  end
    
  %% Calculate autocorrelation  
  R(k+1)=(MeanBinary^2/VarBinary) -...
       (MeanBinary/VarBinary)*( sum(Binary_pointer<=N-m*k)/(n-k) + sum(Binary_pointer>=m*k+1)/(n-k)) +...
       (1/VarBinary)*MdiagP/(n-k);
end

%%output
k=[0:points]*m*5e-9;
R;