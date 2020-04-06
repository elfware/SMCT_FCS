%function [R k]=BinaryAoutoCorr_V3(TTTRPos,points,m)
clear all
close all
TTTRPos=sort(round(rand(1e5,1)*1e6))*5e-9;
points=1000;
m=1;

%m=1; %Binning  amount m=1 is one step or no binning => resolution 5ns, m=2 is two step binning of 2 => resolution 10ns, m*5 ns is the resolution for binningn size m  
Binary_pointer=round((TTTRPos-TTTRPos(1))/5e-9)+1;
N1=size(Binary_pointer,1); %Amount of ones 
N=Binary_pointer(end);%last pointer which is also the total amount of enterise zeros and ones


%%xcorr test
tic
photonDens=1/diff((Binary_pointer-1)*5e-9);
photonDensTime=(Binary_pointer(1:end-1)-1)*5e-9;
[Rxcorr Lag]=xcorr(photonDens,points,'biased');
toc
figure()
plot(photonDensTime(Lag(points+1:end)+1),Rxcorr(points+1:end)/Rxcorr(points+1))
%%

n=(N/m); %Amount of binns
MeanBinary=N1/n;


if points>n-1
disp(['To large displacment restrict to max points: ' num2str(n)])
points=n;
end

OFFset=[-(m-1):m-1];
ScaleFactor= 10^(ceil(log10(size(OFFset,2)+1)));

tic
for k=0:points;
   [ANoOffset BNoOffset]=find(Binary_pointer(find(Binary_pointer<=N-k*m))<=N-k*m-OFFset);%0.13

   %BNoOffset; %Which m is used

   Binary_pointerNoOffset=Binary_pointer( ANoOffset ); %0.016

   %[Binary_pointerNoOffset BNoOffset];
   [AOffset BOffset]=find(Binary_pointer(find(Binary_pointer>=k*m+1))>=k*m+1+OFFset);%0.13

   %BOffset; %Which m is used
   if m==1
    Binary_pointerOffset=Binary_pointer(  AOffset  )-k*m-OFFset(BOffset);
   else
    Binary_pointerOffset=Binary_pointer(  AOffset  )-k*m-OFFset(BOffset)';
   end
   
   %[Binary_pointerOffset BOffset];
   
   %[[Binary_pointerNoOffset; Binary_pointerOffset] [BNoOffset; BOffset]];
   CIndex=[BNoOffset; BOffset]; %indexed m values 
   CPointers=[Binary_pointerNoOffset; Binary_pointerOffset];
   
   [A B]=sort(CPointers*ScaleFactor+CIndex);%1.16
   D=find(diff(A)==0); % Pointer in sorted list to all elements that creates a 1 in a product 
   CPointers2=CPointers(B(D));% Pointer to all elements that creates a 1 in a product 
   SMdiagP=sum(mod(CPointers2-1,m)<=(m-1)-abs(CIndex(B(D))-m));
   MdiagP=SMdiagP;

%     [A B]=find( (Binary_pointerNoOffset*ScaleFactor+BNoOffset)==(Binary_pointerOffset*ScaleFactor+BOffset)');% Pointer in list to all elements that creates a 1 in a product  
%     CPointers2=Binary_pointerNoOffset(A);% Pointer to all elements that creates a 1 in a product 
%     SMdiagP=sum(mod(CPointers2-1,m)<=(m-1)-abs(BNoOffset(A)-m));
%     MdiagP=SMdiagP;

      
%    MdiagP=0;
%    for OFFset=[-(m-1):(m-1)]
%      Binary_pointerNoOffset=Binary_pointer(find((Binary_pointer)<=N-k*m-OFFset));
%      Binary_pointerOffset=Binary_pointer(find((Binary_pointer-k*m-OFFset)>=1))-k*m-OFFset;
%      
%      SMdiagP=sum(mod(find(diff(sort([Binary_pointerNoOffset; Binary_pointerOffset]))==0)-1,m)<=(m-1)-abs(OFFset));
%    
%      MdiagP=MdiagP+SMdiagP;
%    end
  
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
toc

%%output
k=[0:points]*m*5e-9;
R;

hold on
plot(k,R,'-*r')