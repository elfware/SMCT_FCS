function R=BinaryAoutoCorr(TTTRPos,points,Res)

Binary_pointer=round((TTTRPos-TTTRPos(1))/Res)+1;
N1=size(Binary_pointer,1);
N=Binary_pointer(end);

MeanBinary_pointer=N1/N;
VarBinary_pointer=(N*MeanBinary_pointer^2+N1*(1-2*MeanBinary_pointer))/(N-1);

for k=1:points;
  Binary_pointerNoOffset=Binary_pointer(find((Binary_pointer)<=N-k));
  Binary_pointerOffset=Binary_pointer(find((Binary_pointer-k)>=1))-k;
  SXtXtpk=sum(diff(sort([Binary_pointerNoOffset; Binary_pointerOffset]))==0);
  R(k)=(MeanBinary_pointer^2/VarBinary_pointer) + (SXtXtpk-MeanBinary_pointer*( sum(Binary_pointer<=N-k) + sum(Binary_pointer>=k+1) ))/((N-k)*VarBinary_pointer);
end
