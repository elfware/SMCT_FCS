function [frac] = theoryFractionBound(Kd,Ptot,Dtot)
%% Returns fraction bound ([PD]/([D]+[PD]))for the equilibrium [PD] <--> [P]+[D] (Kd)

thesum=Kd+Ptot+Dtot;

frac=sqrt((thesum/2).^2+Ptot.*Dtot)-thesum/2;
frac=frac./Dtot;



end

