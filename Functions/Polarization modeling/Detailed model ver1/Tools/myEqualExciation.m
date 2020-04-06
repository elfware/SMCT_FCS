function [my] = myEqualExciation(PHI,THETA)
% Returns The direction dipole moment mu of a fluorophore, same excöiation in x y and z. Physics
% definition of spherical coordinated PHI and THETA according to BAcker and
% Moerner 2014.

st=sin(THETA);
sp=sin(PHI);

ct=cos(THETA);
cp=cos(PHI);

my =[st.*cp; st.*sp; ct];

end

