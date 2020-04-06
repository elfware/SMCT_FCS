function [my] = myUnpolExciation(PHI,THETA)
% Returns The direction dipole moment mu of a fluorophore, scalde by
% exciation of unpolarized light (smaller in z-direction). Physics
% definition of spherical coordinated PHI and THETA according to BAcker and
% Moerner 2014.

st=sin(THETA);
sp=sin(PHI);

ct=cos(THETA);
cp=cos(PHI);

my =st*[st*cp;st*sp;ct];

end

