function [G] = greenTensorFF(phi,theta)
%This function returns the farfield greens tenor of a dipole as definedin Backer
%and Moerner 2014. phi is the azimuthal angle and theta is the polar angle
%defining where you are looking o the electric filed (where on the torus)
%in the same coordiante system as the dipole.
sp= sin(phi);
st= sin(theta);
cp= cos(phi);
ct= cos(theta);

spSquared=sp^2;
stSquared=st^2;
cpSquared=cp^2;
ctSquared=ct^2;

fctr1=-cp*sp*stSquared;
fctr2=ct*st;
fctr3=-sp*fctr2;
fctr4=-cp*fctr2;

G=[1-cpSquared*stSquared fctr1 fctr4;
    fctr1 1-spSquared*stSquared fctr3;
    fctr4 fctr3 1-ctSquared];

end

