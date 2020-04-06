function [G] = greenTensorBFP(phi,theta)
%This function returns the greens tenors at the back focal plane of the objective
%for a dipole as definedin Backer
%and Moerner 2014. phi is the azimuthal angle and theta is the polar angle
%defining where you are looking o the electric filed (where on the torus)
%in the same coordiante system as the dipole.
sp= sin(phi);
rho= sin(theta);
cp= cos(phi);
%ct= cos(theta);

s2p=sin(2*phi);
spSquared=sp^2;
cpSquared=cp^2;
rhoSquared=rho^2;
fctr1=sqrt(1-rhoSquared);
fctr2=s2p*(fctr1 -1)/2;
G =[spSquared + cpSquared*fctr1, fctr2, -rho*cp;
    fctr2, cpSquared+spSquared*fctr1, -rho*sp;
    0,0,0]*sqrt(1/fctr1);

end

