function [G] = greenTensorBFP_v2(phi,theta)
%This function returns the greens tenors at the back focal plane of the objective
%for a dipole as definedin Backer
%and Moerner 2014. phi is the azimuthal angle and theta is the polar angle
%defining where you are looking o the electric filed (where on the torus)
%in the same coordiante system as the dipole.

% Change log 
% 2016 11 24: The function cna now take phi and theta as vectors (Elias)

L=size(phi,2);

sp= sin(phi);
rho= sin(theta);
cp= cos(phi);
%ct= cos(theta);

s2p=sin(2*phi);
spSquared=sp.^2;
cpSquared=cp.^2;
rhoSquared=rho.^2;
fctr1=sqrt(1-rhoSquared);
fctr2=s2p.*(fctr1 -1)./2;
%  G =[spSquared + cpSquared.*fctr1, fctr2, -rho.*cp;
%      fctr2, cpSquared+spSquared.*fctr1, -rho.*sp;
%      0,0,0].*sqrt(1./fctr1);


G11=reshape( (spSquared + cpSquared.*fctr1).*sqrt(1./fctr1),1,1,L);
G12=reshape( (fctr2).*sqrt(1./fctr1),1,1,L);
G13=reshape( (-rho.*cp).*sqrt(1./fctr1),1,1,L);
G21=reshape( (fctr2).*sqrt(1./fctr1),1,1,L);
G22=reshape( (cpSquared+spSquared.*fctr1).*sqrt(1./fctr1),1,1,L);
G23=reshape( (-rho.*sp).*sqrt(1./fctr1),1,1,L);
G3a=zeros(1,1,L);    

G=[G11 G12 G13;...
   G21 G22 G23;...
   G3a G3a G3a];
end
