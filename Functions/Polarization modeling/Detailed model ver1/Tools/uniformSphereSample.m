function [my] = uniformSphereSample(phiCenter,thetaCenter,Fp)
%This function returns a uniform sample from sphere surface, with "fraction
%of sampling" from surfacec gven på Fp (Fp=1 whole sphere surface)
%The patch is centered at (phicenter,thetacenter)
center=[cos(phiCenter)*sin(thetaCenter);sin(phiCenter)*sin(thetaCenter);cos(thetaCenter)];
cosNegPhiCenter=cos(-phiCenter);
sinNegPhiCenter=sin(-phiCenter);
u=rand(); %random number between 0 and 1
v=(1-Fp)+rand()*Fp; % rnadom numberbetween 1-Fp and 1

phi_patch=2*pi*u;%phi_patch between [0 2pi] uniform over patch area
%2*v-1 should go from [1-2*Fp to 1] since the surfface area of sampled area of the sampled patch is is 4*pi*r^2*Fp which gives theta_P_max = acos(1-2*FP) and theta_P_min=0=acos(1)
theta_patch=acos(2*v-1);  %thetapatch between [0 acos(1-2*Fp] over the patch area. this Interval is consistent with that the patch area = Fp/total surface area of spere
%coordinates in the patch coordiante system
x_patch=cos(phi_patch)*sin(theta_patch);
y_patch=sin(phi_patch)*sin(theta_patch);
z_patch=cos(theta_patch);

baseVec3=center;

%Rotate down center to phi=0 (around z-axis) to get consistent
%definition of baseVec1 and baseVec2.
RotZ=[cosNegPhiCenter -sinNegPhiCenter 0;sinNegPhiCenter cosNegPhiCenter 0;0 0 1];
baseVec3rot=RotZ*baseVec3;
%After the rotation the y axis is always orthogonal to baseVec3rot
baseVec2rot=[0;1;0];
%The third base vector is the cross product of the two already found
baseVec1rot=cross(baseVec3rot,baseVec2rot);
%Rotate back the base vectors
RotZback=[cosNegPhiCenter sinNegPhiCenter 0;-sinNegPhiCenter cosNegPhiCenter 0;0 0 1];
baseVec1=RotZback*baseVec1rot;
baseVec2=RotZback*baseVec2rot;
%get cartesian coordinates in normal "sphere" coordiante system
my = baseVec1*x_patch + baseVec2*y_patch + baseVec3*z_patch;


end

