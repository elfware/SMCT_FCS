function [I,aa,bb,cc,ab,ac,bc] = calcPolIntsNA(theta,phi,epsilon)
% returns the polarizatiion intensities I = [Ix;Iy;Iz] =[IH,IV,Iz]
%With numerical aperture defined by the angle epsilon, theta is te
%azimuthal angle and phi is the inclenation angle

x=cos(phi)*cos(theta);
y=cos(phi)*sin(theta);
z=sin(phi);

sce=sin(phi+epsilon)*cos(phi+epsilon);
sc=sin(phi)*cos(phi);

%From integration incresing phi to phi+epsilon
k1 = (1/2)*(epsilon+sce-sc);
k2=(1/2)*(epsilon-sce+sc);
k3=(1/2)*(cos(phi)^2-cos(phi+epsilon)^2);

aa=zeros(1,3);
bb=zeros(1,3);
cc=zeros(1,3);
ab=zeros(1,3);
ac=zeros(1,3);
bc=zeros(1,3);
% integrate around fluorophore axis, Ix coeeficients
aa(1)=pi*(3*x^4-2*x^2+1);
bb(1)=pi*(z^2+3*x^2*y^2);
cc(1)=pi*(y^2+3*x^2*z^2);
ab(1)=pi*(3*x^2-1)*x*y;
ac(1)=pi*(3*x^2-1)*x*z;
bc(1)=pi*(3*x^2-1)*y*z;

%Iy coefficients
aa(2)=pi*(z^2+3*x^2*y^2);
bb(2)=pi*(3*y^4-2*y^2+1);
cc(2)=pi*(x^2+3*z^2*y^2);
ab(2)=pi*(3*y^2-1)*y*x;
ac(2)=pi*(3*y^2-1)*x*z;
bc(2)=pi*(3*y^2-1)*y*z;

%Iz cpefficients
aa(3)=pi*(y^2+3*x^2*z^2);
bb(3)=pi*(x^2+3*y^2*z^2);
cc(3)=pi*(3*z^4-2*z^2+1);
ab(3)=pi*(3*z^2-1)*x*y;
ac(3)=pi*(3*z^2-1)*x*z;
bc(3)=pi*(3*z^2-1)*y*z;

I=k1*(aa*cos(theta)^2+2*ab*cos(theta)*sin(theta)+bb*sin(theta)^2)+k2*cc+k3*(2*ac*cos(theta)+2*bc*sin(theta));
%divide by area that was integrated over (area of numerical aparture on
%sphere surface)
I=I/(4*pi*sin(epsilon*2));
end

