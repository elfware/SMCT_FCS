function Para=DistansDirectCal(x,y,C,varinace,m)
%% Assuming that the fluorophore measured intensity is
%% I=(Io\sigma^2)*e^(-((x-xo)^2+(y-yo)^2)/2*sigma^2)
%% Taking the ratio of two measurements one no longer depends on Io
%%
%% If m=0:
%% estimation of position only (output is [xp yp]) 
%% Input must be three points long, if longer it will only use the first three
%%
%% If m=1:
%% estimation of position and varinace (output is [variance xp yp]) 
%% the variance input is in this case only a dummy.
%% Input must be four points long, if longer it will only use the first four  

if m==0
    A=[x(2)-x(1) y(2)-y(1);...
       x(3)-x(1) y(3)-y(1)];

    LC=log(C);   
    B=0.5*[x(2).^2-x(1).^2+y(2).^2-y(1).^2-2*varinace*(LC(1)-LC(2));...
	   x(3).^2-x(1).^2+y(3).^2-y(1).^2-2*varinace*(LC(1)-LC(3))];
       
    Para=A\B;%[xp yp]
    
elseif m==1
    LC=log(C);   
    A=[LC(1)-LC(2) x(2)-x(1) y(2)-y(1);...
       LC(1)-LC(3) x(3)-x(1) y(3)-y(1);...
       LC(1)-LC(4) x(4)-x(1) y(4)-y(1)];

    B=0.5*[x(2).^2-x(1).^2+y(2).^2-y(1).^2;...
	   x(3).^2-x(1).^2+y(3).^2-y(1).^2;...
	   x(4).^2-x(1).^2+y(4).^2-y(1).^2];
       
    Para=A\B;%[variance xp yp]
end