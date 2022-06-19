function [a,e,i,OM,om,th] = car2kep(rr,vv,mu)

% car2kep.m-Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a,e,i,OM,om,th]=car2kep(r,v,mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in 
% radians.
%
%INPUT:
%r         [3x1]      Position vector             [km]
%v         [3x1]      Velocity vector             [km/s]
%mu        [1x1]      Gravitational parameter     [km^3/s^2]
%
% OUTPUT:
% a        [1x1]      Semi-major axis      [km]
% e        [1x1]      Eccentricity         [-]
% i        [1x1]      Inclination          [rad]
% OM       [1x1]      RAAN                 [rad]
% om       [1x1]      Pericentre anomaly   [rad]
% th       [1x1]      True anomaly         [rad]

k=[0;0;1];

% Setp 1:

r=norm(rr);      
v=norm(vv); 

% Sep 2:

E=1/2*v^2-mu/r; % Energy
a=-mu/(2*E);    % Semi-major axis

% Step 3:

ee=1/mu.*((v^2-mu/r).*rr-(dot(rr,vv).*vv)); % vector eccentricity
e=norm(ee); % eccentricity
if(abs(ee(1))<0.00000001)
    ee(1)=0;
end
if(abs(ee(2))<0.00000001)
    ee(2)=0;
end
if(abs(ee(3))<0.00000001)
    ee(3)=0;
end
if(abs(e)<0.00000001)
    e=0;
end
 % Step 4:
 
 hh=cross(rr,vv); % vector momentum
 h=norm(hh);
 
 % Step 5:
 
 i=acos(hh(3)/h); % inclination
 if(abs(i)<0.00000001)
     i=0;
 end
 
 % Step 6:
 
 nn=cross(k,hh);
 n=norm(nn);
 if(n<0.000001)
     n=0;
 end
 
 % Step 7:
 
 if(nn(2)>=0)
     OM=acos(nn(1)/n);
 else
     OM=2*pi-acos(nn(1)/n);
 end
 if(n==0)
     OM=0;
 end
 
 % Step 8:
 
 if(ee(3)>=0)
     om=acos(dot(nn,ee)/(n*e));
 else
     om=2*pi-acos(dot(nn,ee)/(n*e));
 end
 if(e==0)
     om=0;
 end
 
 % Step 9:
 
 vr=dot(rr,vv)/r;
 if(vr>=0)
     th=acos(dot(ee,rr)/(e*r));
 else
     th=2*pi-acos(dot(ee,rr)/(e*r));
 end
 if(e==0)
     th=0;
 end
end