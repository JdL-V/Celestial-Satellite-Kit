function [alpha,delta]=car2Dec_Ra(r)

% car2Dec_Ra.m-transform cartesian coordinates in right ascension and declination.
%
% PROTOTYPE:
% [alpha,delta]=car2Dec_Ra(r)
%
% INPUT:
% r         [:x3] satellite position during all the time span              [km]
%
% OUTPUT:
% alpha   [1x1] right ascension [rad]
% delta   [1x1] declination     [rad]

xx = r(:,1);
yy = r(:,2);
zz = r(:,3);
rr = vecnorm(r')';
l = xx./rr;
m = yy./rr;
n = zz./rr;
s = max(size(xx));
delta = asin(n);
alpha = zeros(s,1);
alpha = acos(l./cos(delta));
for i = 1:s
    if(m(i) <= 0)
        alpha(i) = 2*pi-alpha(i);
    end
end
end