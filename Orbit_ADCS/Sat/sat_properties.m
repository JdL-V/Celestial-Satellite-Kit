function [I0] = sat_properties()
% Definition of moments of inertia
A = 1.578e-4;
B = 3.131e-4;
C = 3.093e-4;

I0 = [A 0 0;
      0 B 0;
      0 0 C];
end

