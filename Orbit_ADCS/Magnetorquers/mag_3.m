function [I, mu] = mag_2(dt, N, A, J, BBx, B)
    N = [250; 850; 850];
    A = [33^2; 11^2; 11^2]*1e6;
    I = [1; 1; 1];

    B_sat

    mu = cross(N.*I.*A, B_sat)
end