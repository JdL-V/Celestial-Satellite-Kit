function [T_mx, T_my, T_mz] = mag_3(dt, tspan, J, BBx, BBy, BBz, B, w, N, A, I)
    npoints = length(tspan);
    T_m = zeros(npoints, 3);
    w = 0;
    for k = 2:npoints
        A56 = [1, 0,                0;
               0, cos(w*tspan(k)), -sin(w*tspan(k));
               0, sin(w*tspan(k)),  cos(w*tspan(k))];
            
        A15 = [BBx(:,k), BBy(:,k)/norm(BBy(:,k)), BBz(:,k)/norm(BBz(:,k))]';
        B_sat = A56*A15*B(:,k);
        
        T_mx(k,:) = cross(N.*I.*A.*[1 0 0], B_sat);
        T_my(k,:) = cross(N.*I.*A.*[0 1 0], B_sat);
        T_mz(k,:) = cross(N.*I.*A.*[0 0 1], B_sat);
    end
end