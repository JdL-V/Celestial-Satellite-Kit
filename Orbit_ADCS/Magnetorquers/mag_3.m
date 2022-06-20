function [mu] = mag_3(dt, tspan, J, BBx, BBy, BBz, B, w, N, A, I)
    npoints = length(tspan);
    mu = zeros(npoints, 3);
    for k = 2:npoints
        A56 = [1, 0,                0;
               0, cos(w*tspan(k)), -sin(w*tspan(k));
               0, sin(w*tspan(k)),  cos(w*tspan(k))];
            
        A15 = [BBx(:,k), BBy(:,k)/norm(BBy(:,k)), BBz(:,k)/norm(BBz(:,k))]';
        B_sat = A56*A15*B;
        
        mu(k,:) = cross(N.*I.*A, B_sat);
    end
end