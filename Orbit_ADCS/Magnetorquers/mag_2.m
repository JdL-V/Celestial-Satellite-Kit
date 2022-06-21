function [I, mu] = mag_2(dt, N, A, J, BBx, B)
    D = -1.2591e-8;
    npoints = size(BBx, 2);
    I = zeros(1, npoints - 2);
    mu = zeros(1, npoints - 2);
%     I2 = zeros(3, npoints - 2);
    for k = 2:npoints - 1
        dth1 = (acos(dot(BBx(:,k), BBx(:,k-1))));
        dth2 = (acos(dot(BBx(:,k+1), BBx(:,k))));
        Be = norm(cross(BBx(:,k), B(:,k)));
        ddth = 2*(dth2 - dth1)/dt^2;
        I(k-1) = (J*abs(ddth) - D)/(N*A*Be);
        mu(k-1) = (J*abs(ddth) - D)/Be; 
        % A2 = 1e-6*[30 11 11].^2
        % I2(:,k-1) = (J*ddth - D)/(N*A2.*cross(BBx(:,k), B(:,k)));
    end
    % I3 = sum(I2)
end
