function [Panel, BBx, BBy, BBz] = norm2sun(a_e, irr0, npoints, tspan, R, visibility, r, v, XYZ, th, w)
    npanels = 3;
    Panel = zeros(size(th,2),npanels);
    www = [1; 0; 0];
    for k = 1:npoints
        Bx = (R(k,:)/norm(R(k,:)));
        By = cross(Bx, www);
        Bz = cross(Bx, By);
        BBx(:,k) = Bx';
        BBy(:,k) = By'/norm(By);
        BBz(:,k) = Bz'/norm(Bz);

        irr = a_e^2*irr0/norm(R(k,:))^2;
        Panel(k,:) = ones(1,npanels)*irr*visibility(k);
    end