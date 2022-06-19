function [Panel, BBx, BBy, BBz]  = perp2orbit(a_e, irr0, npoints, tspan, R, visibility, r, v, XYZ)
    r0 = r(1,:)';
    v0 = v(1,:)';
    % vector x, normal a la órbita
    Bx = (cross(r0',v0')/(norm(r0)*norm(v0)))'; 
    By = -r0/norm(r0);
    By = By/norm(By);

    % vector z is calculated using a cross product.
    Bz = cross(Bx,By);
    Bz = Bz/norm(Bz);
    Panel = zeros(size(th,2),4);
    BBx = zeros(3,npoints);
    BBy = zeros(3,npoints);
    BBz = zeros(3,npoints);
    % The loop is computed in order to solve for each timestep.
    for k = 1:npoints
        
        sun_dir = (R(k,:)/norm(R(k,:)))'; % vector dirección solar.
        
        % matriz de giro con velocidad angular w (rad/s).
        A56 = [1, 0,                0;
                0, cos(w*tspan(k)), -sin(w*tspan(k));
                0, sin(w*tspan(k)),  cos(w*tspan(k))];


        BBx(:,k) = (cross(r(k,:)',v(k,:)')/(norm(r(k,:))*norm(v(k,:))))';
        BBy(:,k) = By;
        BBz(:,k) = Bz;

        A15 = [Bx,By,Bz]';

        sun_body = A56*A15*sun_dir;

        Dot_Zp = [0,0,1]*sun_body;
        Dot_Zn = [0,0,-1]*sun_body;
        Dot_Yp = [0,1,0]*sun_body;
        Dot_Yn = [0,-1,0]*sun_body;

        irr = a_e^2*irr0/norm(R(k,:))^2;
        Panel(k,1) = irr*Dot_Zp*(Dot_Zp>cosd(75))*visibility(k);
        Panel(k,2) =  irr*Dot_Zn*(Dot_Zn>cosd(75))*visibility(k);
        Panel(k,3) = irr*Dot_Yp*(Dot_Yp>cosd(75))*visibility(k);
        Panel(k,4) =  irr*Dot_Yn*(Dot_Yn>cosd(75))*visibility(k);
    end