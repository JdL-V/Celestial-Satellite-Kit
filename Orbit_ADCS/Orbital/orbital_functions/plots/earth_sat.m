function earth_sat(a, r, tspan, npoints, alpha)
    om_e = (2*pi)/([23 56 4.09053]*[3600 60 1]');
    h = figure();
    ax = axes('XLim',[-2*a 2*a], 'YLim',[-2*a 2*a], 'ZLim',[-2*a 2*a]);
    orb = plot3(r(:,1), r(:,2), r(:,3));
    grid on
    hold on
    axis equal
    % axis ([-a a -a a -a a ])
    globe = Celestial_body('earth',1);
    focus = globe;
    t1 = hgtransform('Parent', ax);
    t2 = hgtransform('Parent', ax);
    t3 = hgtransform('Parent', ax);
    set(globe, 'Parent', t1);
    set(focus,'Parent', t2);
    set(orb,'Parent', t3);
    sat = plot3(nan,nan,nan,'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',7);
    globe.FaceAlpha = alpha;

    r_max=norm(max(r));
    % dummyPlot = plot3(r_max*[1/2,-1,-1,1/2,1/2,-1,-1,1/2], r_max*[1,1,-1,-1,1,1,-1,-1], r_max*[-1,-1,-1,-1,1,1,1,1], 'LineStyle', 'none');
    xlabel('x')
    ylabel('y')
    zlabel('z')

    rot_earth = tspan.*om_e;
    R_equator = makehgtform('xrotate', deg2rad(-23.4));
    for j=1:round(npoints/350):npoints
        % ax.CameraTarget=[r(j,1) r(j,2) r(j,3)];
        Txyz = makehgtform('translate',[r(j,1) r(j,2) r(j,3)]);
        Rz = makehgtform('zrotate',rot_earth(j));
    %     Ry= makehgtform('scale',9*j);
        set(t1,'Matrix',Txyz*R_equator*Rz)
        set(t2,'Matrix',Rz)
    %     set(t3,'Matrix',Ry)
        set(sat,'XData',r(j,1),'YData',r(j,2),'ZData',r(j,3));
        drawnow
    end

    % figure
    % orbit_spacecraft = plot3(r(:,1),r(:,2),r(:,3),'Color','b','LineWidth',2);
    % hold on
    % sc = plot3(nan,nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','k');

    % axis equal

    % %%
    % for k = 1:round(npoints/1000):npoints
    %     set(sc,'XDATA',r(k,1))
    %     set(sc,'YDATA',r(k,2))
    %     set(sc,'ZDATA',r(k,3))
    %     drawnow
    % end