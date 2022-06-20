function sun_orbit(r_e, R, npoints, visibility)
    figure
    orbit_earth = plot3(r_e(:,1),r_e(:,2),r_e(:,3),'Color','r','LineWidth',2);
    hold on
    orbit_spacecraft = plot3(R(:,1),R(:,2),R(:,3),'Color','b','LineWidth',2);

    % sun = plot3(0,0,0,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    earth = plot3(nan,nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    sc = plot3(nan,nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    % link = plot3(nan,nan,nan,'r','LineWidth',2);

    axis equal
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    zlim([min([min(R(:,3)) -abs(max(R(:,1))-min(R(:,1)))/8 -abs(max(R(:,2))-min(R(:,2)))/8]), max([max(R(:,3)) abs(max(R(:,1))-min(R(:,1)))/8 abs(max(R(:,2))-min(R(:,2)))/8])])
    %%
    for k = 1:round(npoints/1000):npoints
        set(earth,'XDATA',r_e(k,1))
        set(earth,'YDATA',r_e(k,2))
        set(earth,'ZDATA',r_e(k,3))
        
        set(sc,'XDATA',R(k,1))
        set(sc,'YDATA',R(k,2))
        set(sc,'ZDATA',R(k,3))
        
    %     if visibility(k)
    %         set(link,'XDATA',[0 R(k,1)])
    %         set(link,'YDATA',[0 R(k,2)])
    %         set(link,'ZDATA',[0 R(k,3)])
    %     else
    %         set(link,'XDATA',nan)
    %         set(link,'YDATA',nan)
    %         set(link,'ZDATA',nan)
    %     end
        
        drawnow
    end