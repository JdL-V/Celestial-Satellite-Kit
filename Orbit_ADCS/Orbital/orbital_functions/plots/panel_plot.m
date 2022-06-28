function panel_plot(BBx, BBy, BBz, tspan, visibility, T, Panel)
    figure
    plot(tspan, BBx)
    xlim([0,max(tspan)])

    figure
    plot(tspan, BBy)
    xlim([0,max(tspan)])

    figure
    plot(tspan, BBz)
    xlim([0,max(tspan)])

    figure
    plot(tspan/T,visibility,'k','LineWidth',1.5)
    yticks([0 1])
    ylim([-0.25 1.25])
    yticklabels({'no','yes'})
    xl = xlabel('$n_{orbits}$', 'interpreter', 'latex');
    yl = ylabel('Light', 'interpreter', 'latex');
    set(get(gca,'ylabel'),'Rotation',0)
    set(yl, 'Units', 'Normalized', 'Position', [-0.1, 0.91, 0]);
    set(xl, 'Units', 'Normalized', 'Position', [0.85, -0.075, 0]);
    box on


    figure
    plot(tspan/T, Panel)
    xlim([0,max(tspan)/T])
    xl = xlabel('$n_{orbits}$', 'interpreter', 'latex');
    yl = ylabel('G $\left[\mathrm{\frac{W}{m^2}}\right]$', 'interpreter', 'latex');
    set(get(gca,'ylabel'),'Rotation',0)
    set(yl, 'Units', 'Normalized', 'Position', [-0.1, 0.91, 0]);
    set(xl, 'Units', 'Normalized', 'Position', [0.85, -0.075, 0]);
    legend('Panel 1','Panel 2','Panel 3','Panel 4')
    box on
