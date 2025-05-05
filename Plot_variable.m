function Plot_variable(u_cc,v_cc,p_cc,X_cc,Y_cc,Nx,Ny,Delta_x,Delta_y,color,k2,t,customColormap,Re,xmin,xmax,ymin,ymax,Cylinder_display)

    % Plot
    u_cc_plot = u_cc;
    v_cc_plot = v_cc;
    p_cc_plot = p_cc;
    u_cc_plot(color.cylinder) = nan;
    v_cc_plot(color.cylinder) = nan;
    p_cc_plot(color.cylinder) = nan;
    
    w = calculate_vorticity(u_cc,v_cc,Nx,Ny,Delta_x,Delta_y,color);
    w_plot = w;
    w_plot(color.cylinder) = nan;

    xmin_plot = -5;
    xmax_plot = 10;
    ymin_plot = -5;
    ymax_plot = 5;
    
    figure(2);
    set(gcf,'Position',[100 100 1000 800])
    subplot(3,2,k2)
    contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc_plot(2:end-1,2:end-1),50,'LineStyle','none')
    colorbar
    colormap(customColormap)
    xlabel('X (m)')
    ylabel('Y (m)')
    % zlabel('u')
    axis equal
    title(['u: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
    xlim([xmin_plot,xmax_plot])
    ylim([ymin_plot,ymax_plot])
    % clim([-2,2])
    drawnow
    
    figure(3)
    set(gcf,'Position',[100 100 1000 800])
    subplot(3,2,k2)
    contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),v_cc_plot(2:end-1,2:end-1),50,'LineStyle','none')
    colorbar
    colormap(customColormap)
    xlabel('X (m)')
    ylabel('Y (m)')
    % zlabel('u')
    axis equal
    title(['v: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
    xlim([xmin_plot,xmax_plot])
    ylim([ymin_plot,ymax_plot])
    drawnow
    
    figure(4)
    set(gcf,'Position',[100 100 1000 800])
    subplot(3,2,k2)
    contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),w_plot(2:end-1,2:end-1),100,'LineStyle','none')
    colorbar
    colormap(customColormap)
    xlabel('X (m)')
    ylabel('Y (m)')
    % zlabel('u')
    axis equal
    title(['Vorticity: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
    clim([-1,1])
    xlim([xmin_plot,xmax_plot])
    ylim([ymin_plot,ymax_plot])
    drawnow

    figure(5)
    set(gcf,'Position',[100 100 1000 800])
    subplot(3,2,k2)
    contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),p_cc_plot(2:end-1,2:end-1),100,'LineStyle','none')
    colorbar
    colormap(customColormap)
    xlabel('X (m)')
    ylabel('Y (m)')
    % zlabel('u')
    axis equal
    title(['Pressure: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
    xlim([xmin_plot,xmax_plot])
    ylim([ymin_plot,ymax_plot])
    drawnow
    
    figure(6)
    set(gcf,'Position',[100 100 1000 800])
    subplot(3,2,k2)
    l = streamslice(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc(2:end-1,2:end-1),v_cc(2:end-1,2:end-1),10); hold on
    surface(X_cc, Y_cc, zeros(size(Cylinder_display)), Cylinder_display, ...
        'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'FaceAlpha', 1, 'CDataMapping', 'direct');
    % colormap([0 0 0])  % Set colormap to black
    % caxis([0 1])       % Set color range so 1 maps to black, NaN to transparent
    hold off
    axis equal
    set(l,'LineWidth',1)
    set(l, 'Color', "k")
    axis equal
    title(['Streamline: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
    xlim([xmin_plot,xmax_plot])
    ylim([ymin_plot,ymax_plot])
    drawnow

end