%% Maps and Hydrographs
close all
set(gcf,'units','inches','position',[0,0,7,8])
size_font = 10;
ybegin = 1;
yend = ny_max;
xbegin = 1;
xend = nx_max;
if flag_waterquality == 1
    
    ax2 = subplot(1,3,1);
    min_washed = 1e-4;
    % Tot_Washed
    z = Tot_Washed*1000; 
    z(z<=min_washed)=nan;
    z(isinf(z)) = nan;
    zmax = max(max(max(z)));
    zmin = min(min(min(z)));

    F = z([ybegin:1:yend],[xbegin:1:xend]);
    surf(x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
    shading INTERP;
    %     title(sprintf('Time(min) = %d',t_title),'Interpreter','Latex','FontSize',size_font)
    view(0,90);
    colormap(ax2,linspecer)
    cbh = colorbar;
    caxis([zmin zmax]);
    ylabel(cbh,'Total Mass Washed (g)','Interpreter','Latex','FontSize',size_font)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',size_font)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',size_font)
    zlabel ('Total Mass Washed (g)','Interpreter','Latex','FontSize',size_font)
    grid off
    box on ; 

    ax5 = subplot(1,3,2)
    % Bf
    final_mass = Pol_mass_map(:,:,end)/cell_area*1000; % g/m2
    z = final_mass; 
%     z(z<=0.01)=nan;
    z(isinf(z)) = nan;
    zmax = max(max(max(z)));
    zmin = min(min(min(z)));

    F = z([ybegin:1:yend],[xbegin:1:xend],1);
    F(idx2(:,:)) = 0;
    surf(x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
    shading INTERP;

    %     title(sprintf('Time(min) = %d',t_title),'Interpreter','Latex','FontSize',size_font)
    view(0,90);
    colormap(ax5,linspecer)
    cbh = colorbar;
    caxis([zmin zmax]);
    shading INTERP;
    ylabel(cbh,'Final Mass ($g/m^2$)','Interpreter','Latex','FontSize',size_font)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',size_font)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',size_font)
    zlabel ('Final Mass ($g/m^2$)','Interpreter','Latex','FontSize',size_font)
    grid off
    box on ; 

    ax5 = subplot(1,3,3)
    % Bf
    initial_mass_pol = Pol_mass_map(:,:,1)/cell_area*1000; % g/m2
    z = initial_mass_pol; 
%     z(z<=0.01)=nan;
    z(isinf(z)) = nan;
    zmax = max(max(max(z)));
    zmin = min(min(min(z)));

    if zmax == zmin
        zmin = 0.8*zmax;
    end

    F = z([ybegin:1:yend],[xbegin:1:xend],1);
    F(idx2(:,:)) = 0;
    surf(x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
    %     title(sprintf('Time(min) = %d',t_title),'Interpreter','Latex','FontSize',size_font)
    view(0,90);
    colormap(ax5,linspecer)
    cbh = colorbar;
    caxis([zmin zmax]);
    ylabel(cbh,'Initial Mass ($g/m^2$)','Interpreter','Latex','FontSize',size_font)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',size_font)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',size_font)
    zlabel ('Initial Mass ($g/m^2$)','Interpreter','Latex','FontSize',size_font)
    grid off
    box on ;     


end