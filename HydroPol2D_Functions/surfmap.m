function [ax] = surfmap(M)
    close all
    figure('units','normalized','outerposition',[0 0 1 1]);
    % Creates a map of a given matrix
    ax = surf(M);  view(0,90); shading interp; colormap jet; colorbar;
    set(gca,'Fontname','Garamond');
    set(gca,'Fontsize',14);
end