function [ax] = array_plots(array_3D)
n_plots = size(array_3D,3);
for i = 1:n_plots
    subplot(n_plots,1,i)
    ax = surf(array_3D(:,:,i));
    view(0,90); set(gca,'FontName','Garamond','Fontsize',12);
    colormap('turbo'); colorbar;
    axis tight
    shading interp
end
