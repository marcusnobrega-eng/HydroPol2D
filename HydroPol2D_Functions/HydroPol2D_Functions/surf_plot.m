%% Surf Plots
function [handle] = surf_plot(H,time_plot,symbol,units,matrix,shad,plot_title,color_disc,face_alpha,new_figure,angle,X,Y)
% Developer: Marcus Nobrega, Ph.D

if new_figure == 1
    figure
end
% Video
time_frame = time_plot(1); % days
vector_data = matrix;
handle = surf(X,Y,vector_data,'FaceAlpha',face_alpha); view(angle);
if plot_title == 1
    title(['t = ',num2str(time_frame),' [days]'],'interpreter','latex')
end
colormap(turbo(color_disc));
axis tight
if shad == 1
    shading interp
end
clb = colorbar;
if H == 0
    H = 1;
end

if min(min(vector_data)) > 0
    min_data = 0;
else
    min_data = min(min(vector_data));
end

try
    caxis([min_data, H]);
catch ME
    if min(min(vector_data)) > 0
        caxis([0.5*min_data, 2*H]);
    end
end
ylabel(clb,['$ ' symbol '~[ \mathrm{ ' units '}]$'],'Interpreter','latex','FontSize',14);
clb.TickDirection = 'out';
try
    zlim([min(min(vector_data)), H])
catch
    zlim([0.5*min(min(vector_data)), 2*H])
end
set(gca,'TickDir','out');
set(gca,'FontName','Garammond')
set(gca,'FontSize',14);
set(gca,'FontWeight','Bold','LineWidth', 1.5);
handle.EdgeColor = 'none';
grid on
box on
hold on
% Marker Properties
% 'marker','v','MarkerEdgeColor','black','MarkerSize',1,'MarkerFaceColor','red');
xlabel('$x$~$[\mathrm{m}]$','Interpreter','latex');
ylabel('$y$~$[\mathrm{m}]$','Interpreter','latex');
hold off
end
