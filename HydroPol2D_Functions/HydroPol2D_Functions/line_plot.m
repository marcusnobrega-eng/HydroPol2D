function line_plot(x_value,x_symbol,x_units,y_values,symbol,units,topBC,TimeBC,symbol_T,units_T,title_name,n_plots,id_plot)
% Create a hydrograph in left axis with a boundary condition in the reverse
% right axis.
% Adapted for the Semi-Distributed Model
% x_value: time vector 
% x_symbol: symbol of x e.g., 'x'
% y_values: discharge or other vector
% symbol: symbol of the y_values (e.g., 'Q')
% units: units of the y_values (e.g., '\mathrm{mm \cdot h^{-1}}
% topBC: rainfall or other top B.C vector to plot with the y_values
% symbol_T: symbol of the top B.C. (e.g., 'i')
% units_T: units of the top B.C.
% title_name: Name to plot at the top of the graph
% Number of series in the plot
% id_plot: index of the plot to choose the color

c = linspecer(n_plots);

matrix = topBC;
topBC = 0;
% Adjusting TopBC
if size(matrix,3) > 1 || isempty(matrix)
    for i = 1:size(matrix,3)
        data = matrix(:,:,i);
        topBC(1,i) = mean(data(:));
    end
else % Single value per area
    topBC = matrix;
end
if length(topBC) ~=1
    yyaxis left
end
set(gca,'YColor','black');
plot(x_value,y_values,'Color',c(id_plot,:),'linewidth',2)
xlabel(['$ ' x_symbol '~[ \mathrm{ ' x_units '}]$'],'Interpreter','latex')
ylabel(['$ ' symbol '~[ \mathrm{ ' units '}]$'],'Interpreter','latex');
set(gca,'Fontsize',14)
set(gca,'FontName','garamond');
set(gca,'TickDir','out');
hold on
try
    ylim([min(min((y_values))),1.5*max(max((y_values)))])
catch
    % do nothing
end

if length(topBC) ~= 1
    yyaxis right
    set(gca,'ydir','reverse');
    set(gca,'YColor','black')    
    plot(TimeBC(1:length(topBC)),topBC,'Color',[0 128 128]/256,'linewidth',1.5)
    ylabel(['$ ' symbol_T '~[ \mathrm{ ' units_T '}]$'],'Interpreter','latex');
    try
        ylim([0,10*max(max(topBC))])
    catch
        ylim([0,1])
    end
else

end

title(title_name,'Interpreter','latex','FontSize',14);

end