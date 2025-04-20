function [ax] = array_plots(array_3D)
    % Get the number of plots (the third dimension of the input array)
    n_plots = size(array_3D, 3);
    
    % Set figure to full-screen
    figure('Position', get(0, 'ScreenSize'));
        close all

    for i = 1:n_plots
        % Create subplot for each 2D slice of the 3D array
        subplot(ceil(n_plots/2), min(n_plots,2), i);
        
        array_3D = double(array_3D);
        array_3D(array_3D == 0) = nan;

        % Plot the surface
        ax = surf(array_3D(:,:,i));

        % Set view to 2D (top-down)
        view(0, 90);
        
        % Set the font to Montserrat and size to 12
        set(gca, 'FontName', 'Montserrat', 'FontSize', 12);
        
        % Set labels in bold
        xlabel('X-axis', 'FontWeight', 'bold');
        ylabel('Y-axis', 'FontWeight', 'bold');

        % Apply turbo colormap and add a colorbar
        colormap('turbo');
        colorbar;

        % Adjust axis properties
        axis tight;
        shading interp;
        axis equal

        % Set ticks outside for all axes
        set(gca, 'TickDir', 'out');
        
        % Apply ticks outside to colorbar
        cbar = colorbar;
        set(cbar, 'TickDirection', 'out');
    end
end
