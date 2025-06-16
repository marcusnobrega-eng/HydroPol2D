function [ax] = surfmap(M)
% =========================================================================
% 📊 SURFMAP: Nicely formatted surface plot of a 2D matrix
% -------------------------------------------------------------------------
% Author.....: Maria Castro (adapted)
% Purpose....: Generate a top-down view of a matrix using smooth color
%              gradients and scientific formatting for visualization.
% =========================================================================

    M = double(M);  % Ensure input is double

    % === 🛠️ Flip matrix to correct orientation (row 1 = top) ===
    M = flipud(M);  % Flip vertically to match geospatial image convention

    % === 🎨 Create figure (50% of screen, centered) ===
    screenSize = get(0, 'ScreenSize');
    figWidth  = 0.5 * screenSize(3);
    figHeight = 0.5 * screenSize(4);
    figX = (screenSize(3) - figWidth) / 2;
    figY = (screenSize(4) - figHeight) / 2;

    figure('Position', [figX, figY, figWidth, figHeight]);

    % === 📈 Surface plot ===
    ax = surf(M);
    view(2);               % Top-down
    shading interp;        % Smooth gradients
    colormap('turbo');    % Scientific-friendly colormap
    colorbar;              % Add colorbar

    % === 🧭 Axis formatting ===
    axis tight;
    axis on;
    set(gca, 'YDir', 'normal');  % Ensure correct direction
    grid on;

    % === ✍️ Aesthetics ===
    set(gca, 'FontName', 'Garamond');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.2);
    title('Surface Map', 'FontSize', 16);

end
