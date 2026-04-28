function [t, i, idf] = Huff_Curves(td, dt, K, a, b, c, RP, flag_plot, path)
% Developer: Marcus Nobrega, Ph.D
% Updated to save outputs in the user-defined "path"
%
% Goal:
%   Create a Huff hyetograph from an IDF curve
%
% Inputs:
%   td        : Rainfall duration (min)
%   dt        : Time-step (min)
%   K,a,b,c   : IDF parameters such that
%               i(mm/hr) = K*(RP^a)/((b+td)^c)
%   RP        : Return period (years)
%   flag_plot : Boolean [0,1]
%               1 = plot and export results
%               0 = do not plot
%   path      : Folder path where output files will be saved
%
% Outputs:
%   t   : Vector with time from dt to td, with dt intervals (min)
%   i   : Selected Huff rainfall intensity for each t (mm/h)
%   idf : IDF rainfall intensity for each t (mm/h)
%
% Example:
%   [t, i, idf] = Huff_Curves(60, 5, 819, 0.138, 10, 0.75, 10, 1, 'C:\Results')
%
% Notes:
%   - The following files are exported to "path" when flag_plot = 1:
%       Quartiles_Huff.pdf
%       Intensity_Huff_All_Quartiles.pdf
%       Intensity_Huff.pdf
%   - If the folder in "path" does not exist, it will be created.
%   - The function linspecer must be available in your MATLAB path if used.

%% 0. Check input path
if nargin < 9 || isempty(path)
    path = pwd; % Use current folder if no path is provided
end

if ~exist(path, 'dir')
    mkdir(path);
    disp(['Folder created successfully: ', path]);
end

%% 1. Set duration
t = dt:dt:td;
nsteps = length(t);

%% 2. Calculate IDF intensity
idf = K.*(RP^a)./((b + t).^c);      % Rainfall intensity for all t durations (mm/h)
rainfal_volume = idf(end)*td/60;    % Total rainfall depth (mm)

delta_t = zeros(1, length(idf));
for k = 1:length(idf)
    if k == 1
        delta_t(k) = t(k);
    else
        delta_t(k) = t(k) - t(k-1);
    end
end

p = idf .* delta_t / 60; %#ok<NASGU> % Incremental rainfall depth (mm), kept for consistency

%% 3. Huff Curves Method
% Huff Curves represented as a 6th-order polynomial for each quartile:
% p/pt = c6(t/td)^6 + c5(t/td)^5 + c4(t/td)^4 + c3(t/td)^3 + c2(t/td)^2 +
%        c1(t/td)^1 + c0
%
% Rows correspond to the 1st, 2nd, 3rd, and 4th quartiles
huff_coefficients = [...
    -0.9633   3.8869   -7.8950   10.0890   -8.0108   3.8936   -0.0032;
   -39.4360 125.1800 -146.0400   73.6040  -13.9360   1.6243   -0.0068;
    46.5420 -131.5500 132.6300  -57.3150   10.7960  -0.1107    0.0050;
   -25.2890  67.5400  -64.9260   28.0310   -5.2061   0.8535   -0.0042];

t_prime = t / td;
p_prime_quartiles = zeros(4, nsteps);

for k = 1:nsteps
    for j = 1:size(huff_coefficients, 2)
        expn = size(huff_coefficients, 2) - j;
        p_prime_quartiles(:, k) = p_prime_quartiles(:, k) + huff_coefficients(:, j) * (t_prime(k)^expn);
    end
end

% Boundary limits for cumulative fraction
p_prime_quartiles(p_prime_quartiles < 0) = 0;
p_prime_quartiles(p_prime_quartiles > 1) = 1;

%% 4. Calculate rainfall intensities for each quartile
rainfall_intensity = zeros(size(p_prime_quartiles));

for k = 1:nsteps
    if k == 1
        rainfall_intensity(:, k) = rainfal_volume * p_prime_quartiles(:, k) / (dt/60); % mm/h
    else
        rainfall_intensity(:, k) = rainfal_volume * ...
            (p_prime_quartiles(:, k) - p_prime_quartiles(:, k-1)) / (dt/60); % mm/h
    end
end

rainfall_intensity(rainfall_intensity < 0) = 0; % Boundary condition

%% 5. Select the recommended Huff quartile rainfall
if td <= 120
    rain = rainfall_intensity(1, :);
elseif td <= 12*60
    rain = rainfall_intensity(2, :);
elseif td <= 24*60
    rain = rainfall_intensity(3, :);
else
    rain = rainfall_intensity(4, :);
end

%% 6. Outputs
i = rain(:);
t = t(:);
idf = idf(:);

%% 7. Plotting and exporting
if flag_plot == 1
    close all
    colors = linspecer(4);

    % ---------------------------
    % Figure 1: Huff cumulative rainfall quartiles
    % ---------------------------
    figure
    set(gcf, 'Units', 'inches', 'Position', [3,3,6.5,4])

    plot(t_prime, p_prime_quartiles(1,:), 'Color', colors(1,:), 'LineStyle', '-',  'LineWidth', 2, 'Marker', '*'); hold on
    plot(t_prime, p_prime_quartiles(2,:), 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2, 'Marker', '.');
    plot(t_prime, p_prime_quartiles(3,:), 'Color', colors(3,:), 'LineStyle', '-.', 'LineWidth', 2, 'Marker', '<');
    plot(t_prime, p_prime_quartiles(4,:), 'Color', colors(4,:), 'LineStyle', ':',  'LineWidth', 2, 'Marker', '>');
    hold off

    xlabel('$t/t_d$', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$P/P_t$', 'Interpreter', 'latex', 'FontSize', 12)
    title('Huff Cumulative Rainfall', 'Interpreter', 'latex', 'FontSize', 12)
    legend('Quartile 1', 'Quartile 2', 'Quartile 3', 'Quartile 4', ...
        'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12)

    set(gca, 'TickLength', [0.015 0.01]);
    set(gca, 'TickDir', 'out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    grid on
    box on

    exportgraphics(gcf, fullfile(path, 'Quartiles_Huff.pdf'), 'ContentType', 'vector')

    % ---------------------------
    % Figure 2: All quartile hyetographs
    % ---------------------------
    figure
    set(gcf, 'Units', 'inches', 'Position', [3,3,6.5,4])

    hold on
    plot(t, rainfall_intensity(1,:), 'Color', colors(1,:), 'LineStyle', '-',  'LineWidth', 2, 'Marker', '*');
    plot(t, rainfall_intensity(2,:), 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2, 'Marker', '.');
    plot(t, rainfall_intensity(3,:), 'Color', colors(3,:), 'LineStyle', '-.', 'LineWidth', 2, 'Marker', '<');
    plot(t, rainfall_intensity(4,:), 'Color', colors(4,:), 'LineStyle', ':',  'LineWidth', 2, 'Marker', '>');
    hold off

    xlabel('Elapsed Time [min]', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('Rainfall Intensity [$\mathrm{mm.h^{-1}}$]', 'Interpreter', 'latex', 'FontSize', 12)
    title(['RP = ', num2str(RP, '%.2f'), ' years'], 'Interpreter', 'latex')
    legend('Quartile 1', 'Quartile 2', 'Quartile 3', 'Quartile 4', ...
        'Location', 'best', 'Interpreter', 'latex')

    grid on
    axis tight
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca, 'TickDir', 'out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    box on

    exportgraphics(gcf, fullfile(path, 'Intensity_Huff_All_Quartiles.pdf'), 'ContentType', 'vector')

    % ---------------------------
    % Figure 3: Selected Huff hyetograph
    % ---------------------------
    figure
    set(gcf, 'Units', 'inches', 'Position', [3,3,6.5,4])

    bar(t, rain, 'FaceColor', [0 .5 .5], 'EdgeColor', [0 .55 .55], 'LineWidth', 1.5);

    xlabel('Elapsed Time [min]', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('Rainfall Intensity [$\mathrm{mm.h^{-1}}$]', 'Interpreter', 'latex', 'FontSize', 12)
    title(['RP = ', num2str(RP, '%.2f'), ' years'], 'Interpreter', 'latex')
    legend('Huff Hyetograph', 'Location', 'best', 'Interpreter', 'latex')

    grid on
    axis tight
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca, 'TickDir', 'out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    box on

    exportgraphics(gcf, fullfile(path, 'Intensity_Huff.pdf'), 'ContentType', 'vector')

    disp(['Figures exported successfully to: ', path]);
end

end