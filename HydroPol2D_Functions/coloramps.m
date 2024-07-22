%% ----------------- Coloramps --------------------- %%
function [Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete] = coloramps()
%% Spectrum
RGB=[0.1127         0    0.3515
     0.2350         0    0.6663
     0.3536         0    1.0000
     0.4255         0    1.0000
     0.4384         0    1.0000
     0.3888         0    1.0000
     0.2074         0    1.0000
          0         0    1.0000
          0    0.4124    1.0000
          0    0.6210    1.0000
          0    0.7573    0.8921
          0    0.8591    0.6681
          0    0.9642    0.4526
          0    1.0000    0.1603
          0    1.0000         0
          0    1.0000         0
          0    1.0000         0
          0    1.0000         0
     0.4673    1.0000         0
     0.8341    1.0000         0
     1.0000    0.9913         0
     1.0000    0.8680         0
     1.0000    0.7239         0
     1.0000    0.5506         0
     1.0000    0.3346         0
     1.0000         0         0
     1.0000         0         0
     1.0000         0         0
     1.0000         0         0
     0.9033         0         0
     0.7412         0         0
     0.5902         0         0];
Spectrum=interp1(linspace(1, 255, 32),RGB,[1:1:255]);

%% Water Depths
% Define the colors used in the color ramp
colors = [
    0.27 0.47 0.68;   % blue
    0.41 0.68 0.84;   % light blue
    0.95 0.95 0.72;   % light yellow
    0.98 0.74 0.43;   % orange
    0.85 0.33 0.10;   % dark orange
    1.00 0.00 0.00;   % red
    0.50 0.00 0.25;   % light purple
    0.50 0.00 0.50    % purple
];

% Define the values corresponding to each color
values = [
    0.0  % blue
    0.5  % light blue
    1.0  % light yellow
    1.5  % orange
    2.0  % dark orange
    2.5  % red
    3.0  % light purple
    4.0  % purple
];

% Generate a smoothly varying colormap
depth_ramp = interp1(values, colors, linspace(min(values), max(values), 255), 'pchip');

% Example usage:
% surf(peaks);

%% Terrain

% Define the colors
color1 = [0, 100, 0] / 255;  % Dark Green
color2 = [152, 251, 152] / 255;  % Light Green
color3 = [255, 255, 0] / 255;  % Yellow
color4 = [184, 134, 11] / 255;  % Dark Yellow
color5 = [139, 69, 19] / 255;  % Brown
color6 = [101, 67, 33] / 255;  % Dark Brown
color7 = [255, 0, 0] / 255;  % Red
color8 = [139, 0, 0] / 255;  % Dark Red

% Define the number of colors in the color ramp
numColors = 256;

% Create a matrix with the color values
colors = [linspace(color1(1), color2(1), numColors/8)', ...
          linspace(color1(2), color2(2), numColors/8)', ...
          linspace(color1(3), color2(3), numColors/8)'; ...
          linspace(color2(1), color3(1), numColors/8)', ...
          linspace(color2(2), color3(2), numColors/8)', ...
          linspace(color2(3), color3(3), numColors/8)'; ...
          linspace(color3(1), color4(1), numColors/8)', ...
          linspace(color3(2), color4(2), numColors/8)', ...
          linspace(color3(3), color4(3), numColors/8)'; ...
          linspace(color4(1), color5(1), numColors/8)', ...
          linspace(color4(2), color5(2), numColors/8)', ...
          linspace(color4(3), color5(3), numColors/8)'; ...
          linspace(color5(1), color6(1), numColors/8)', ...
          linspace(color5(2), color6(2), numColors/8)', ...
          linspace(color5(3), color6(3), numColors/8)'; ...
          linspace(color6(1), color7(1), numColors/8)', ...
          linspace(color6(2), color7(2), numColors/8)', ...
          linspace(color6(3), color7(3), numColors/8)'; ...
          linspace(color7(1), color8(1), numColors/8)', ...
          linspace(color7(2), color8(2), numColors/8)', ...
          linspace(color7(3), color8(3), numColors/8)'];

% Smooth the color ramp
terrain_ramp = smoothdata(colors, 'gaussian', 10);

% % Plot the color ramp
% figure;
% colormap(terrain_ramp);
% colorbar;

%% Blues
% Define the range of colors (shades of blue)
num_colors = 100; % Number of distinct colors
start_color = [0, 0.2, 0.6]; % Lighter blue (RGB values)
end_color = [0.8, 0.9, 1]; % Darker blue (RGB values)

% Adjust contrast by changing the exponent (gamma correction)
gamma = 2; % Higher values increase contrast, lower values decrease it

% Interpolate between the start and end colors with gamma correction
r = linspace(start_color(1), end_color(1), num_colors).^gamma;
g = linspace(start_color(2), end_color(2), num_colors).^gamma;
b = linspace(start_color(3), end_color(3), num_colors).^gamma;

% Create the colormap matrix
blue_ramp = [r' g' b'];

%% Blues 2
% Define a custom colormap with distinguishable shades of blue
custom_map = [
    0.031, 0.365, 0.639;  % Deep blue
    0.098, 0.569, 0.737;  % Darker blue
    0.306, 0.675, 0.816;  % Dark blue
    0.529, 0.765, 0.882;  % Medium blue
    0.725, 0.851, 0.941;  % Soft blue
    0.882, 0.937, 0.988;  % Light blue  
    0.973, 0.988, 1.000;  % Very light blue
];

% Number of intervals for colormap
num_intervals = size(custom_map, 1);

% Number of points in colormap
num_points = 256;

% Interpolate colormap to get 256 points
blues_2 = interp1(linspace(0, 1, num_intervals), custom_map, linspace(0, 1, num_points));

%% Colors
blue_colors(1,:) = [31, 102, 169]/255;
blue_colors(2,:) = [52, 148, 204]/255;
blue_colors(3,:) = [141, 197, 228]/255;

pallete.blue_colors = blue_colors;

red_colors(1,:) = [159, 0, 0]/255;
red_colors(2,:) = [196, 103, 102]/255;
red_colors(3,:) = [216, 165, 166]/255;

pallete.red_colors = red_colors;

green_colors(1,:) = [31, 111, 112]/255;
green_colors(2,:) = [84, 162, 161]/255;
green_colors(3,:) = [159, 200, 200]/255;

pallete.green_colors = green_colors;



end
