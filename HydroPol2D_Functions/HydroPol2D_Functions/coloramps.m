%% ----------------- Coloramps --------------------- %%
function [Spectrum,depth_ramp,terrain_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS,WSE_RAS] = coloramps()
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

%% Maria Pallete
temp =      [0,255,255;
            0,247,251;
            0,239,248;
            0,230,244;
            0,222,240;
            0,214,236;
            0,206,233;
            0,197,229;
            0,189,225;
            0,181,221;
            0,173,218;
            0,165,214;
            0,156,210;
            0,148,206;
            0,140,203;
            0,132,199;
            0,123,195;
            0,115,191;
            0,107,188;
            0,99,184;
            0,90,180;
            0,82,176;
            0,74,173;
            0,66,169;
            0,58,165;
            0,49,161;
            0,41,158;
            0,33,154;
            0,25,150;
            0,16,146;
            0,8,143;
            0,0,139]/255;

Depth_RAS =interp1(linspace(1, 255, size(temp,1)),temp,[1:1:255]);
% colormap(Depth_RAS); colorbar


%% Terrain RAS
temp =      [102,205,170;
            136,216,132;
            171,228,93;
            206,239,55;
            241,250,16;
            222,239,0;
            165,210,0;
            107,181,0;
            49,153,0;
            8,129,0;
            66,138,0;
            123,146,0;
            181,154,0;
            239,163,0;
            237,139,0;
            210,101,0;
            184,64,0;
            157,26,0;
            141,3,3;
            147,12,12;
            152,22,22;
            158,31,31;
            164,41,41;
            158,59,59;
            150,78,78;
            141,98,98;
            133,117,117;
            140,140,140;
            169,167,167;
            198,195,195;
            227,223,223;
            255,250,250]/255;


Terrain_RAS =interp1(linspace(1, 255, size(temp,1)),temp,[1:1:255]);
% colormap(Terrain_RAS); colorbar

%% Velocity RAS

temp =  [0,0,139;
        0,0,167;
        0,0,195;
        0,0,223;
        0,0,252;
        30,50,232;
        65,107,205;
        100,165,178;
        135,223,151;
        163,241,119;
        191,245,84;
        217,249,49;
        244,253,14;
        255,242,0;
        255,220,0;
        255,198,0;
        255,177,0;
        255,154,0;
        255,131,0;
        255,108,0;
        255,84,0;
        251,67,0;
        240,60,0;
        229,53,0;
        218,47,0;
        206,40,0;
        195,33,0;
        184,27,0;
        173,20,0;
        162,13,0;
        150,7,0;
        139,0,0]/255;

Velocity_RAS =interp1(linspace(1, 255, size(temp,1)),temp,[1:1:255]);
% colormap(Velocity_RAS); colorbar

% WSE
temp = [127	255	0
        103	231	0
        78	206	0
        53	181	0
        28	156	0
        4	132	0
        41	148	0
        91	173	0
        140	198	0
        189	222	0
        239	247	0
        255	243	0
        255	226	0
        255	209	0
        255	191	0
        255	174	0
        255	149	0
        255	117	0
        255	85	0
        255	53	0
        255	21	0
        255	12	13
        255	50	53
        255	87	92
        255	124	131
        255	161	171
        255	186	205
        255	149	215
        255	111	225
        255	74	235
        255	37	245
        255	0	255
        ]/255;
WSE_RAS = interp1(linspace(1, 255, size(temp,1)),temp,[1:1:255]);
end