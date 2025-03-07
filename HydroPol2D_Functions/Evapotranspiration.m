%% Process Potential Evapotranspiration
% Author: Marcus Nobrega, Ph.D., & Matheus Schroeder dos Santos M.Sc.
% Date: 02/24/2024
% Objective: Develop a 2-D array containing the potential evapotranspiration (ETP)
% for all cells in a given watershed.
% Reference: Marco Antônio Fonseca Conceição. (2006). "Calculation Guide for Reference 
% Evapotranspiration Using the Penman-Monteith-FAO Method." EMBRAPA Technical Circular.

function [ETP,Ep] = Evapotranspiration(DEM,Temp,Temp_max,Temp_min,Day_year,lat,U_2_input,U_R_input,Krs,albedo,G_input)

%% Input Data

theta_S_B_input = 4.903*(10^-9); % Stefan-Boltzmann constant (MJ m^-2/day)

% Get the dimensions of the Digital Elevation Model (DEM)
n_rows = size(DEM,1);
n_cols = size(DEM,2);

% Assign input values to variables
latitude = lat;
G = G_input;
T = Temp;
U_2 = U_2_input;
U_R = U_R_input;
Tmax = Temp_max;
Tmin = Temp_min;

%% Process Data
% Identify invalid DEM values
neg_DEM = DEM < -200;   % Threshold for invalid negative elevations
high_DEM = DEM > 99999; % Threshold for unrealistically high elevations
inf_nan_MAPS = isinf(DEM) + isnan(DEM) + neg_DEM + high_DEM; % Logical array for invalid values
idx = inf_nan_MAPS > 0; % Index of invalid cells

% Apply mask to all input variables where DEM is invalid
G(idx) = nan;
T(idx) = nan;
U_2(idx) = nan;
Tmax(idx) = nan;
Tmin(idx) = nan;

% Assign day of the year
J = Day_year * ones(n_rows, n_cols);
J(idx) = nan;

% Apply mask to albedo
albedo(idx) = nan;

% Apply mask to Stefan-Boltzmann constant
theta_S_B = theta_S_B_input * ones(n_rows, n_cols);
theta_S_B(idx) = nan;

%% Compute Evapotranspiration
phi = latitude * pi / 180; % Convert latitude to radians
dec_sol = 0.409 * sin((2 * pi / 365) * J - 1.39); % Solar declination (radians)

% Compute the sunrise hour angle (ws)
X = (1 - ((tan(phi)).^2).*((tan(dec_sol)).^2));
X(X<= 0) = 0.00001;

ws = (pi / 2) - atan((-tan(phi) .* tan(dec_sol)) ./ (X.^0.5)); % Hour angle at sunrise (radians)

% Compute relative Earth-Sun distance
dr = 1 + 0.033 * cos((2 * pi / 365) * J);

% Compute solar radiation at the top of the atmosphere (MJ m^-2/day)
Ra = (118.08 / pi) .* dr .* (ws .* sin(phi) .* sin(dec_sol) + cos(phi) .* cos(dec_sol) .* sin(ws));

% Compute incident solar radiation (MJ m^-2/day) based on temperature variation
Rs_input = Krs .* Ra .* sqrt((Tmax - Tmin));

% Compute incident solar radiation under clear skies (MJ m^-2/day)
Rso = (0.75 + 2 * (10^-5) .* DEM) .* Ra;

% Compute net shortwave radiation (MJ m^-2/day)
Rns = (1 - albedo) .* Rs_input;

% Compute vapor pressure (kPa)
e_s = 0.6108 * exp((17.27 * T) ./ (T + 237.3)); % Saturation vapor pressure
e_a = 0.61 * exp((17.27 * Tmin) ./ (Tmin + 237.3)); % Actual vapor pressure (alternative if U_R is unavailable)

% Compute net longwave radiation (MJ m^-2/day)
Rnl = theta_S_B .* ((((Tmax + 273.16).^4) + ((Tmin + 273.16).^4)) / 2) .* (0.34 - 0.14 * sqrt(e_a)) .* ...
      (1.35 * (Rs_input ./ Rso) - 0.35);

% Compute daily net radiation (MJ m^-2/day)
Rn = Rns - Rnl;

% Compute the slope of the vapor pressure curve (kPa/°C)
delta = (4098 * (0.6108 * exp((17.27 * T) ./ (T + 237.3)))) ./ ((T + 237.3).^2);

% Compute local atmospheric pressure (kPa)
Patm = 101.3 * (((293 - 0.0065 * DEM) / 293).^5.26);

% Compute psychrometric constant (kPa/°C)
gama = 0.665 * (10^-3) * Patm;

% Compute reference evapotranspiration (ETP) using the Penman-Monteith equation (mm/day)
ETP = (0.408 * delta .* (Rn - G) + (gama * 900 .* U_2 .* (e_s - e_a)) ./ (T + 273)) ./ ...
      (delta + gama .* (1 + 0.34 .* U_2));

Ep = (0.408 * delta .* (Rn - G) + (gama * 900 .* U_2 .* (e_s - e_a)) ./ (T + 273)) ./ ...
      (delta + gama .* (1 + 0)); % No surface resistance

% Apply mask to ETP for invalid data points
ETP(idx) = nan;
Ep(idx) = nan;

end