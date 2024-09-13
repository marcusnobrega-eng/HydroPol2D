% ETP Model - Penman Monteith Model
% Developer: Matheus Santos & Marcus Nobrega
% Goal - Solve PM model spatially
%
% Input:
%
% Output:
% ETP: Interpolated value of ETP for each cell
function [ETP] = ETP_model(k,day_of_year,x_coordinate,y_coordinate,x_grid,y_grid,maxtemp_stations,mintemp_stations,avgtemp_stations,u2_stations,ur_stations,G_stations,DEM_etp,lat,Krs,alfa_albedo_input,idx_mask)
%% Interpolation

% Checking Missing Data
n_stations = length(x_coordinate); % First try
% Maximum Temperature
var_obs = maxtemp_stations(k,1:n_stations)';
idx1 = isnan(var_obs);
% Minimum Temperature
var_obs = mintemp_stations(k,1:n_stations)';
idx2 = isnan(var_obs);
% Average Temperature
var_obs = avgtemp_stations(k,1:n_stations)';
idx3 = isnan(var_obs);
% U2
var_obs = u2_stations(k,1:n_stations)';
idx4 = isnan(var_obs);
% UR
var_obs = ur_stations(k,1:n_stations)';
idx5 = isnan(var_obs);
% G
var_obs = G_stations(k,1:n_stations)';
idx6 = isnan(var_obs);

% Array with missed values
idx = logical(idx1 + idx2 + idx3 + idx4 + idx5 + idx6);

% Number of effective stations
n_stations = sum(idx == 0);

% Deleting Data from Missed Stations
x_coordinate(idx) = [];
y_coordinate(idx) = [];
maxtemp_stations(:,idx') = [];
mintemp_stations(:,idx')  = [];
avgtemp_stations(:,idx')  = [];
u2_stations(:,idx')  = [];
ur_stations(:,idx')  = [];
G_stations(:,idx')  = [];

% Maximum Temperature
var_obs = maxtemp_stations(k,1:n_stations)';
[max_temp,~,~] = IDW_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid);
% Minimum Temperature
var_obs = mintemp_stations(k,1:n_stations)';
[min_temp,~,~] = IDW_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid);
% Average Temperature
var_obs = avgtemp_stations(k,1:n_stations)';
[avg_temp,~,~] = IDW_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid);
% U2
var_obs = u2_stations(k,1:n_stations)';
[u2,~,~] = IDW_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid);
% UR
var_obs = ur_stations(k,1:n_stations)';
[ur,~,~] = IDW_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid);
% G
var_obs = G_stations(k,1:n_stations)';
[G,~,~] = IDW_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid);

%% Run ETP code
[ETP] = Evapotranspiration(DEM_etp, avg_temp,max_temp,min_temp,day_of_year,lat,u2, ur, Krs, alfa_albedo_input, G);
ETP(idx_mask) = nan;
end