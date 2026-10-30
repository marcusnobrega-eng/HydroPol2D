% VAL-ET-CANOPY-001: canopy evaporation potential; negative ET is no ET.
%
% HydroPol2D-Python evaporates canopy water at the reference ET (FAO-56),
% or at the prescribed evaporation + transpiration when ET comes from maps
% (core/engine.py). MATLAB used the open-water rate Ep, which is set only
% in station mode, so with ET maps the canopy never emptied after its first
% fill. Python also clips negative reference ET and open-water evaporation
% at 0; MATLAB only replaced NaN, so a negative Penman-Monteith value on a
% cold, humid day added water to the canopy and to ponded cells.
%   A. Map mode, full canopy (LAI 3, 0.6 mm), no rain, maps E = 1 and
%      T = 3 mm/day, one hour: the canopy loses 4/24 mm.
%   B. Station mode, ETP 5 and Ep 7 mm/day: the canopy loses 5/24 mm.
%   C. update_spatial_BC with a station day whose Penman-Monteith values are
%      negative (mock/ETP_model.m): ETP and Ep are clipped at 0.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

map_loss = canopy_loss(functions_dir, 1, struct('ETP', 0, 'Ep', 0), 1, 3);
station_loss = canopy_loss(functions_dir, 0, struct('ETP', 5, 'Ep', 7), 0, 0);
map_error = abs(map_loss - 4 / 24);
station_error = abs(station_loss - 5 / 24);

% run() executes a script from its own folder, so the test double must sit
% next to a copy of update_spatial_BC.m.
mock_dir = tempname; mkdir(mock_dir);
copyfile(fullfile(case_dir, 'mock', 'ETP_model.m'), mock_dir);
copyfile(fullfile(functions_dir, 'update_spatial_BC.m'), mock_dir);
[etp, ep] = station_day(mock_dir);
rmdir(mock_dir, 's');
negative_clipped = all(etp(2:end) == 0) && all(ep(2:end) == 0) && etp(1) == 2 && ep(1) == 3;

passed = map_error < 1e-12 && station_error < 1e-12 && negative_clipped;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-ET-CANOPY-001", "Canopy potential; negative ET clipped", ...
    map_loss, station_loss, map_error, station_error, min(etp), min(ep), status, true, ...
    'VariableNames', {'case_id','case_name','map_mode_canopy_loss_mm', ...
    'station_mode_canopy_loss_mm','map_mode_error_mm','station_mode_error_mm', ...
    'min_ETP_mm_day','min_Ep_mm_day','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-ET-CANOPY-001 failed.');
end
fprintf('VAL-ET-CANOPY-001 passed.\n');

function loss = canopy_loss(functions_dir, map_mode, H, evaporation, transpiration)
ny = 2; nx = 2; k = 2; %#ok<NASGU>
flags = struct('flag_subgrid', 0, 'flag_overbanks', 0, 'flag_ETP', 1, ...
    'flag_abstraction', 1, 'flag_input_ETP_map', map_mode); %#ok<NASGU>
BC_States = struct('Liquid_Precipitation', zeros(ny, nx), ...
    'current_evaporation_map', evaporation * ones(ny, nx), ...
    'current_transpiration_map', transpiration * ones(ny, nx));
if ~map_mode, BC_States = rmfield(BC_States, {'current_evaporation_map','current_transpiration_map'}); end %#ok<NASGU>
LAI_raster = struct('Z', 3 * ones(ny, nx)); %#ok<NASGU>
Hydro_States = struct('S', 0.6 * ones(ny, nx), 'ETP', H.ETP * ones(ny, nx), 'Ep', H.Ep * ones(ny, nx));
time_step = 60; %#ok<NASGU>          % min
zero_matrix = zeros(ny, nx); elevation = zeros(ny, nx); idx_nan = false(ny, nx); %#ok<NASGU>
DEM_raster = struct('cellsize', 10); depths = struct('d_t', zeros(ny, nx)); errors = zeros(1, 4); %#ok<NASGU>
run(fullfile(functions_dir, 'interception_module.m'));
loss = 0.6 - Hydro_States.S(1);
end

function [etp, ep] = station_day(functions_dir)
ny = 2; nx = 3; k = 1; %#ok<NASGU>
flags = struct('flag_rainfall', 0, 'flag_spatial_rainfall', 0, 'flag_input_rainfall_map', 0, ...
    'flag_satellite_rainfall', 0, 'flag_real_time_satellite_rainfall', 0, 'flag_inflow', 0, ...
    'flag_stage_hydrograph', 0, 'flag_ETP', 1, 'flag_input_ETP_map', 0, 'flag_snow_modeling', 0); %#ok<NASGU>
Elevation_Properties = struct('elevation_cell', zeros(ny, nx)); %#ok<NASGU>
idx_nan = false(ny, nx); %#ok<NASGU>
Hydro_States = struct();
BC_States = struct(); %#ok<NASGU>
ETP_Parameters = struct('climatologic_spatial_duration', [0; 1440], ...
    'time_ETP', datetime(2020, 1, 1) + days([0; 1]), 'coordinates_stations', [0 0], ...
    'maxtemp_stations', 0, 'mintemp_stations', 0, 'avgtemp_stations', 0, 'u2_stations', 0, ...
    'ur_stations', 0, 'G_stations', 0, 'DEM_etp', zeros(ny, nx), 'lat', 45, 'Krs', 0.16, ...
    'alfa_albedo_input', 0.23); %#ok<NASGU>
Spatial_Rainfall_Parameters = struct('x_grid', zeros(nx, ny), 'y_grid', zeros(nx, ny)); %#ok<NASGU>
Inflow_Parameters = struct('n_stream_gauges', 0); %#ok<NASGU>
Rainfall_Parameters = struct(); Wshed_Properties = struct(); %#ok<NASGU>
t_previous = 0; t = 10; %#ok<NASGU>
run(fullfile(functions_dir, 'update_spatial_BC.m'));
etp = Hydro_States.ETP(:); ep = Hydro_States.Ep(:);
end
