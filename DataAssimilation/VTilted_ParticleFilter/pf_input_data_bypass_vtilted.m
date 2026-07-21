%PF_INPUT_DATA_BYPASS_VTILTED Inputs for the full HydroPol2D PF prototype.
% This script is run only by the isolated DA workflow. It uses the normal
% bypass interface and does not change shared model inputs or source rasters.

if ~exist('model_root', 'var') || ~isfolder(model_root)
    error('The DA runner must define model_root before loading PF inputs.');
end

run(fullfile(model_root, 'Config', 'input_data_bypass_script.m'));

duration_min = str2double(getenv('HYDROPOL2D_DA_TOTAL_DURATION_MIN'));
rain_mm_h = str2double(getenv('HYDROPOL2D_DA_RAIN_MM_H'));
rain_duration_min = str2double(getenv('HYDROPOL2D_DA_RAIN_DURATION_MIN'));
if ~isfinite(duration_min) || duration_min <= 0, duration_min = 30; end
if ~isfinite(rain_mm_h) || rain_mm_h < 0, rain_mm_h = 100; end
if ~isfinite(rain_duration_min) || rain_duration_min < 0, rain_duration_min = duration_min; end
rain_duration_min = min(rain_duration_min, duration_min);

InputData_Bypass.general.date_begin = datetime(2025, 5, 1, 0, 0, 0);
InputData_Bypass.general.date_end = InputData_Bypass.general.date_begin + minutes(duration_min);
InputData_Bypass.general.routing_time = duration_min;
InputData_Bypass.general.record_time_maps = max(5, min(15, duration_min));
InputData_Bypass.general.record_time_hydrographs = max(5, min(15, duration_min));
InputData_Bypass.general.slope_outlet = 0.02;
InputData_Bypass.general.n_outlets_data = 1;

InputData_Bypass.flags.flag_rainfall = double(rain_mm_h > 0 && rain_duration_min > 0);
InputData_Bypass.flags.flag_spatial_rainfall = 0;
InputData_Bypass.flags.flag_input_rainfall_map = 0;
InputData_Bypass.flags.flag_satellite_rainfall = 0;
InputData_Bypass.flags.flag_inflow = 0;
InputData_Bypass.flags.flag_stage_hydrograph = 0;
InputData_Bypass.flags.flag_ETP = 0;
InputData_Bypass.flags.flag_input_ETP_map = 0;
InputData_Bypass.flags.flag_infiltration = 1;
InputData_Bypass.flags.flag_groundwater_modeling = 1;
InputData_Bypass.flags.flag_waterquality = 0;
InputData_Bypass.flags.flag_reservoir = 0;
InputData_Bypass.flags.flag_boundary = 0;
InputData_Bypass.flags.flag_subgrid = 0;
InputData_Bypass.flags.flag_overbanks = 0;
InputData_Bypass.flags.flag_export_maps = 0;
InputData_Bypass.flags.flag_dashboard = 0;

% Classes 1 and 3 are the symmetric hillslopes. Class 2 is the channel.
LULC = table( ...
    {'Hillslope'; 'Channel'; 'Hillslope'}, [1; 2; 3], ...
    [0.035; 0.070; 0.035], zeros(3, 1), zeros(3, 1), ...
    zeros(3, 1), zeros(3, 1), zeros(3, 1), ones(3, 1), ...
    999 * ones(3, 1), [0.45; 0.25; 0.45], [0.85; 0.70; 0.85], ...
    'VariableNames', {'LC','Index','roughness','h_0_mm','d_0_mm', ...
    'C1','C2','C3','C4','index_impervious','root_depth_m','Kc'});
InputData_Bypass.LULC = struct('table', LULC);

SOIL = table( ...
    {'Sandy loam'; 'Clay loam'; 'Sandy loam'}, [1; 2; 3], ...
    [18.0; 4.0; 18.0], [1.89; 1.31; 1.89], [6.0; 1.9; 6.0], ...
    [0.43; 0.47; 0.43], [0.07; 0.10; 0.07], [0.23; 0.26; 0.23], ...
    [0.15; 0.12; 0.15], [0.8; 0.8; 0.8], 1.5 * ones(3, 1), ...
    ones(3, 1), ones(3, 1), ones(3, 1), 0.10 * ones(3, 1), ...
    0.25 * ones(3, 1), 0.50 * ones(3, 1), ...
    'VariableNames', {'Soil_type','Index','ksat_mm_h','n_vg','alpha_vg_1_m', ...
    'theta_sat','theta_r','theta_i','Sy','ksat_gw_mm_h','Soil_Depth_m', ...
    'Ks_multiplier_near_surface','Ks_multiplier_root_zone', ...
    'Ks_multiplier_transmission','Ltop_m','dh_max_m','l_vg'});
InputData_Bypass.SOIL = struct('table', SOIL);

if rain_duration_min >= duration_min
    rainfall_time = [0; duration_min];
    rainfall_intensity = [rain_mm_h; rain_mm_h];
else
    rainfall_time = [0; rain_duration_min; duration_min];
    rainfall_intensity = [rain_mm_h; rain_mm_h; 0];
end
InputData_Bypass.Rainfall_Parameters = struct( ...
    'time_rainfall', rainfall_time, ...
    'intensity_rainfall', rainfall_intensity, ...
    'time_step_rainfall', max(min(diff(rainfall_time)), 1), ...
    'rainfall_duration', rain_duration_min, ...
    'n_obs_rainfall', numel(rainfall_time));
