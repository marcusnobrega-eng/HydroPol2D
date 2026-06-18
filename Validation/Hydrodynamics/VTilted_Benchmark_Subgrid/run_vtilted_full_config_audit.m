clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
model_root = fullfile(repo_root, 'HydroPol2D_Model');
functions_dir = fullfile(model_root, 'HydroPol2D_Functions');
topo_path = fullfile(model_root, 'topotoolbox-master');
static_dir = fullfile(model_root, 'Validation', 'Phase1_VTilted_Catchment', 'Static');
routing_mode = lower(strtrim(string(getenv('HYDROPOL2D_VTILT_ROUTING'))));
if strlength(routing_mode) == 0
    routing_mode = "local_inertial";
    out_dir = fullfile(case_dir, 'Outputs', 'FullConfigAudit');
else
    routing_mode = canonical_routing_mode(routing_mode);
    out_dir = fullfile(case_dir, 'Outputs', 'FullConfigRouting', char(routing_mode));
end

addpath(genpath(functions_dir));
addpath(genpath(topo_path));
addpath('/Users/mngomes/Downloads', '-end');

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

Paths = make_validation_paths(fullfile(out_dir, 'HydroPol2D_Output'), true);
InputPaths = make_vtilted_input_paths(static_dir, topo_path, functions_dir);

input_data_bypass_script_path = fullfile(case_dir, 'input_data_bypass_script_vtilted_compat.m');
if exist(input_data_bypass_script_path, 'file') ~= 2
    error('Attached input-data bypass script not found: %s', input_data_bypass_script_path);
end

use_inputpaths_bypass = 1;
use_inputdata_bypass = 1;
clean_output_folder = true;
run_postprocessing = false;
enable_logging = false;
model_folder = '';
GD = [];
resultsDir = Paths.Results;
export_root_dir = Paths.Root;

HydroPol2D_preprocessing;

routing_mode = canonical_routing_mode_from_env();
out_dir = output_dir_for_mode(case_dir, routing_mode);
write_preprocessing_audit(out_dir, InputPaths, input_data_bypass_script_path, ...
    DEM_raster, Wshed_Properties, flags, running_control, Rainfall_Parameters, ...
    Input_Rainfall, LULC_Properties, outlet_index, elevation, C_a, slope_outlet);

HydroPol2D_Main_While;

routing_mode = canonical_routing_mode_from_env();
out_dir = output_dir_for_mode(case_dir, routing_mode);
write_hydrograph_outputs(out_dir, running_control, outlet_states, ...
    outlet_runoff_volume, depths, C_a, BC_States, Wshed_Properties, ...
    Rainfall_Parameters, elevation);

function mode = canonical_routing_mode_from_env()
mode = lower(strtrim(string(getenv('HYDROPOL2D_VTILT_ROUTING'))));
if strlength(mode) == 0
    mode = "local_inertial";
else
    mode = canonical_routing_mode(mode);
end
end

function mode = canonical_routing_mode(mode)
switch lower(strtrim(string(mode)))
    case {"full_momentum", "fullmomentum", "fm"}
        mode = "full_momentum";
    case {"local_inertial", "localinertial", "li"}
        mode = "local_inertial";
    case {"cellular_automata", "ca"}
        mode = "cellular_automata";
    case {"kinematic", "kinematic_wave", "kw"}
        mode = "kinematic";
    case {"diffusive", "diffusive_wave", "dw"}
        mode = "diffusive";
    otherwise
        error('Unsupported HYDROPOL2D_VTILT_ROUTING mode: %s', mode);
end
end

function out_dir = output_dir_for_mode(case_dir, routing_mode)
if routing_mode == "local_inertial" && strlength(string(getenv('HYDROPOL2D_VTILT_ROUTING'))) == 0
    out_dir = fullfile(case_dir, 'Outputs', 'FullConfigAudit');
else
    out_dir = fullfile(case_dir, 'Outputs', 'FullConfigRouting', char(routing_mode));
end
end

function InputPaths = make_vtilted_input_paths(static_dir, topo_path, functions_dir)
InputPaths = struct();
InputPaths.case_root = fileparts(static_dir);
InputPaths.topo_path = topo_path;
InputPaths.hydropol2d_tools = functions_dir;
InputPaths.DEM_path = fullfile(static_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_dir, 'SOIL.tif');
InputPaths.DTB_path = fullfile(static_dir, 'DTB.tif');
InputPaths.GW_table_path = fullfile(static_dir, 'GW_table.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_dir, 'Initial_SM.tif');
InputPaths.LAI_path = fullfile(static_dir, 'LAI.tif');
InputPaths.Albedo_path = fullfile(static_dir, 'Albedo.tif');
InputPaths.Subgrid_DEM_path = '';
InputPaths.RiverWidths_path = '';
InputPaths.RiverDepths_path = '';
InputPaths.Warmup_Depth_path = '';
InputPaths.Initial_Buildup_path = '';
InputPaths.B1_path = '';
InputPaths.B2_path = '';
InputPaths.W1_path = '';
InputPaths.W2_path = '';
InputPaths.Rainfall_Rasters_Folder = '';
InputPaths.Transpiration_Rasters_Folder = '';
InputPaths.Evaporation_Rasters_Folder = '';
InputPaths.Rainfall_Raster_Files = {};
InputPaths.Transpiration_Raster_Files = {};
InputPaths.Evaporation_Raster_Files = {};
InputPaths.Inflow_Hydrograph_CSV = '';
InputPaths.Stage_Hydrograph_CSV = '';
InputPaths.Observed_Gauges_CSV = '';
InputPaths.ETP_input_spreadsheet = '';
InputPaths.Rainfall_Timeseries_File = '/Users/mngomes/Downloads/Rainfall_Intensity_Data.xlsx';
InputPaths.Outlet_Cells_CSV = '';
end

function Paths = make_validation_paths(root_dir, clean_output)
if exist(root_dir, 'dir') && clean_output
    delete_dir_contents(root_dir);
end
Paths = struct();
Paths.Root = root_dir;
Paths.Results = fullfile(root_dir, 'Modeling_Results');
Paths.Temp = fullfile(root_dir, 'Temporary_Files');
Paths.Logs = fullfile(root_dir, 'Logs');
Paths.FigPDF = fullfile(Paths.Results, 'Figures_PDF');
Paths.FigFIG = fullfile(Paths.Results, 'Figures_FIG');
Paths.Tables = fullfile(Paths.Results, 'Tables_CSV');
Paths.RastersWD = fullfile(Paths.Results, 'Rasters_Water_Depths');
Paths.RastersWSE = fullfile(Paths.Results, 'Rasters_WSE');
Paths.RastersStatic = fullfile(Paths.Results, 'Rasters_Static');
Paths.RastersVelocity = fullfile(Paths.Results, 'Rasters_Velocity');
Paths.RastersHazard = fullfile(Paths.Results, 'Rasters_Hazard');
Paths.WQMaps = fullfile(Paths.Results, 'Rasters_WQ');
Paths.HRMaps = fullfile(Paths.Results, 'Rasters_Human_Risk');
Paths.Anim = fullfile(Paths.Results, 'GIFs_MP4');
Paths.Shapes = fullfile(Paths.Results, 'Shapefiles');
fields = fieldnames(Paths);
for i = 1:numel(fields)
    if ~exist(Paths.(fields{i}), 'dir')
        mkdir(Paths.(fields{i}));
    end
end
end

function delete_dir_contents(p)
items = dir(p);
for i = 1:numel(items)
    name = items(i).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue
    end
    target = fullfile(p, name);
    if items(i).isdir
        rmdir(target, 's');
    else
        delete(target);
    end
end
end

function write_preprocessing_audit(out_dir, InputPaths, input_script, DEM_raster, ...
    Wshed_Properties, flags, running_control, Rainfall_Parameters, Input_Rainfall, ...
    LULC_Properties, outlet_index, elevation, C_a, slope_outlet)
[row_out, col_out] = find(outlet_index);
rain_times = Rainfall_Parameters.time_rainfall(:);
rain_intensities = Rainfall_Parameters.intensity_rainfall(:);

audit = table();
audit.item = [
    "input_data_script";
    "dem_path";
    "rainfall_file";
    "dem_rows";
    "dem_cols";
    "resolution_m";
    "routing_time_min";
    "time_step_model_min";
    "min_time_step_s";
    "max_time_step_s";
    "flag_full_momentum";
    "flag_inertial";
    "flag_CA";
    "flag_kinematic";
    "flag_diffusive";
    "flag_infiltration";
    "flag_spatial_rainfall";
    "flag_input_rainfall_map";
    "flag_subgrid";
    "slope_outlet";
    "n_outlet_cells";
    "outlet_rows";
    "outlet_cols";
    "roughness_min";
    "roughness_max";
    "rain_times_min";
    "rain_intensity_mm_h";
    "domain_area_m2"
    ];
audit.value = [
    string(input_script);
    string(InputPaths.DEM_path);
    string(InputPaths.Rainfall_Timeseries_File);
    string(size(DEM_raster.Z, 1));
    string(size(DEM_raster.Z, 2));
    string(Wshed_Properties.Resolution);
    string(running_control.routing_time);
    string(running_control.time_step_model);
    string(running_control.min_time_step);
    string(running_control.max_time_step);
    string(flags.flag_full_momentum);
    string(flags.flag_inertial);
    string(flags.flag_CA);
    string(flags.flag_kinematic);
    string(flags.flag_diffusive);
    string(flags.flag_infiltration);
    string(flags.flag_spatial_rainfall);
    string(flags.flag_input_rainfall_map);
    string(flags.flag_subgrid);
    string(slope_outlet);
    string(numel(row_out));
    strjoin(string(row_out.'), ' ');
    strjoin(string(col_out.'), ' ');
    string(min(LULC_Properties.roughness(:), [], 'omitnan'));
    string(max(LULC_Properties.roughness(:), [], 'omitnan'));
    strjoin(string(rain_times.'), ' ');
    strjoin(string(rain_intensities.'), ' ');
    string(sum(C_a(~isnan(elevation)), 'all', 'omitnan'))
    ];
writetable(audit, fullfile(out_dir, 'Preprocessing_Config_Audit.csv'));

rain = table(rain_times, rain_intensities, ...
    'VariableNames', {'Time_min','Rainfall_Intensity_mm_h'});
writetable(rain, fullfile(out_dir, 'Rainfall_Used.csv'));

outlets = table(row_out, col_out, elevation(outlet_index), ...
    'VariableNames', {'row','col','elevation_m'});
writetable(outlets, fullfile(out_dir, 'Outlet_Cells.csv'));

if isfield(Input_Rainfall, 'time')
    manifest = table(Input_Rainfall.time(:), string(Input_Rainfall.labels_Directory(:)), ...
        'VariableNames', {'Time_min','File'});
    writetable(manifest, fullfile(out_dir, 'Rainfall_Raster_Manifest.csv'));
end
end

function write_hydrograph_outputs(out_dir, running_control, outlet_states, ...
    outlet_runoff_volume, depths, C_a, BC_States, Wshed_Properties, ...
    Rainfall_Parameters, elevation)
time_min = running_control.time_hydrograph(:);
q_m3s = outlet_states.outlet_hydrograph(:);
keep = isfinite(time_min) & isfinite(q_m3s);
time_min = time_min(keep);
q_m3s = q_m3s(keep);
hydrograph = table(time_min, q_m3s, 'VariableNames', {'Time_min','Q_m3s'});
writetable(hydrograph, fullfile(out_dir, 'FullConfig_Hydrograph.csv'));

final_storage_m3 = sum(C_a .* depths.d_t ./ 1000, 'all', 'omitnan');
rain_volume_m3 = rainfall_volume_from_hyetograph(Rainfall_Parameters, ...
    sum(C_a(~isnan(elevation)), 'all', 'omitnan'), running_control.routing_time);
outlet_volume_m3 = outlet_runoff_volume / 1000 * Wshed_Properties.drainage_area;
mass_residual_m3 = final_storage_m3 + outlet_volume_m3 - rain_volume_m3;
mass_error_pct = 100 * mass_residual_m3 / max(rain_volume_m3, eps);
summary = table(rain_volume_m3, outlet_volume_m3, final_storage_m3, ...
    mass_residual_m3, mass_error_pct, ...
    'VariableNames', {'rain_volume_m3','outlet_volume_m3','final_storage_m3', ...
    'mass_residual_m3','mass_error_pct'});
writetable(summary, fullfile(out_dir, 'FullConfig_Mass_Summary.csv'));
end

function rain_volume_m3 = rainfall_volume_from_hyetograph(Rainfall_Parameters, area_m2, routing_time_min)
t = double(Rainfall_Parameters.time_rainfall(:));
p = double(Rainfall_Parameters.intensity_rainfall(:));
if isempty(t) || isempty(p)
    rain_volume_m3 = 0;
    return
end
[t, order] = sort(t);
p = p(order);
if numel(t) < 2
    rain_depth_m = 0;
else
    t_next = min(t(2:end), routing_time_min);
    duration_h = max(t_next - t(1:end-1), 0) ./ 60;
    rain_depth_m = sum(p(1:end-1) .* duration_h, 'omitnan') ./ 1000;
end
rain_volume_m3 = rain_depth_m .* area_m2;
end
