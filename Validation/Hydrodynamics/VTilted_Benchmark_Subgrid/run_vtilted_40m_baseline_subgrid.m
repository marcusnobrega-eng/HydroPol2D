clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
model_root = fullfile(repo_root, 'HydroPol2D_Model');
functions_dir = fullfile(model_root, 'HydroPol2D_Functions');
topo_path = fullfile(model_root, 'topotoolbox-master');
static_root = fullfile(case_dir, 'Static_40m_Comparison');

scenario = canonical_scenario_from_env();
out_dir = fullfile(case_dir, 'Outputs', 'VTilted_40m_Comparison', char(scenario));

addpath(genpath(functions_dir));
addpath(genpath(topo_path));
addpath('/Users/mngomes/Downloads', '-end');

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

Paths = make_validation_paths(fullfile(out_dir, 'HydroPol2D_Output'), true);
InputPaths = make_input_paths(static_root, scenario, topo_path, functions_dir);

input_data_bypass_script_path = fullfile(case_dir, 'input_data_bypass_script_vtilted_40m_compat.m');

use_inputpaths_bypass = 1;
use_inputdata_bypass = 1;
clean_output_folder = true;
run_postprocessing = false;
enable_logging = false;
model_folder = '';
GD = [];
resultsDir = Paths.Results;
export_root_dir = Paths.Root;
SubgridTables = [];

HydroPol2D_preprocessing;
scenario = canonical_scenario_from_env();
out_dir = fullfile(case_dir, 'Outputs', 'VTilted_40m_Comparison', char(scenario));
input_data_bypass_script_path = fullfile(case_dir, 'input_data_bypass_script_vtilted_40m_compat.m');
write_preprocessing_audit(out_dir, scenario, InputPaths, input_data_bypass_script_path, ...
    DEM_raster, Wshed_Properties, flags, running_control, Rainfall_Parameters, ...
    Input_Rainfall, LULC_Properties, outlet_index, elevation, C_a, slope_outlet);

HydroPol2D_Main_While;
scenario = canonical_scenario_from_env();
out_dir = fullfile(case_dir, 'Outputs', 'VTilted_40m_Comparison', char(scenario));
if ~exist('SubgridTables', 'var')
    SubgridTables = [];
end
write_hydrograph_outputs(out_dir, running_control, outlet_states, outlet_runoff_volume, ...
    depths, C_a, BC_States, Wshed_Properties, Rainfall_Parameters, elevation, flags, SubgridTables);

function scenario = canonical_scenario_from_env()
scenario = lower(strtrim(string(getenv('HYDROPOL2D_VTILT_40M_SCENARIO'))));
if strlength(scenario) == 0
    scenario = "reference20m";
end
switch scenario
    case {"reference20m", "ref20", "20m"}
        scenario = "reference20m";
    case {"baseline40m", "ordinary40m", "40m"}
        scenario = "baseline40m";
    case {"subgrid40m", "sg40m"}
        scenario = "subgrid40m";
    otherwise
        error('Unsupported HYDROPOL2D_VTILT_40M_SCENARIO: %s', scenario);
end
end

function InputPaths = make_input_paths(static_root, scenario, topo_path, functions_dir)
switch scenario
    case "reference20m"
        static_dir = fullfile(static_root, 'Static_20m_Crop');
        outlet_csv = fullfile(static_root, 'Outlet', 'Outlet_20m_Crop.csv');
        subgrid_dem = '';
    case "baseline40m"
        static_dir = fullfile(static_root, 'Static_40m_Resampled');
        outlet_csv = fullfile(static_root, 'Outlet', 'Outlet_40m.csv');
        subgrid_dem = '';
    case "subgrid40m"
        static_dir = fullfile(static_root, 'Static_40m_Subgrid');
        outlet_csv = fullfile(static_root, 'Outlet', 'Outlet_40m.csv');
        subgrid_dem = fullfile(static_dir, 'Subgrid_DEM.tif');
end

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
InputPaths.Subgrid_DEM_path = subgrid_dem;
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
InputPaths.Outlet_Cells_CSV = outlet_csv;
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

function write_preprocessing_audit(out_dir, scenario, InputPaths, input_script, DEM_raster, ...
    Wshed_Properties, flags, running_control, Rainfall_Parameters, Input_Rainfall, ...
    LULC_Properties, outlet_index, elevation, C_a, slope_outlet)
[row_out, col_out] = find(outlet_index);
rain_times = Rainfall_Parameters.time_rainfall(:);
rain_intensities = Rainfall_Parameters.intensity_rainfall(:);

audit = table();
audit.item = [
    "scenario";
    "input_data_script";
    "dem_path";
    "subgrid_dem_path";
    "outlet_cells_csv";
    "rainfall_file";
    "dem_rows";
    "dem_cols";
    "resolution_m";
    "routing_time_min";
    "time_step_model_min";
    "min_time_step_s";
    "max_time_step_s";
    "flag_inertial";
    "flag_subgrid";
    "flag_overbanks";
    "flag_infiltration";
    "flag_spatial_rainfall";
    "flag_input_rainfall_map";
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
    string(scenario);
    string(input_script);
    string(InputPaths.DEM_path);
    string(InputPaths.Subgrid_DEM_path);
    string(InputPaths.Outlet_Cells_CSV);
    string(InputPaths.Rainfall_Timeseries_File);
    string(size(DEM_raster.Z, 1));
    string(size(DEM_raster.Z, 2));
    string(Wshed_Properties.Resolution);
    string(running_control.routing_time);
    string(running_control.time_step_model);
    string(running_control.min_time_step);
    string(running_control.max_time_step);
    string(flags.flag_inertial);
    string(flags.flag_subgrid);
    string(flags.flag_overbanks);
    string(flags.flag_infiltration);
    string(flags.flag_spatial_rainfall);
    string(flags.flag_input_rainfall_map);
    string(slope_outlet);
    string(numel(row_out));
    strjoin(string(row_out.'), ' ');
    strjoin(string(col_out.'), ' ');
    string(min(LULC_Properties.roughness(:), [], 'omitnan'));
    string(max(LULC_Properties.roughness(:), [], 'omitnan'));
    strjoin(string(rain_times.'), ' ');
    strjoin(string(rain_intensities.'), ' ');
    string(Wshed_Properties.drainage_area)
    ];
writetable(audit, fullfile(out_dir, 'Preprocessing_Config_Audit.csv'));

outlets = table(row_out, col_out, elevation(outlet_index), ...
    'VariableNames', {'row','col','elevation_m'});
writetable(outlets, fullfile(out_dir, 'Outlet_Cells.csv'));

rain = table(rain_times, rain_intensities, ...
    'VariableNames', {'Time_min','Rainfall_Intensity_mm_h'});
writetable(rain, fullfile(out_dir, 'Rainfall_Used.csv'));

if isfield(Input_Rainfall, 'time')
    manifest = table(Input_Rainfall.time(:), string(Input_Rainfall.labels_Directory(:)), ...
        'VariableNames', {'Time_min','File'});
    writetable(manifest, fullfile(out_dir, 'Rainfall_Raster_Manifest.csv'));
end
end

function write_hydrograph_outputs(out_dir, running_control, outlet_states, ...
    outlet_runoff_volume, depths, C_a, BC_States, Wshed_Properties, ...
    Rainfall_Parameters, elevation, flags, SubgridTables)
time_min = running_control.time_hydrograph(:);
q_m3s = outlet_states.outlet_hydrograph(:);
keep = isfinite(time_min) & isfinite(q_m3s);
time_min = time_min(keep);
q_m3s = q_m3s(keep);
hydrograph = table(time_min, q_m3s, 'VariableNames', {'Time_min','Q_m3s'});
writetable(hydrograph, fullfile(out_dir, 'Hydrograph.csv'));

if flags.flag_subgrid == 1 && exist('SubgridTables', 'var') && ~isempty(SubgridTables)
    if isfield(SubgridTables, 'sfincs_exact') && SubgridTables.sfincs_exact
        eta_final = SubgridTables.z_zmin + max(depths.d_t ./ 1000, 0);
        final_storage_m3 = sum(hp2d_sfincs_cell_volume_from_zs(SubgridTables, eta_final), ...
            'all', 'omitnan');
    else
        final_storage_m3 = sum(hp2d_subgrid_lookup_depth(SubgridTables.volume_cell, ...
            max(depths.d_t ./ 1000, 0), SubgridTables.dz, SubgridTables.maxDepth), ...
            'all', 'omitnan');
    end
else
    final_storage_m3 = sum(C_a .* depths.d_t ./ 1000, 'all', 'omitnan');
end

rain_volume_m3 = rainfall_volume_from_hyetograph(Rainfall_Parameters, ...
    Wshed_Properties.drainage_area, running_control.routing_time);
outlet_volume_m3 = outlet_runoff_volume / 1000 * Wshed_Properties.drainage_area;
mass_residual_m3 = final_storage_m3 + outlet_volume_m3 - rain_volume_m3;
mass_error_pct = 100 * mass_residual_m3 / max(rain_volume_m3, eps);
summary = table(rain_volume_m3, outlet_volume_m3, final_storage_m3, ...
    mass_residual_m3, mass_error_pct, ...
    'VariableNames', {'rain_volume_m3','outlet_volume_m3','final_storage_m3', ...
    'mass_residual_m3','mass_error_pct'});
writetable(summary, fullfile(out_dir, 'Mass_Summary.csv'));
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
