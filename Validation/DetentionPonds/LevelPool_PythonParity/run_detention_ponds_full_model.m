% VAL-POND-001 full-model run (one variant): normal HydroPol2D preprocessing and main
% loop in bypass mode, with two detention ponds on a tilted plane.
%
% Called by run_detention_ponds_check.m with ValidationCase.variant set
% ('groundwater' or 'no_groundwater'). Kept as a separate script because
% preprocessing clears the workspace (only the variables it lists survive, among
% them ValidationCase, Cases, icase, case_dir and summary_dir).
%
% Domain: 12 x 16 cells of 10 m, falling 0.01 to the south and 0.004 toward the
% middle column. Pond A (DTM, unlined, 10 mm/h bottom Ksat) is a 3 x 4 basin dug
% 0.6 m into rows 5-7, columns 7-10 (polygon footprint); it has an orifice and a
% spillway onto the ditch south of it (columns 8-9, 0.6 m deep, to the outlet) and a
% pump that sends water out of the model.
% Pond B (power law, lined) covers rows 2-3, columns 13-14 (cell-list footprint),
% with one orifice. Rain 60 mm/h for 20 min; 'no_groundwater': station ETP and
% infiltration, groundwater off; 'groundwater': the legacy groundwater module (local
% update, no lateral flow) with infiltration off and a dry soil.

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
addpath(fullfile(repo_root, 'HydroPol2D_Functions'), '-begin');
hydropol2d_add_runtime_paths(hydropol2d_find_root(case_dir));

variant = ValidationCase.variant;
output_root = fullfile(case_dir, 'Outputs', 'FullModel', variant);
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
static_dir = fullfile(output_root, 'Static');
if exist(output_root, 'dir'); rmdir(output_root, 's'); end
mkdir(static_dir);
if ~exist(summary_dir, 'dir'); mkdir(summary_dir); end

% Static rasters.
ny = 12; nx = 16; dx = 10;
[C, R] = meshgrid(1:nx, 1:ny);
DEM = 100 + 0.01 * dx * (ny - R) + 0.004 * dx * abs(C - 8.5);
pond_a_cells = R >= 5 & R <= 7 & C >= 7 & C <= 10;
DEM(pond_a_cells) = DEM(pond_a_cells) - 0.6 + 0.02 * (C(pond_a_cells) - 8.5);
ditch = R >= 8 & (C == 8 | C == 9);       % takes the pond's releases to the outlet
DEM(ditch) = DEM(ditch) - 0.6;
georef = maprefcells([0, nx * dx], [0, ny * dx], [ny, nx], 'ColumnsStartFrom', 'north');
% Float64 rasters: with float32 heads hp2d_boussinesq_adaptive can stall (a single
% suggested step one ulp short of the interval is retaken forever; see README.md).
write_raster = @(name, Z) geotiffwrite(fullfile(static_dir, [name '.tif']), double(Z), georef, ...
    'CoordRefSysCode', 3857);
write_raster('DEM', DEM);
write_raster('LULC', ones(ny, nx));
write_raster('SOIL', ones(ny, nx));
write_raster('DTB', 2.0 * ones(ny, nx));
write_raster('GW_table', DEM - 1.2);
write_raster('LAI', zeros(ny, nx));
write_raster('Albedo', 0.23 * ones(ny, nx));
write_raster('Initial_SM', 30 * ones(ny, nx));
etp_path = fullfile(output_root, 'ETP_input_data.xlsx');
write_station_forcing(etp_path);

% Pond geometry for the bypass configuration (map coordinates of the cell centres:
% x = (column - 0.5) dx, y = (ny - row + 0.5) dx).
ValidationCase.pond = struct( ...
    'a_vertices_m', [62.5 52.5; 97.5 52.5; 97.5 77.5; 62.5 77.5], ...
    'a_invert_m', min(DEM(pond_a_cells)), ...
    'overflow_x_m', (9 - 0.5) * dx, 'overflow_y_m', (ny - 8 + 0.5) * dx, ...
    'b_cells', [2 13; 2 14; 3 13; 3 14], ...
    'b_invert_m', min(DEM(R >= 2 & R <= 3 & C >= 13 & C <= 14)) - 0.3);

Paths = make_output_paths(output_root);
InputPaths = make_input_paths(static_dir, etp_path);
input_data_bypass_script_path = fullfile(case_dir, 'Config', 'input_data_bypass_script.m');
use_inputpaths_bypass = 1;
use_inputdata_bypass = 1;
clean_output_folder = true;
run_postprocessing = false;
enable_logging = false;
model_folder = '';
resultsDir = Paths.Results;
export_root_dir = Paths.Root;

HydroPol2D_preprocessing;
running_control.flag_system_mass_ledger = true;
HydroPol2D_Main_While;

% ------------------------------------------------------------------------
% Checks
% ------------------------------------------------------------------------
variant = ValidationCase.variant;
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
L = system_mass_ledger;
P = Detention_Ponds;
if ~(isstruct(P) && isfield(P, 'storage_m3'))
    % Baseline without ponds: the model's own ledger residual on this domain.
    Row = table(string(variant), L.cumulative_precipitation_m3, L.cumulative_outlet_m3, ...
        L.cumulative_residual_m3, 'VariableNames', {'variant','rain_m3','outlet_m3', ...
        'ledger_cumulative_residual_m3'});
    out_path = fullfile(summary_dir, 'FullModel_Baseline_Without_Ponds.csv');
    if isfile(out_path)
        T = readtable(out_path, 'TextType', 'string');
        T(T.variant == Row.variant, :) = [];
        Row = [T; Row];
    end
    writetable(Row, out_path);
    disp(Row);
    return
end
rain_m3 = L.cumulative_precipitation_m3;
losses_m3 = L.cumulative_canopy_evaporation_m3 + L.cumulative_etr_m3 + ...
    L.cumulative_open_water_evaporation_m3 + L.cumulative_snow_sublimation_m3 + ...
    L.cumulative_outlet_m3 + L.cumulative_pond_evaporation_m3 + L.cumulative_pond_seepage_out_m3;
inputs_m3 = rain_m3 + L.cumulative_boundary_inflow_m3 + L.cumulative_prescribed_recharge_m3;
system_residual_m3 = (L.final_storage_m3 - L.initial_storage_m3) - (inputs_m3 - losses_m3);
% The same ledger without the pond terms (pond storage and pond losses left out).
no_pond_residual_m3 = ((L.final_storage_m3 - sum(P.storage_m3)) - ...
    (L.initial_storage_m3 - sum(P.initial_storage_m3))) - ...
    (inputs_m3 - losses_m3 + L.cumulative_pond_evaporation_m3 + L.cumulative_pond_seepage_out_m3);

device_net = accumarray(P.dev_pond, P.cumulative_device_m3, [P.n_ponds 1]);
pond_residual_m3 = P.storage_m3 - P.initial_storage_m3 - P.cumulative_rainfall_m3 - ...
    P.cumulative_lateral_m3 - P.cumulative_surface_in_m3 + P.cumulative_evaporation_m3 + ...
    P.cumulative_seepage_m3 + device_net;
pump = P.dev_kind == 2;
pump_external_error_m3 = abs(sum(P.cumulative_device_m3(pump)) - sum(P.cumulative_external_m3));
pond_cell_surface_mm = max(abs(double(gather(depths.d_t(P.cells)))));
pond_cell_soil_mm = max(abs(double(gather(Soil_Properties.I_t(P.cells)))));
seepage_m3 = sum(P.cumulative_seepage_m3);
seepage_out_m3 = sum(P.cumulative_seepage_out_m3);
if strcmp(variant, 'groundwater')
    seepage_route_error_m3 = seepage_out_m3;           % all of it went to the aquifer
else
    seepage_route_error_m3 = abs(seepage_out_m3 - seepage_m3);   % all of it left
end

Row = table(string(variant), rain_m3, L.cumulative_outlet_m3, sum(P.initial_storage_m3), ...
    sum(P.storage_m3), sum(P.cumulative_lateral_m3), sum(P.cumulative_rainfall_m3), ...
    sum(P.cumulative_evaporation_m3), seepage_m3, seepage_out_m3, ...
    sum(P.cumulative_surface_in_m3), P.cumulative_device_m3(1), P.cumulative_device_m3(2), ...
    P.cumulative_device_m3(3), P.cumulative_device_m3(4), sum(P.cumulative_external_m3), ...
    nnz(diff([0; P.history.device_discharge_m3_s(:, 3) > 0]) == 1), ...
    system_residual_m3, L.cumulative_residual_m3, no_pond_residual_m3, max(abs(pond_residual_m3)), ...
    pump_external_error_m3, seepage_route_error_m3, pond_cell_surface_mm, pond_cell_soil_mm, ...
    'VariableNames', {'variant','rain_m3','outlet_m3','pond_initial_storage_m3', ...
    'pond_final_storage_m3','critical_inflow_m3','pond_rain_m3','pond_evaporation_m3', ...
    'pond_seepage_m3','pond_seepage_out_m3','pond_cell_surface_in_m3', ...
    'low_orifice_m3','overflow_m3','lift_pump_m3','b_orifice_m3','external_device_m3', ...
    'pump_starts','system_residual_m3','ledger_cumulative_residual_m3', ...
    'residual_without_pond_terms_m3','pond_residual_m3','pump_external_error_m3', ...
    'seepage_route_error_m3','pond_cell_surface_mm','pond_cell_soil_mm'});
% The ledger closes to the model's own residual on this domain: with infiltration and
% ETP (N) ~3e-6 of the rain without ponds, from the legacy soil bookkeeping; without them
% (G) to round-off (2e-13 m3 without ponds). The pond terms are 1-2% of the rain, so
% leaving any of them out fails by orders of magnitude.
ledger_tolerance = 1e-5;
if flags.flag_infiltration == 0
    ledger_tolerance = 1e-9;
end
Row.passed = abs(Row.system_residual_m3) <= ledger_tolerance * max(Row.rain_m3, 1) && ...
    Row.pond_residual_m3 <= 1e-9 && Row.pump_external_error_m3 <= 1e-12 && ...
    Row.seepage_route_error_m3 <= 1e-12 && Row.pond_cell_surface_mm == 0 && ...
    Row.pond_cell_soil_mm == 0 && all(P.cumulative_lateral_m3 > 0) && ...
    Row.pond_seepage_m3 > 0 && (Row.pond_evaporation_m3 > 0 || flags.flag_ETP ~= 1) && ...
    all(abs(P.cumulative_device_m3) > 0);
out_path = fullfile(summary_dir, 'FullModel_Summary.csv');
if isfile(out_path)
    T = readtable(out_path, 'TextType', 'string');
    T(T.variant == Row.variant, :) = [];
    T = [T; Row];
else
    T = Row;
end
writetable(T, out_path);
copyfile(fullfile(Paths.Tables, 'Detention_Ponds_Time_Series.csv'), ...
    fullfile(summary_dir, sprintf('FullModel_%s_Pond_Time_Series.csv', variant)));
disp(Row);

function InputPaths = make_input_paths(static_dir, etp_path)
InputPaths = struct();
InputPaths.DEM_path = fullfile(static_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_dir, 'SOIL.tif');
InputPaths.DTB_path = fullfile(static_dir, 'DTB.tif');
InputPaths.GW_table_path = fullfile(static_dir, 'GW_table.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_dir, 'Initial_SM.tif');
InputPaths.LAI_path = fullfile(static_dir, 'LAI.tif');
InputPaths.Albedo_path = fullfile(static_dir, 'Albedo.tif');
InputPaths.Initial_SWE_path = '';
InputPaths.Initial_Snow_Depth_path = '';
InputPaths.ETP_input_spreadsheet = etp_path;
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
InputPaths.Inflow_Hydrograph_CSV = '';
InputPaths.Stage_Hydrograph_CSV = '';
InputPaths.Observed_Gauges_CSV = '';
InputPaths.Rainfall_Timeseries_File = '';
InputPaths.Outlet_Cells_CSV = '';
end

function Paths = make_output_paths(root_dir)
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

function write_station_forcing(path)
% Warm, dry summer days so the ponds evaporate (Penman open-water Ep). The ground
% heat flux is 0, as for daily data (ETP_model reads it, as HydroPol2D-Python does).
times = datetime(2025, 7, 15, 0, 0, 0) + hours([0; 24; 48]);
data = [32 18 25 3 35 0; 32 18 25 3 35 0; 32 18 25 3 35 0];
cells = cell(6, 8);
cells(1,:) = {'Index','Date','Tmax_C','Tmin_C','Tavg_C','U2_m_s','RH_pct','G_MJ_m2_day'};
cells(2,:) = {NaN,NaT,NaN,NaN,NaN,0,NaN,0};
cells(3,:) = {NaN,NaT,NaN,NaN,NaN,NaN,NaN,NaN};
for i = 1:numel(times)
    cells{3+i,1} = i;
    cells{3+i,2} = times(i);
    cells(3+i,3:8) = num2cell(data(i,:));
end
writecell(cells, path);
end
