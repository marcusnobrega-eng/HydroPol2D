%% ========================================================================
% HydroPol2D SNISB dam-break bypass configuration
% ========================================================================
% This bypass file is used by run_snisb_hydropol2d_batch.m. It keeps the
% HydroPol2D setup intentionally narrow for the first national-dam workflow:
%   - local inertial solver by default; full momentum can be selected with
%     Overrides.RoutingModel = 'full_momentum'
%   - no rainfall
%   - no infiltration / ETP
%   - inflow boundary from SNISB dam-break inlet-cell hydrograph CSV
%   - 2-day routing period
%   - map outputs every 1 minute
%
% The inflow reader accepts either:
%   1) the native SNISB long CSV:
%      time_s,inlet_cell_id,x_m,y_m,inflow_m3s
%   2) the HydroPol2D wide CSV:
%      time_min,Q_1,x_1,y_1,Q_2,x_2,y_2,...
% ========================================================================

InputData_Bypass = struct();

if ~exist('InputPaths','var') || ~isstruct(InputPaths)
    error(['InputPaths was not found. Run input_paths_bypass(...) before ' ...
           'input_data_bypass_snisb.m.']);
end

%% General controls
InputData_Bypass.general = struct();

InputData_Bypass.general.time_step_model      = 5/60;    % [min] 5 s
InputData_Bypass.general.min_time_step        = 0.01;    % [s]
InputData_Bypass.general.max_time_step        = 10;      % [s]
InputData_Bypass.general.time_step_increments = 0.001;   % [s]
InputData_Bypass.general.time_step_change     = 0.0001;

InputData_Bypass.general.alfa_min = 0.50;
InputData_Bypass.general.alfa_max = 0.50;

simulation_minutes = 2 * 24 * 60;     % [min] production default: 2 days
record_time_maps_minutes = 1;          % [min] production default
routing_model = 'local_inertial';
if exist('Overrides','var') && isstruct(Overrides)
    if isfield(Overrides,'SimulationMinutes') && ~isempty(Overrides.SimulationMinutes)
        simulation_minutes = double(Overrides.SimulationMinutes);
    end
    if isfield(Overrides,'RecordTimeMapsMinutes') && ~isempty(Overrides.RecordTimeMapsMinutes)
        record_time_maps_minutes = double(Overrides.RecordTimeMapsMinutes);
    end
    if isfield(Overrides,'RoutingModel') && ~isempty(Overrides.RoutingModel)
        routing_model = lower(strtrim(char(Overrides.RoutingModel)));
    end
end
if ~any(strcmp(routing_model, {'local_inertial','full_momentum'}))
    error('Unsupported SNISB routing model: %s', routing_model);
end

InputData_Bypass.general.date_begin  = datetime(2025,1,1,0,0,0);
InputData_Bypass.general.date_end    = InputData_Bypass.general.date_begin + minutes(simulation_minutes);
InputData_Bypass.general.routing_time = simulation_minutes;

InputData_Bypass.general.dt_rainfall_maps_min      = 60;
InputData_Bypass.general.dt_transpiration_maps_min = 1440;
InputData_Bypass.general.dt_evaporation_maps_min   = 1440;
InputData_Bypass.general.rainfall_filename_example = 'rain_2025_01_01_00_00.tif';

InputData_Bypass.general.slope_outlet                 = 0.02;
InputData_Bypass.general.n_outlets_data               = 1;
if isfield(InputPaths,'Outlet_Cells_CSV') && exist(InputPaths.Outlet_Cells_CSV,'file') == 2
    try
        outlet_table_for_count = readtable(InputPaths.Outlet_Cells_CSV);
        InputData_Bypass.general.n_outlets_data = max(1, height(outlet_table_for_count));
    catch
        InputData_Bypass.general.n_outlets_data = 1;
    end
end
InputData_Bypass.general.record_time_maps             = record_time_maps_minutes;
InputData_Bypass.general.record_time_hydrographs      = 1;    % [min]
InputData_Bypass.general.Pol_min                      = 0.001;
InputData_Bypass.general.depth_wse                    = 0.01;
InputData_Bypass.general.flag_wse                     = 0;
InputData_Bypass.general.record_time_spatial_rainfall = 60;
InputData_Bypass.general.time_save_ETP                = 1440;
InputData_Bypass.general.record_time_spatial_ETP      = 1440;
InputData_Bypass.general.Krs_ETP                      = 0.13;
InputData_Bypass.general.albedo                       = 0.25;

InputData_Bypass.general.alfa_1  = 1.0;
InputData_Bypass.general.alfa_2  = 1.0;
InputData_Bypass.general.beta_1  = 1.0;
InputData_Bypass.general.beta_2  = 1.0;
InputData_Bypass.general.Manning = 0.035;
if exist('Overrides','var') && isstruct(Overrides) && ...
        isfield(Overrides,'Manning') && ~isempty(Overrides.Manning)
    InputData_Bypass.general.Manning = double(Overrides.Manning);
end

InputData_Bypass.general.ADD    = 5.0;
InputData_Bypass.general.min_Bt = 0.0;
InputData_Bypass.general.Bmin   = 0.0;
InputData_Bypass.general.Bmax   = 1.0;

InputData_Bypass.general.min_area            = 0.5;
InputData_Bypass.general.tau                 = 0.2;
InputData_Bypass.general.K_value             = 10;
InputData_Bypass.general.sl                  = 0.001;
InputData_Bypass.general.slope_DTM           = 0.05;
InputData_Bypass.general.resolution_resample = 30;
InputData_Bypass.general.topo_path           = InputPaths.topo_path;

% Required by input_data_script even when design storms are disabled.
InputData_Bypass.general.RP                = 10;
InputData_Bypass.general.Rainfall_Duration = 180;
InputData_Bypass.general.K                 = 1000;
InputData_Bypass.general.a                 = 0.20;
InputData_Bypass.general.b                 = 10;
InputData_Bypass.general.c                 = 0.80;
InputData_Bypass.general.dt_design         = 5;

%% Flags
InputData_Bypass.flags = struct();

InputData_Bypass.flags.flag_rainfall                     = 0;
InputData_Bypass.flags.flag_spatial_rainfall             = 0;
InputData_Bypass.flags.flag_ETP                          = 0;
InputData_Bypass.flags.flag_input_rainfall_map           = 0;
InputData_Bypass.flags.flag_rainfall_multiple_runs       = 0;
InputData_Bypass.flags.flag_data_source                  = 0;
InputData_Bypass.flags.flag_inflow                       = 1;
InputData_Bypass.flags.flag_satellite_rainfall           = 0;
InputData_Bypass.flags.flag_alternated_blocks            = 0;
InputData_Bypass.flags.flag_huff                         = 0;
InputData_Bypass.flags.flag_stage_hydrograph             = 0;
InputData_Bypass.flags.flag_input_ETP_map                = 0;

InputData_Bypass.flags.flag_timestep                     = 1;
InputData_Bypass.flags.flag_infiltration                 = 1;
InputData_Bypass.flags.flag_critical                     = 0;
InputData_Bypass.flags.flag_D8                           = 0;
InputData_Bypass.flags.flag_CA                           = 0;
InputData_Bypass.flags.flag_inertial                     = double(strcmp(routing_model, 'local_inertial'));
InputData_Bypass.flags.flag_full_momentum                = double(strcmp(routing_model, 'full_momentum'));
InputData_Bypass.flags.flag_waterbalance                 = 0;
InputData_Bypass.flags.flag_waterquality                 = 0;
InputData_Bypass.flags.flag_reservoir                    = 0;
InputData_Bypass.flags.flag_wq_model                     = 0;
InputData_Bypass.flags.flag_groundwater_modeling         = 0;
InputData_Bypass.flags.flag_real_time_satellite_rainfall = 0;
InputData_Bypass.flags.flag_dam_break                    = 0;
InputData_Bypass.flags.flag_human_instability            = 0;
InputData_Bypass.flags.flag_boundary                     = 0;
InputData_Bypass.flags.flag_numerical_scheme             = 1;
InputData_Bypass.flags.flag_outlet_type                  = 1;
InputData_Bypass.flags.flag_adaptive_timestepping        = 1;
InputData_Bypass.flags.flag_neglect_infiltration_river   = 0;
InputData_Bypass.flags.flag_subgrid                      = 0;
InputData_Bypass.flags.flag_spatial_albedo               = 0;
InputData_Bypass.flags.flag_river_rasters                = 0;
InputData_Bypass.flags.flag_baseflow                     = 0;
InputData_Bypass.flags.flag_kinematic                    = 0;
InputData_Bypass.flags.flag_diffusive                    = 0;
InputData_Bypass.flags.flag_DTM                          = 0;
InputData_Bypass.flags.flag_abstraction                  = 0;
InputData_Bypass.flags.flag_overbanks                    = 0;
InputData_Bypass.flags.flag_snow_modeling                = 0;
InputData_Bypass.flags.flag_WQ_Rasters                   = 0;
InputData_Bypass.flags.flag_GPU                          = 0;
InputData_Bypass.flags.flag_single                       = 0;
InputData_Bypass.flags.flag_warmup                       = 0;
InputData_Bypass.flags.flag_initial_buildup              = 0;
InputData_Bypass.flags.flag_resample                     = 0;
InputData_Bypass.flags.flag_smoothening                  = 0;
InputData_Bypass.flags.flag_trunk                        = 0;
InputData_Bypass.flags.flag_fill_DEM                     = 0;
InputData_Bypass.flags.flag_smooth_cells                 = 0;
InputData_Bypass.flags.flag_reduce_DEM                   = 0;
InputData_Bypass.flags.flag_export_maps                  = 1;
InputData_Bypass.flags.flag_river_heigth_compensation    = 0;
InputData_Bypass.flags.flag_dashboard                    = 0;
InputData_Bypass.flags.flag_elapsed_time                 = 1;
InputData_Bypass.flags.flag_obs_gauges                   = 0;

if exist('Overrides','var') && isstruct(Overrides)
    if isfield(Overrides,'EnableInflow') && ~isempty(Overrides.EnableInflow)
        InputData_Bypass.flags.flag_inflow = double(logical(Overrides.EnableInflow));
    end
    if isfield(Overrides,'EnableRainfall') && ~isempty(Overrides.EnableRainfall)
        InputData_Bypass.flags.flag_rainfall = double(logical(Overrides.EnableRainfall));
    end
    if isfield(Overrides,'EnableInfiltration') && ~isempty(Overrides.EnableInfiltration)
        InputData_Bypass.flags.flag_infiltration = double(logical(Overrides.EnableInfiltration));
    end
    if isfield(Overrides,'EnableWarmup') && ~isempty(Overrides.EnableWarmup)
        InputData_Bypass.flags.flag_warmup = double(logical(Overrides.EnableWarmup));
    end
    if isfield(Overrides,'EnablePerimeterOutlet') && ~isempty(Overrides.EnablePerimeterOutlet)
        InputData_Bypass.flags.flag_boundary = double(logical(Overrides.EnablePerimeterOutlet));
    end
end

%% MapBiomas LULC roughness table
LULC_table = table();
LULC_table.LC = { ...
    'Unclassified'; 'Forest Formation'; 'Savanna Formation'; 'Mangrove'; ...
    'Forest Plantation'; 'Wetland'; 'Grassland'; 'Pasture'; 'Agriculture'; ...
    'Temporary Crops'; 'Sugar Cane'; 'Mosaic of Uses'; 'Beach/Dune/Sand'; ...
    'Urban Area'; 'Other Non Vegetated'; 'Rocky Outcrop'; 'Mining'; ...
    'Aquaculture'; 'Salt Flat'; 'River/Lake/Ocean'; 'Soybean'; 'Rice'; ...
    'Other Temporary Crops'; 'Coffee'; 'Citrus'; 'Other Perennial Crops'; ...
    'Wooded Sandbank'; 'Herbaceous Sandbank'; 'Cotton'; 'Photovoltaic'};
LULC_table.Index = [0;3;4;5;9;11;12;15;18;19;20;21;23;24;25;29;30;31;32;33;39;40;41;46;47;48;49;50;62;75];
LULC_table.roughness = [0.050;0.120;0.080;0.150;0.100;0.100;0.055;0.050;0.045;0.045;0.060;0.055;0.030;0.030;0.035;0.040;0.035;0.035;0.035;0.035;0.045;0.060;0.045;0.080;0.080;0.080;0.100;0.080;0.045;0.030];
LULC_table.h_0_mm = zeros(height(LULC_table),1);
LULC_table.d_0_mm = zeros(height(LULC_table),1);
LULC_table.C1 = zeros(height(LULC_table),1);
LULC_table.C2 = zeros(height(LULC_table),1);
LULC_table.C3 = zeros(height(LULC_table),1);
LULC_table.C4 = zeros(height(LULC_table),1);
LULC_table.index_impervious = [24; nan(height(LULC_table)-1,1)];
InputData_Bypass.LULC.table = LULC_table;

%% OpenLandMap USDA texture soil table
SOIL_table = table();
SOIL_table.Soil_type = { ...
    'Clay'; 'Silty Clay'; 'Sandy Clay'; 'Clay Loam'; 'Silty Clay Loam'; ...
    'Sandy Clay Loam'; 'Silty Loam'; 'Loam'; 'Sandy Loam'; 'Silt'; ...
    'Loamy Sand'; 'Sand'; 'Water'};
SOIL_table.Index = [1;2;3;4;5;6;7;8;9;10;11;12;0];
SOIL_table.ksat_mm_h     = [0.3;0.5;0.6;1.0;1.0;1.5;7.6;3.4;10.9;3.4;29.9;117.8;0.3];
SOIL_table.n_vg          = [1.09;1.23;1.31;1.31;1.23;1.48;1.68;1.56;1.89;1.37;1.75;2.68;1.09];
SOIL_table.alpha_vg_1_m  = [0.8;1.0;1.5;1.9;1.5;3.0;2.0;3.6;6.0;1.6;11.0;14.5;0.8];
SOIL_table.theta_sat     = [0.385;0.423;0.321;0.309;0.432;0.330;0.432;0.399;0.387;0.481;0.390;0.430;0.385];
SOIL_table.theta_r       = [0.068;0.089;0.075;0.095;0.089;0.065;0.067;0.078;0.100;0.034;0.049;0.045;0.068];
SOIL_table.theta_i       = [0.36;0.41;0.27;0.29;0.38;0.26;0.33;0.26;0.21;0.35;0.12;0.10;0.36];  % field capacity (theta_fc)
SOIL_table.Sy            = [0.03;0.04;0.06;0.08;0.09;0.11;0.14;0.18;0.22;0.12;0.25;0.30;0.03];
SOIL_table.ksat_gw_mm_h  = [6;10;12;20;20;30;152;68;218;68;598;2356;6];
SOIL_table.Soil_Depth_m  = ones(height(SOIL_table),1);
InputData_Bypass.SOIL.table = SOIL_table;

%% Inflow hydrograph
if InputData_Bypass.flags.flag_inflow == 1
    if ~isfield(InputPaths,'Inflow_Hydrograph_CSV') || isempty(InputPaths.Inflow_Hydrograph_CSV)
        error('InputPaths.Inflow_Hydrograph_CSV is required when flag_inflow = 1.');
    end
    InputData_Bypass.Inflow = read_snisb_or_wide_inflow(InputPaths.Inflow_Hydrograph_CSV);
else
    InputData_Bypass.Inflow = struct( ...
        'time_inflow', [], ...
        'inflow_discharge', [], ...
        'x_coords', {{}}, ...
        'y_coords', {{}});
end

InputData_Bypass.Stage = struct( ...
    'time_stage', [], ...
    'stage', [], ...
    'x_coords', {{}}, ...
    'y_coords', {{}});

%% Lumped rainfall hyetograph
Rainfall_Parameters = struct();
rainfall_file = '';
if isfield(InputPaths,'Rainfall_Timeseries_File') && ~isempty(InputPaths.Rainfall_Timeseries_File)
    rainfall_file = InputPaths.Rainfall_Timeseries_File;
end
if isempty(rainfall_file) && exist('Overrides','var') && isstruct(Overrides) && ...
        isfield(Overrides,'Rainfall_Timeseries_File') && ~isempty(Overrides.Rainfall_Timeseries_File)
    rainfall_file = Overrides.Rainfall_Timeseries_File;
end

if InputData_Bypass.flags.flag_rainfall == 1
    if isempty(rainfall_file) || exist(rainfall_file,'file') ~= 2
        error('Rainfall_Timeseries_File is required when flag_rainfall = 1. Tried: %s', rainfall_file);
    end
    T_rain = readtable(rainfall_file);
    if width(T_rain) < 2
        error('Rainfall time-series file must contain at least two columns: [time_min, intensity_mm_h].');
    end
    Rainfall_Parameters.time_rainfall = double(table2array(T_rain(:,1)));
    Rainfall_Parameters.intensity_rainfall = double(table2array(T_rain(:,2)));
    Rainfall_Parameters.time_rainfall = Rainfall_Parameters.time_rainfall(:);
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(:);
    valid = isfinite(Rainfall_Parameters.time_rainfall) & ...
            isfinite(Rainfall_Parameters.intensity_rainfall) & ...
            Rainfall_Parameters.intensity_rainfall >= 0;
    Rainfall_Parameters.time_rainfall = Rainfall_Parameters.time_rainfall(valid);
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(valid);
    validate_constant_timestep(Rainfall_Parameters.time_rainfall, 'Rainfall hyetograph');
    Rainfall_Parameters.time_step_rainfall = Rainfall_Parameters.time_rainfall(2) - Rainfall_Parameters.time_rainfall(1);
    Rainfall_Parameters.rainfall_duration = Rainfall_Parameters.time_rainfall(end);
    Rainfall_Parameters.n_obs_rainfall = numel(Rainfall_Parameters.time_rainfall);
else
    Rainfall_Parameters.time_rainfall = [0; 1];
    Rainfall_Parameters.intensity_rainfall = [0; 0];
    Rainfall_Parameters.time_step_rainfall = 1;
    Rainfall_Parameters.rainfall_duration = 1;
    Rainfall_Parameters.n_obs_rainfall = 2;
end
InputData_Bypass.Rainfall_Parameters = Rainfall_Parameters;

InputData_Bypass.Reservoir_Data = struct( ...
    'index', [], 'x_us', [], 'y_us', [], ...
    'k1', [], 'h1', [], 'k2', [], 'x_ds1', [], 'y_ds1', [], ...
    'k3', [], 'h2', [], 'k4', [], 'x_ds2', [], 'y_ds2', []);

InputData_Bypass.Obs = table();

%% Local readers
function Inflow = read_snisb_or_wide_inflow(csv_path)
    if exist(csv_path,'file') ~= 2
        error('SNISB inflow hydrograph CSV not found: %s', csv_path);
    end

    T = readtable(csv_path);
    vars = T.Properties.VariableNames;

    if all(ismember({'time_s','inlet_cell_id','x_m','y_m','inflow_m3s'}, vars))
        Inflow = read_snisb_long_inflow(T);
    elseif ismember('time_min', vars) && any(startsWith(vars, 'Q_'))
        Inflow = read_hydropol_wide_inflow(T);
    else
        error(['Unsupported inflow CSV format: %s\n' ...
               'Expected either SNISB long columns or HydroPol2D wide columns.'], csv_path);
    end
end

function Inflow = read_snisb_long_inflow(T)
    time_s = double(T.time_s(:));
    ids = unique(T.inlet_cell_id, 'stable');
    times = unique(time_s, 'stable');
    times = sort(times(:));

    n_times = numel(times);
    n_ids = numel(ids);
    Q = zeros(n_times, n_ids);
    x_coords = cell(1, n_ids);
    y_coords = cell(1, n_ids);

    for i = 1:n_ids
        id = ids(i);
        idx = T.inlet_cell_id == id;
        Ti = T(idx,:);
        [tf, loc] = ismember(double(Ti.time_s(:)), times);
        Q(loc(tf), i) = double(Ti.inflow_m3s(tf));
        x_coords{i} = double(Ti.x_m(1));
        y_coords{i} = double(Ti.y_m(1));
    end

    validate_constant_timestep(times / 60, 'SNISB inflow hydrograph');
    Inflow.time_inflow = times / 60;
    Inflow.inflow_discharge = Q;
    Inflow.x_coords = x_coords;
    Inflow.y_coords = y_coords;
end

function Inflow = read_hydropol_wide_inflow(T)
    vars = T.Properties.VariableNames;
    q_names = vars(startsWith(vars, 'Q_'));
    x_names = vars(startsWith(vars, 'x_'));
    y_names = vars(startsWith(vars, 'y_'));
    n = numel(q_names);
    if ~(n == numel(x_names) && n == numel(y_names))
        error('Wide inflow CSV must have matching Q_i, x_i, and y_i columns.');
    end

    time_min = double(T.time_min(:));
    validate_constant_timestep(time_min, 'HydroPol2D wide inflow hydrograph');

    Inflow.time_inflow = time_min;
    Inflow.inflow_discharge = double(table2array(T(:, q_names)));
    Inflow.x_coords = cell(1, n);
    Inflow.y_coords = cell(1, n);
    for i = 1:n
        Inflow.x_coords{i} = double(table2array(T(1, x_names(i))));
        Inflow.y_coords{i} = double(table2array(T(1, y_names(i))));
    end
end

function validate_constant_timestep(time_min, label)
    time_min = double(time_min(:));
    if numel(time_min) < 2
        error('%s must contain at least two time steps.', label);
    end
    dt = diff(time_min);
    if any(dt <= 0)
        error('%s times must be strictly increasing.', label);
    end
    if any(abs(dt - dt(1)) > 1.0e-7)
        error('%s must have a constant time step.', label);
    end
end
