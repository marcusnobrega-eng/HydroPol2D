% V-tilted routing compatibility inputs for historical diagnostics.
% The base configuration is bundled with HydroPol2D; this script changes
% only the controlled forcing and requested routing mode.

run(fullfile(model_root, 'Config', 'input_data_bypass_script.m'));

InputData_Bypass.general.date_begin = datetime(2025, 5, 1, 0, 0, 0);
InputData_Bypass.general.date_end = InputData_Bypass.general.date_begin + minutes(240);
InputData_Bypass.general.routing_time = 240;
InputData_Bypass.general.record_time_maps = 10;
InputData_Bypass.general.record_time_hydrographs = 10;
InputData_Bypass.general.slope_outlet = 0.02;
InputData_Bypass.general.n_outlets_data = 2;
InputData_Bypass.general.max_time_step = 1;
InputData_Bypass.general.time_step_change = 0.001;
InputData_Bypass.general.alfa_min = 0.4;
InputData_Bypass.general.alfa_max = 0.4;

% This is a rain-on-grid benchmark with no imposed inflow or stage boundary.
InputData_Bypass.flags.flag_rainfall = 1;
InputData_Bypass.flags.flag_inflow = 0;
InputData_Bypass.flags.flag_stage_hydrograph = 0;
InputData_Bypass.flags.flag_resample = 0;
InputData_Bypass.flags.flag_timestep = 2;
InputData_Bypass.flags.flag_human_instability = 0;
InputData_Bypass.general.resolution_resample = 20;

% The synthetic V-tilted raster uses classes 1 (left hillslope), 2
% (channel strip), and 3 (right hillslope), rather than the global LULC
% identifiers in the generic bypass configuration.
vtilted_lulc = table();
vtilted_lulc.LC = {'Left hillslope'; 'Channel strip'; 'Right hillslope'};
vtilted_lulc.Index = [1; 2; 3];
vtilted_lulc.roughness = [0.015; 0.150; 0.015];
vtilted_lulc.h_0_mm = zeros(3, 1);
vtilted_lulc.d_0_mm = zeros(3, 1);
vtilted_lulc.C1 = zeros(3, 1);
vtilted_lulc.C2 = zeros(3, 1);
vtilted_lulc.C3 = zeros(3, 1);
vtilted_lulc.C4 = zeros(3, 1);
vtilted_lulc.index_impervious = 999 * ones(3, 1);
vtilted_lulc.root_depth_m = ones(3, 1);
InputData_Bypass.LULC.table = vtilted_lulc;

InputData_Bypass.Rainfall_Parameters.time_rainfall = [0; 15; 30; 45; 60; 75; 90];
InputData_Bypass.Rainfall_Parameters.intensity_rainfall = 10.8 * ones(7, 1);
InputData_Bypass.Rainfall_Parameters.time_step_rainfall = 15;
InputData_Bypass.Rainfall_Parameters.rainfall_duration = 90;
InputData_Bypass.Rainfall_Parameters.n_obs_rainfall = 7;

routing_mode = lower(strtrim(string(getenv('HYDROPOL2D_VTILT_ROUTING'))));
if strlength(routing_mode) == 0
    routing_mode = "local_inertial";
end

InputData_Bypass.flags.flag_full_momentum = 0;
InputData_Bypass.flags.flag_inertial = 0;
InputData_Bypass.flags.flag_CA = 0;
InputData_Bypass.flags.flag_kinematic = 0;
InputData_Bypass.flags.flag_diffusive = 0;

switch routing_mode
    case {"full_momentum", "fullmomentum", "fm"}
        InputData_Bypass.flags.flag_full_momentum = 1;
    case {"local_inertial", "localinertial", "li"}
        InputData_Bypass.flags.flag_inertial = 1;
    case {"cellular_automata", "ca"}
        InputData_Bypass.flags.flag_CA = 1;
    case {"kinematic", "kinematic_wave", "kw"}
        InputData_Bypass.flags.flag_kinematic = 1;
    case {"diffusive", "diffusive_wave", "dw"}
        InputData_Bypass.flags.flag_diffusive = 1;
    otherwise
        error('Unsupported HYDROPOL2D_VTILT_ROUTING mode: %s', routing_mode);
end
