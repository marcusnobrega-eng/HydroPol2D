% V-tilted 40 m compatibility inputs for retired lookup-subgrid diagnostics.

run(fullfile(model_root, 'Config', 'input_data_bypass_script.m'));

scenario = lower(strtrim(string(getenv('HYDROPOL2D_VTILT_40M_SCENARIO'))));
if strlength(scenario) == 0
    scenario = "reference20m";
end

InputData_Bypass.general.date_begin = datetime(2025, 5, 1, 0, 0, 0);
InputData_Bypass.general.date_end = InputData_Bypass.general.date_begin + minutes(240);
InputData_Bypass.general.routing_time = 240;
InputData_Bypass.general.record_time_maps = 10;
InputData_Bypass.general.record_time_hydrographs = 10;
InputData_Bypass.general.resolution_resample = 40;
InputData_Bypass.general.slope_outlet = 0.02;
if ~isfield(InputData_Bypass, 'subgrid')
    InputData_Bypass.subgrid = struct();
end
InputData_Bypass.subgrid.subgrid_dz_m = 0.005;
InputData_Bypass.subgrid.subgrid_max_depth_m = 2.0;

InputData_Bypass.Rainfall_Parameters.time_rainfall = [0; 15; 30; 45; 60; 75; 90];
InputData_Bypass.Rainfall_Parameters.intensity_rainfall = 10.8 * ones(7, 1);
InputData_Bypass.Rainfall_Parameters.time_step_rainfall = 15;
InputData_Bypass.Rainfall_Parameters.rainfall_duration = 90;
InputData_Bypass.Rainfall_Parameters.n_obs_rainfall = 7;

InputData_Bypass.flags.flag_full_momentum = 0;
InputData_Bypass.flags.flag_inertial = 1;
InputData_Bypass.flags.flag_CA = 0;
InputData_Bypass.flags.flag_kinematic = 0;
InputData_Bypass.flags.flag_diffusive = 0;
InputData_Bypass.flags.flag_resample = 0;
InputData_Bypass.flags.flag_spatial_rainfall = 0;
InputData_Bypass.flags.flag_input_rainfall_map = 0;
InputData_Bypass.flags.flag_satellite_rainfall = 0;
InputData_Bypass.flags.flag_real_time_satellite_rainfall = 0;
InputData_Bypass.flags.flag_infiltration = 0;
InputData_Bypass.flags.flag_ETP = 0;
InputData_Bypass.flags.flag_groundwater_modeling = 0;
InputData_Bypass.flags.flag_export_maps = 1;

switch scenario
    case {"reference20m", "baseline40m"}
        InputData_Bypass.flags.flag_subgrid = 0;
        InputData_Bypass.flags.flag_overbanks = 0;
        InputData_Bypass.flags.flag_river_rasters = 0;
    case "subgrid40m"
        InputData_Bypass.flags.flag_subgrid = 1;
        InputData_Bypass.flags.flag_overbanks = 0;
        InputData_Bypass.flags.flag_river_rasters = 0;
    otherwise
        error('Unsupported HYDROPOL2D_VTILT_40M_SCENARIO: %s', scenario);
end
