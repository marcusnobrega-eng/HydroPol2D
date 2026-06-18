run('/Users/mngomes/Downloads/input_data_bypass_script.m');

scenario = lower(strtrim(string(getenv('HYDROPOL2D_VTILT_40M_SCENARIO'))));
if strlength(scenario) == 0
    scenario = "reference20m";
end

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

InputData_Bypass.general.record_time_maps = 10;
InputData_Bypass.general.record_time_hydrographs = 10;
InputData_Bypass.general.resolution_resample = 40;
InputData_Bypass.subgrid.subgrid_dz_m = 0.005;
InputData_Bypass.subgrid.subgrid_max_depth_m = 2.0;

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
