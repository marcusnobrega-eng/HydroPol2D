% HydroPol2D Solver

function [Qmod, Cmod,Qmod_gauges,Cmod_gauges] = HydroPol2D_Routing_Solver(t,time_obs,easting_obs_cell,northing_obs_cell,min_soil_moisture,Reservoir_Data,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties)


z2_save = 0;

    
    % Main While
    flags.flag_automatic_calibration = 1; % Activating saving modeled data in HydroPol2D_Main_While
%     flags.flag_obs_gauges = 0;
    % Run HydroPol2D Routing Solver

    HydroPol2D_Main_While

end