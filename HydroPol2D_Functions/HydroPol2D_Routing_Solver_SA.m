% HydroPol2D Solver

function [Qmod, Cmod,Dmod,Flooded_Area,Risk_Area] = HydroPol2D_Routing_Solver_SA(min_soil_moisture,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties,idx_rivers,Lateral_Groundwater_Flux,Reservoir_Data,Input_Rainfall)


z2_save = 0;
Qmod = nan;
Cmod = nan;
Dmod = nan;
Flooded_Area = nan;
Risk_Area = nan;     
    % Main While
%     flags.flag_automatic_calibration = 1; % Activating saving modeled data in HydroPol2D_Main_While
    flags.flag_obs_gauges = 0; 
    
    % Run HydroPol2D Routing Solver

    HydroPol2D_Main_While

    % Save Data
    Qmod = outlet_states.outlet_hydrograph; 
    Dmod = outlet_states.depth_outlet; 
    if flags.flag_waterquality == 1
        Cmod = WQ_States.outet_pollutograph; 
    end

end