% Run the linear reservoir model
GW_Depth = (BC_States.h_t - (elevation - Soil_Properties.Soil_Depth));
UZ_max_storage = (Soil_Properties.Soil_Depth - GW_Depth).*(Soil_Properties.teta_sat - Soil_Properties.teta_i); % meters of water
[recharge_rate, Soil_Moisture] = simulate_groundwater_recharge(0*Hydro_States.f/1000/3600, Soil_Properties.I_t/1000, Soil_Properties.k, time_step*60, min_soil_moisture/1000, UZ_max_storage, LULC_Properties.idx_imp);
Soil_Properties.I_t = Soil_Moisture*1000;

% Groundwater_Module
if flags.flag_baseflow == 1
    % Run the 2D Boussinesq Model
    [BC_States.h_t, ~, ~, q_exf, q_river] = Boussinesq_2D_explicit(time_step*60,Wshed_Properties.Resolution,Wshed_Properties.Resolution,BC_States.h_t,(Elevation_Properties.elevation_cell - Soil_Properties.Soil_Depth), ...
        Soil_Properties.Sy,recharge_rate,Soil_Properties.ksat_gw/1000/3600,idx_rivers,LULC_Properties.River_K_coeff*Soil_Properties.ksat/1000/3600, ...
        (Elevation_Properties.elevation_cell + depths.d_t/1000),Elevation_Properties.elevation_cell,0.5,Soil_Properties.Soil_Depth,Wshed_Properties.domain,[],[],Wshed_Properties.perimeter);
    BC_States.h_0 = BC_States.h_t; % Refreshing previous state
    % Considering excess of saturation and river exchanges
    depths.d_t = depths.d_t + time_step*60*(q_river*1000 + q_exf*1000); %  mm considering river iteractions and excess of saturation
end
