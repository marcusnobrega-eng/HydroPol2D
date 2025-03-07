% Delta storage
cell = [71,20];


i_a_cell = Hydro_States.i_a(cell(1),cell(2))

ksat_cell = Soil_Properties.ksat(cell(1),cell(2));
teta_sat_cell = Soil_Properties.teta_sat(cell(1),cell(2))
psi_cell = Soil_Properties.psi(cell(1),cell(2))
I_t_cell = Soil_Properties.I_t(cell(1),cell(2))
I_p_cell = Soil_Properties.I_p(cell(1),cell(2))

n_cell = LULC_Properties.roughness(cell(1),cell(2))
depth_cell_p = depths.d_p(cell(1),cell(2))
h_GW_cell_t = BC_States.h_t(cell(1),cell(2)) - elevation(cell(1),cell(2)) + Soil_Properties.Soil_Depth(cell(1),cell(2))
h_GW_cell_p = BC_States.h_t(cell(1),cell(2)) - elevation(cell(1),cell(2)) + Soil_Properties.Soil_Depth(cell(1),cell(2))
ETR_cell = Hydro_States.ETR(cell(1),cell(2))
E_cell = BC_States.delta_E(cell(1),cell(2))
f_cell = Hydro_States.f(cell(1),cell(2))
q_river_cell = q_river(cell(1),cell(2))
q_exf_cell =  q_exf(cell(1),cell(2)) 
depths_end_cell = depths.d_p(cell(1),cell(2)) - inf_volume_cell(cell(1),cell(2))
inf_vol_cell = inf_volume_cell(cell(1),cell(2))
delta_p_cell = BC_States.delta_p_agg(cell(1),cell(2))/(time_step/60) % mm/h
slope_left = (elevation(cell(1),cell(2)-1) - elevation(cell(1),cell(2)))/Wshed_Properties.Resolution
slope_right = (elevation(cell(1),cell(2)+1) - elevation(cell(1),cell(2)))/Wshed_Properties.Resolution
slope_up  = (elevation(cell(1)+1,cell(2)) - elevation(cell(1),cell(2)))/Wshed_Properties.Resolution
slope_down = (elevation(cell(1)-1,cell(2)) - elevation(cell(1),cell(2)))/Wshed_Properties.Resolution
flow_right = outflow_bates(cell(1),cell(2),1)
flow_up = outflow_bates(cell(1),cell(2),2)