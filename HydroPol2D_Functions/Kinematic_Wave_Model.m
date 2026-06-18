function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = Kinematic_Wave_Model(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2,flag_reservoir,z,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width,River_Depth,Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev)
%KINEMATIC_WAVE_MODEL Conservative D4 kinematic-wave routing.
% The kinematic closure computes Manning face fluxes from bed slope only,
% using the upstream cell depth on each face. This keeps the formulation
% distinct from the diffusive wave, which responds to water-surface slope.

[qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
    Routing_Wave_Model_D4('kinematic',flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2,flag_reservoir,z,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width,River_Depth,Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev);
end
