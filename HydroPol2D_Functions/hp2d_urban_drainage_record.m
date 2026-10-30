function Drainage_States = hp2d_urban_drainage_record(Drainage_States, index, time_min)
%HP2D_URBAN_DRAINAGE_RECORD Store the network state at a hydrograph record
%   (the drainage history HydroPol2D-Python keeps at its diagnostic times).

N = Drainage_States.Network;
S = Drainage_States.State;
Hs = Drainage_States.History;
row = numel(Hs.index) + 1;
if ~isempty(Hs.index) && Hs.index(end) == index
    row = numel(Hs.index);
end
[filling, velocity] = hp2d_urban_drainage_filling(N, S.node_head_m, S.pipe_q_m3_s);
F = hp2d_urban_drainage_node_fluxes(N, S.pipe_q_m3_s, S.inlet_q_m3_s);
Hs.index(row, 1) = index;
Hs.time_min(row, 1) = double(gather(time_min));
Hs.total_storage_m3(row, 1) = sum(S.node_volume_m3);
Hs.cumulative_outfall_outflow_m3(row, 1) = S.cumulative_outfall_outflow_m3;
Hs.cumulative_outfall_inflow_m3(row, 1) = S.cumulative_outfall_inflow_m3;
Hs.cumulative_capture_m3(row, 1) = sum(S.inlet_cumulative_capture_m3);
Hs.cumulative_surcharge_m3(row, 1) = sum(S.inlet_cumulative_surcharge_m3);
Hs.node_head_m(row, :) = S.node_head_m.';
Hs.node_storage_m3(row, :) = S.node_volume_m3.';
Hs.node_inflow_m3_s(row, :) = F.total_inflow_m3_s.';
Hs.node_outflow_m3_s(row, :) = F.total_outflow_m3_s.';
Hs.node_surface_capture_m3_s(row, :) = F.surface_capture_m3_s.';
Hs.node_surface_surcharge_m3_s(row, :) = F.surface_surcharge_m3_s.';
Hs.node_outfall_inflow_m3_s(row, :) = F.outfall_inflow_m3_s.';
Hs.node_outfall_outflow_m3_s(row, :) = F.outfall_outflow_m3_s.';
Hs.pipe_discharge_m3_s(row, :) = S.pipe_q_m3_s.';
Hs.pipe_velocity_m_s(row, :) = velocity.';
Hs.pipe_filling_ratio(row, :) = filling.';
Hs.inlet_discharge_m3_s(row, :) = S.inlet_q_m3_s.';
Hs.inlet_cumulative_capture_m3(row, :) = S.inlet_cumulative_capture_m3.';
Hs.inlet_cumulative_surcharge_m3(row, :) = S.inlet_cumulative_surcharge_m3.';
Drainage_States.History = Hs;
end
