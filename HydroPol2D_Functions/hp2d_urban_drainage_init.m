function S = hp2d_urban_drainage_init(N)
%HP2D_URBAN_DRAINAGE_INIT Drainage state of a fresh run: node volumes at
%   the initial heads, still pipes and zeroed counters
%   (initialize_urban_drainage_state in HydroPol2D-Python).

n_nodes = numel(N.node_ids);
n_pipes = numel(N.pipe_ids);
n_inlets = numel(N.inlet_ids);
S = struct();
S.node_volume_m3 = hp2d_urban_drainage_storage(N, N.node_initial_head_m);
S.node_head_m = hp2d_urban_drainage_head(N, S.node_volume_m3);
S.pipe_q_m3_s = zeros(n_pipes, 1);
S.inlet_q_m3_s = zeros(n_inlets, 1);
S.inlet_cumulative_capture_m3 = zeros(n_inlets, 1);
S.inlet_cumulative_surcharge_m3 = zeros(n_inlets, 1);
S.node_max_head_m = S.node_head_m;
S.node_surcharge_duration_s = zeros(n_nodes, 1);
S.pipe_max_abs_q_m3_s = zeros(n_pipes, 1);
S.pipe_max_filling_ratio = zeros(n_pipes, 1);
S.pipe_reversal_count = zeros(n_pipes, 1);
S.inlet_peak_abs_q_m3_s = zeros(n_inlets, 1);
S.cumulative_outfall_outflow_m3 = 0;
S.cumulative_outfall_inflow_m3 = 0;
S.surface_steps = 0;
S.total_substeps = 0;
S.max_substeps = 0;
S.initial_storage_m3 = sum(S.node_volume_m3);
end
