function files = hp2d_urban_drainage_write_outputs(Drainage_States, folder)
%HP2D_URBAN_DRAINAGE_WRITE_OUTPUTS Write the drainage records and the
%   end-of-run summaries as CSV files in folder (the drainage variables of
%   HydroPol2D-Python's results file):
%     Urban_Drainage_Totals.csv  storage and cumulative outfall, capture and
%                                surcharge volumes at each record
%     Urban_Drainage_Nodes.csv   head, depth, storage and flux components
%     Urban_Drainage_Pipes.csv   discharge, velocity and filling
%     Urban_Drainage_Inlets.csv  exchange (positive into the network)
%     Urban_Drainage_Node_Summary.csv, _Pipe_Summary.csv, _Inlet_Summary.csv

if ~exist(folder, 'dir')
    mkdir(folder);
end
N = Drainage_States.Network;
S = Drainage_States.State;
Hs = Drainage_States.History;
files = {};
t = Hs.time_min(:);
nt = numel(t);

Totals = table(t, Hs.total_storage_m3(:), Hs.cumulative_outfall_outflow_m3(:), ...
    Hs.cumulative_outfall_inflow_m3(:), Hs.cumulative_capture_m3(:), Hs.cumulative_surcharge_m3(:), ...
    'VariableNames', {'time_min', 'storage_m3', 'cumulative_outfall_outflow_m3', ...
    'cumulative_outfall_inflow_m3', 'cumulative_capture_m3', 'cumulative_surcharge_m3'});
files{end + 1} = write(Totals, folder, 'Urban_Drainage_Totals.csv');

if nt > 0
    nn = numel(N.node_ids);
    node_id = repmat(N.node_ids(:).', nt, 1);
    time_n = repmat(t, 1, nn);
    depth = Hs.node_head_m - N.node_invert_m.';
    Nodes = table(time_n(:), node_id(:), Hs.node_head_m(:), depth(:), Hs.node_storage_m3(:), ...
        Hs.node_inflow_m3_s(:), Hs.node_outflow_m3_s(:), Hs.node_surface_capture_m3_s(:), ...
        Hs.node_surface_surcharge_m3_s(:), Hs.node_outfall_inflow_m3_s(:), Hs.node_outfall_outflow_m3_s(:), ...
        'VariableNames', {'time_min', 'node_id', 'head_m', 'depth_m', 'storage_m3', 'inflow_m3_s', ...
        'outflow_m3_s', 'surface_capture_m3_s', 'surface_surcharge_m3_s', ...
        'outfall_inflow_m3_s', 'outfall_outflow_m3_s'});
    files{end + 1} = write(Nodes, folder, 'Urban_Drainage_Nodes.csv');
    np = numel(N.pipe_ids);
    pipe_id = repmat(N.pipe_ids(:).', nt, 1);
    time_p = repmat(t, 1, np);
    Pipes = table(time_p(:), pipe_id(:), Hs.pipe_discharge_m3_s(:), Hs.pipe_velocity_m_s(:), ...
        Hs.pipe_filling_ratio(:), 'VariableNames', {'time_min', 'pipe_id', 'discharge_m3_s', ...
        'velocity_m_s', 'filling_ratio'});
    files{end + 1} = write(Pipes, folder, 'Urban_Drainage_Pipes.csv');
    ni = numel(N.inlet_ids);
    if ni > 0
        inlet_id = repmat(N.inlet_ids(:).', nt, 1);
        time_i = repmat(t, 1, ni);
        Inlets = table(time_i(:), inlet_id(:), Hs.inlet_discharge_m3_s(:), ...
            Hs.inlet_cumulative_capture_m3(:), Hs.inlet_cumulative_surcharge_m3(:), ...
            'VariableNames', {'time_min', 'inlet_id', 'discharge_m3_s', ...
            'cumulative_capture_m3', 'cumulative_surcharge_m3'});
        files{end + 1} = write(Inlets, folder, 'Urban_Drainage_Inlets.csv');
    end
end

NodeSummary = table(N.node_ids(:), N.node_invert_m, N.node_rim_m, S.node_head_m, ...
    S.node_head_m - N.node_invert_m, S.node_max_head_m, max(S.node_max_head_m - N.node_rim_m, 0), ...
    S.node_surcharge_duration_s, 'VariableNames', {'node_id', 'invert_elevation_m', ...
    'ground_elevation_m', 'final_head_m', 'final_depth_m', 'maximum_head_m', ...
    'maximum_surcharge_m', 'surcharge_duration_s'});
files{end + 1} = write(NodeSummary, folder, 'Urban_Drainage_Node_Summary.csv');
PipeSummary = table(N.pipe_ids(:), N.pipe_length_m, N.pipe_diameter_m, N.pipe_manning_n, ...
    N.pipe_start_invert_m, N.pipe_end_invert_m, S.pipe_q_m3_s, S.pipe_max_abs_q_m3_s, ...
    S.pipe_max_filling_ratio, S.pipe_reversal_count, 'VariableNames', {'pipe_id', 'length_m', ...
    'diameter_m', 'manning_n', 'upstream_invert_m', 'downstream_invert_m', ...
    'final_discharge_m3_s', 'maximum_discharge_m3_s', 'maximum_filling_ratio', 'reversal_count'});
files{end + 1} = write(PipeSummary, folder, 'Urban_Drainage_Pipe_Summary.csv');
if ~isempty(N.inlet_ids)
    InletSummary = table(N.inlet_ids(:), N.node_ids(N.inlet_node), S.inlet_peak_abs_q_m3_s, ...
        S.inlet_cumulative_capture_m3, S.inlet_cumulative_surcharge_m3, ...
        'VariableNames', {'inlet_id', 'node_id', 'peak_discharge_m3_s', ...
        'cumulative_capture_m3', 'cumulative_surcharge_m3'});
    files{end + 1} = write(InletSummary, folder, 'Urban_Drainage_Inlet_Summary.csv');
end
end

function path = write(T, folder, name)
path = fullfile(folder, name);
writetable(T, path);
end
