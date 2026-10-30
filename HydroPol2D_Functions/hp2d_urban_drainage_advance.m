function [S, surface_m3, step] = hp2d_urban_drainage_advance(N, S, surface_m3, dt_s)
%HP2D_URBAN_DRAINAGE_ADVANCE Advance the drainage network over one surface
%   step dt_s [s] in stable sub-steps (_compiled_urban_drainage_advance in
%   HydroPol2D-Python). Each sub-step is min(remaining, stable step); the
%   run stops with an error when the stable step falls below min_dt_s or
%   the surface step needs more than maximum_substeps.
%
%   S: state from hp2d_urban_drainage_init. surface_m3: water volume of each
%   inlet cell (N.inlet_unique_cell order), returned updated.
%   step: outflow_m3 / inflow_m3 through the outfalls in this surface step
%   (node_outflow_m3 / node_inflow_m3 per boundary node), capture_m3 /
%   surcharge_m3 through the inlets, and substeps.

H = N.hydraulics;
remaining = dt_s;
count = 0;
outflow = 0;
inflow = 0;
node_outflow = zeros(numel(N.node_ids), 1);
node_inflow = zeros(numel(N.node_ids), 1);
capture_total = 0;
surcharge_total = 0;
while remaining > 1e-12 && count < H.maximum_substeps
    [stable, ~, label] = hp2d_urban_drainage_stable_dt(N, S.node_head_m, S.pipe_q_m3_s);
    sub_dt = min(remaining, stable);
    if (stable < H.min_dt_s && remaining > stable + 1e-12) || sub_dt <= 0
        error('HydroPol2D:UrbanDrainage', ...
            'Urban-drainage stability requires %.6g s, below the configured minimum %.6g s at %s.', ...
            stable, H.min_dt_s, label);
    end
    previous_q = S.pipe_q_m3_s;
    R = hp2d_urban_drainage_step(N, S.node_volume_m3, previous_q, surface_m3, S.node_head_m, sub_dt);
    surface_m3 = R.surface_m3;
    S.node_volume_m3 = R.node_volume_m3;
    S.node_head_m = R.node_head_m;
    S.pipe_q_m3_s = R.pipe_q_m3_s;
    S.inlet_q_m3_s = R.inlet_q_m3_s;
    captured = max(R.inlet_q_m3_s, 0) * sub_dt;
    surcharged = max(-R.inlet_q_m3_s, 0) * sub_dt;
    S.inlet_cumulative_capture_m3 = S.inlet_cumulative_capture_m3 + captured;
    S.inlet_cumulative_surcharge_m3 = S.inlet_cumulative_surcharge_m3 + surcharged;
    capture_total = capture_total + sum(captured);
    surcharge_total = surcharge_total + sum(surcharged);
    S.node_max_head_m = max(S.node_max_head_m, R.node_head_m);
    above_rim = N.dynamic & R.node_head_m > N.node_rim_m;
    S.node_surcharge_duration_s(above_rim) = S.node_surcharge_duration_s(above_rim) + sub_dt;
    S.pipe_max_abs_q_m3_s = max(S.pipe_max_abs_q_m3_s, abs(R.pipe_q_m3_s));
    S.pipe_max_filling_ratio = max(S.pipe_max_filling_ratio, ...
        hp2d_urban_drainage_filling(N, R.node_head_m, R.pipe_q_m3_s));
    reversed = previous_q .* R.pipe_q_m3_s < 0 & abs(previous_q) > 1e-10 & abs(R.pipe_q_m3_s) > 1e-10;
    S.pipe_reversal_count = S.pipe_reversal_count + reversed;
    S.inlet_peak_abs_q_m3_s = max(S.inlet_peak_abs_q_m3_s, abs(R.inlet_q_m3_s));
    outflow = outflow + R.outfall_outflow_m3;
    inflow = inflow + R.outfall_inflow_m3;
    node_outflow = node_outflow + R.node_outflow_m3;
    node_inflow = node_inflow + R.node_inflow_m3;
    remaining = remaining - sub_dt;
    count = count + 1;
end
if remaining > 1e-12
    error('HydroPol2D:UrbanDrainage', ...
        'Urban drainage exceeded %d substeps during one %g s surface step.', H.maximum_substeps, dt_s);
end
S.cumulative_outfall_outflow_m3 = S.cumulative_outfall_outflow_m3 + outflow;
S.cumulative_outfall_inflow_m3 = S.cumulative_outfall_inflow_m3 + inflow;
S.surface_steps = S.surface_steps + 1;
S.total_substeps = S.total_substeps + count;
S.max_substeps = max(S.max_substeps, count);
step = struct('outflow_m3', outflow, 'inflow_m3', inflow, 'capture_m3', capture_total, ...
    'surcharge_m3', surcharge_total, 'substeps', count, ...
    'node_outflow_m3', node_outflow, 'node_inflow_m3', node_inflow);
end
