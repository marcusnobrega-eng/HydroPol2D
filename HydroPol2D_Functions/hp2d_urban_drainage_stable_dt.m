function [dt_s, limiting, label] = hp2d_urban_drainage_stable_dt(N, node_head_m, pipe_q_m3_s)
%HP2D_URBAN_DRAINAGE_STABLE_DT Largest stable drainage sub-step [s]
%   (_pipe_stable_dt_s in HydroPol2D-Python): the rigid-column limit of
%   each wet pipe, dt = C sqrt(L S / (g A)), with A the wetted area at the
%   mean end depth and S the smaller storage area of its dynamic end nodes,
%   capped at max_dt_s. limiting is the pipe that sets the step; label names
%   it for error messages.

g = 9.80665;
tiny = 1e-12;
H = N.hydraulics;
[~, ~, depth_start, depth_end] = hp2d_urban_drainage_pipe_ends(N, node_head_m, pipe_q_m3_s);
area = hp2d_urban_drainage_section(0.5 * (depth_start + depth_end), N.pipe_diameter_m);
[~, storage_area] = hp2d_urban_drainage_storage(N, node_head_m);
s_start = storage_area(N.pipe_start);
s_end = storage_area(N.pipe_end);
dyn_start = N.dynamic(N.pipe_start);
dyn_end = N.dynamic(N.pipe_end);
effective = s_end;
effective(dyn_start) = s_start(dyn_start);
both = dyn_start & dyn_end;
effective(both) = min(s_start(both), s_end(both));
wet = area > tiny;
limits = H.max_dt_s * ones(numel(area), 1);
limits(wet) = H.courant * sqrt(N.pipe_length_m(wet) .* max(effective(wet), tiny) ...
    ./ (g * max(area(wet), tiny)));
[dt_s, limiting] = min(limits);
dt_s = min(dt_s, H.max_dt_s);
label = sprintf('pipe ''%s'' between nodes ''%s'' and ''%s''', N.pipe_ids(limiting), ...
    N.node_ids(N.pipe_start(limiting)), N.node_ids(N.pipe_end(limiting)));
end
