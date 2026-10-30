function [filling, velocity_m_s] = hp2d_urban_drainage_filling(N, node_head_m, pipe_q_m3_s)
%HP2D_URBAN_DRAINAGE_FILLING Pipe filling ratio (mean end depth over the
%   diameter) and signed velocity [m/s] for reporting (pipe_filling_ratio
%   and pipe_hydraulic_state in HydroPol2D-Python). The velocity is zero in
%   a dry pipe and clipped to the maximum pipe velocity.

[~, ~, depth_start, depth_end] = hp2d_urban_drainage_pipe_ends(N, node_head_m, pipe_q_m3_s);
depth = 0.5 * (depth_start + depth_end);
filling = depth ./ N.pipe_diameter_m;
if nargout > 1
    area = hp2d_urban_drainage_section(depth, N.pipe_diameter_m);
    velocity_m_s = zeros(size(depth));
    wet = depth > 1e-6 & area > 1e-12;
    velocity_m_s(wet) = pipe_q_m3_s(wet) ./ area(wet);
    vmax = N.hydraulics.maximum_velocity_m_s;
    velocity_m_s = min(max(velocity_m_s, -vmax), vmax);
end
end
