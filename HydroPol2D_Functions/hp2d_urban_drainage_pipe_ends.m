function [head_start, head_end, depth_start, depth_end] = hp2d_urban_drainage_pipe_ends(N, node_head_m, pipe_q_m3_s)
%HP2D_URBAN_DRAINAGE_PIPE_ENDS Head [m] and depth [m] at both ends of each
%   pipe. The heads are the node heads and the depths are clipped to
%   [0, diameter] above the pipe's end inverts.
%
%   A free outfall end: with free_outfall_depth "critical_normal" (default)
%   it is at min(critical, normal) depth of the pipe's current discharge,
%   the SWMM free-outfall rule (no slope toward the outfall: critical
%   depth). With "invert" it is at the node invert, as in HydroPol2D-Python,
%   so the pipe flows at half its upstream depth (Python review F2).

head = node_head_m(:);
head_start = head(N.pipe_start);
head_end = head(N.pipe_end);
D = N.pipe_diameter_m;
if N.hydraulics.free_outfall_depth == "critical_normal" && (any(N.free_at_start) || any(N.free_at_end))
    free = find(N.free_at_start | N.free_at_end);
    q = abs(pipe_q_m3_s(free));
    y = free_outfall_depth(q, D(free), N.pipe_manning_n(free), N.pipe_free_slope(free));
    at_start = N.free_at_start(free);
    head_start(free(at_start)) = N.pipe_start_invert_m(free(at_start)) + y(at_start);
    head_end(free(~at_start)) = N.pipe_end_invert_m(free(~at_start)) + y(~at_start);
end
depth_start = min(max(head_start - N.pipe_start_invert_m, 0), D);
depth_end = min(max(head_end - N.pipe_end_invert_m, 0), D);
end

function y = free_outfall_depth(q, D, n, slope)
g = 9.80665;
y = zeros(size(q));
flowing = q > 0;
if ~any(flowing)
    return
end
q = q(flowing); D = D(flowing); n = n(flowing); slope = slope(flowing);
% Critical depth: Q^2 T = g A^3 (Froude number 1).
low = zeros(size(q)); high = D;
for iteration = 1:50
    middle = 0.5 * (low + high);
    [A, ~, T] = hp2d_urban_drainage_section(middle, D);
    super = q.^2 .* T > g * A.^3;
    low(super) = middle(super);
    high(~super) = middle(~super);
end
critical = 0.5 * (low + high);
% Normal depth: A R^(2/3) = Q n / sqrt(S) below the peak-conveyance depth.
normal = D;
sloped = slope > 0;
if any(sloped)
    target = q(sloped) .* n(sloped) ./ sqrt(slope(sloped));
    Ds = D(sloped);
    peak = 0.9382 * Ds;
    [A, R] = hp2d_urban_drainage_section(peak, Ds);
    solved = Ds;
    below_peak = A .* R.^(2/3) > target;
    low = zeros(size(Ds)); high = peak;
    for iteration = 1:50
        middle = 0.5 * (low + high);
        [A, R] = hp2d_urban_drainage_section(middle, Ds);
        short = A .* R.^(2/3) < target;
        low(short) = middle(short);
        high(~short) = middle(~short);
    end
    solved(below_peak) = 0.5 * (low(below_peak) + high(below_peak));
    normal(sloped) = solved;
end
y(flowing) = min(critical, normal);
end
