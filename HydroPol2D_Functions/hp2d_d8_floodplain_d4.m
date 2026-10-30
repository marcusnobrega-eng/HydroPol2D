function [volume_m3, qx_m3_s, qy_m3_s, hx_m, hy_m] = hp2d_d8_floodplain_d4( ...
    volume_m3, qx_previous_m3_s, qy_previous_m3_s, stage_m, D8, dt_s, flag_critical)
%HP2D_D8_FLOODPLAIN_D4 D4 local-inertial floodplain of the hybrid D8 scheme.
%   The D4 faces carry the floodplain flow; the faces that a cardinal D8
%   river link already joins are closed (D8.inactive_x / inactive_y), so the
%   channel is not conveyed twice. The water level is the D8 stage (the river
%   cell tables on river cells), the face depth h = max(level) - max(terrain),
%   and the face flow follows Bates et al. (2010) per unit width:
%       q = (q_prev - g h dt S) / (1 + g dt n^2 |q_prev| / h^(7/3))
%   with q = 0 where h <= D8.min_depth_m, n the Manning n of the cell that
%   stores the face, and |q| <= h sqrt(g h) when flag_critical = 1. Donor
%   cells are then volume-limited and the continuity update is conservative.
%   There is no outlet face: water leaves through the D8 outlets. This is
%   local_inertial_step(..., water_surface_m, inactive_internal_x/y,
%   include_normal_flow=False) in HydroPol2D-Python.
%
%   Arrays are ny x nx, as in Local_Inertial_Model_D4: qx(i,j) crosses the
%   east face of cell (i,j), positive eastward; qy(i,j) crosses its north
%   face, positive northward (to row i-1). hx, hy are the face depths [m]
%   (0 on closed faces). volume_m3 is the cell storage [m3].

g = 9.81;
dx = D8.Resolution;
[ny, nx] = deal(D8.ny, D8.nx);
valid = D8.valid_grid;
z = D8.terrain_grid;
n = D8.roughness_grid;
V = reshape(max(volume_m3, 0), ny, nx);
V(~valid) = 0;
ws = reshape(stage_m, ny, nx);

hx_m = zeros(ny, nx, 'like', V);
hy_m = zeros(ny, nx, 'like', V);
Sx = zeros(ny, nx, 'like', V);
Sy = zeros(ny, nx, 'like', V);
active_x = false(ny, nx);
active_y = false(ny, nx);
if nx > 1
    hx_m(:, 1:end-1) = max(max(ws(:, 1:end-1), ws(:, 2:end)) - max(z(:, 1:end-1), z(:, 2:end)), 0);
    Sx(:, 1:end-1) = (ws(:, 2:end) - ws(:, 1:end-1)) / dx;
    active_x(:, 1:end-1) = valid(:, 1:end-1) & valid(:, 2:end);
end
if ny > 1
    hy_m(2:end, :) = max(max(ws(1:end-1, :), ws(2:end, :)) - max(z(1:end-1, :), z(2:end, :)), 0);
    Sy(2:end, :) = (ws(1:end-1, :) - ws(2:end, :)) / dx;
    active_y(2:end, :) = valid(1:end-1, :) & valid(2:end, :);
end
active_x = active_x & ~D8.inactive_x;
active_y = active_y & ~D8.inactive_y;

qx_m3_s = face_discharge(qx_previous_m3_s, hx_m, Sx, n, dx, dt_s, D8.min_depth_m, flag_critical, g);
qy_m3_s = face_discharge(qy_previous_m3_s, hy_m, Sy, n, dx, dt_s, D8.min_depth_m, flag_critical, g);
qx_m3_s(~active_x) = 0;
qy_m3_s(~active_y) = 0;
hx_m(~active_x) = 0;
hy_m(~active_y) = 0;

% Donor limiter: no cell gives more than it holds in one step.
outgoing = max(qx_m3_s, 0) + max(qy_m3_s, 0);
outgoing(:, 2:end) = outgoing(:, 2:end) + max(-qx_m3_s(:, 1:end-1), 0);
outgoing(1:end-1, :) = outgoing(1:end-1, :) + max(-qy_m3_s(2:end, :), 0);
scale = ones(ny, nx, 'like', V);
draining = outgoing > 0 & valid;
scale(draining) = min(1, V(draining) ./ (dt_s * outgoing(draining)));
if nx > 1
    donor = scale(:, 1:end-1);
    upwind = scale(:, 2:end);
    reverse = qx_m3_s(:, 1:end-1) < 0;
    donor(reverse) = upwind(reverse);
    qx_m3_s(:, 1:end-1) = qx_m3_s(:, 1:end-1) .* donor;
end
if ny > 1
    donor = scale(2:end, :);
    upwind = scale(1:end-1, :);
    reverse = qy_m3_s(2:end, :) < 0;
    donor(reverse) = upwind(reverse);
    qy_m3_s(2:end, :) = qy_m3_s(2:end, :) .* donor;
end

net_inflow = -qx_m3_s - qy_m3_s;
net_inflow(:, 2:end) = net_inflow(:, 2:end) + qx_m3_s(:, 1:end-1);
net_inflow(1:end-1, :) = net_inflow(1:end-1, :) + qy_m3_s(2:end, :);
V = max(V + dt_s * net_inflow, 0);
V(~valid) = 0;
volume_m3 = V(:);

end

function q = face_discharge(q_previous, h, slope, n, width, dt_s, min_depth, flag_critical, g)
wet = h > min_depth;
h_wet = h .* wet;
q_unit = q_previous / width;
denominator = 1 + g * dt_s .* n.^2 .* abs(q_unit) ./ max(h_wet, 1e-9).^(7/3);
q_unit = (q_unit - g * h_wet * dt_s .* slope) ./ denominator;
q_unit(~wet) = 0;
if flag_critical == 1
    critical = h_wet .* sqrt(g * h_wet);
    q_unit = min(max(q_unit, -critical), critical);
end
q = q_unit * width;
end
