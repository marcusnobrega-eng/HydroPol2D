function R = hp2d_layered_lateral_advance(column, head_m, water_mm, duration_s, dx_m, dy_m, dirichlet_mask, dirichlet_head_m, courant, max_substeps)
%HP2D_LAYERED_LATERAL_ADVANCE D4 Dupuit flow with layered transmissivity on a grid.
%   head_m is [ny x nx]; column slot fields and water_mm are [ny*nx x L] in MATLAB linear
%   order. Explicit sub-steps limited by courant*min(dx,dy)^2*Sy/(4T), with Sy the
%   drainable porosity at the table (floor 1e-3). Face T = mean(T_i, T_j). Outflow is
%   limited by each donor's drainable water. Net volumes go through
%   hp2d_water_table_update. Dirichlet cells keep their head; their exchange is
%   boundary_inflow_m3. Same scheme as HydroPol2D-Python layered_lateral_advance.
if nargin < 7 || isempty(dirichlet_mask), dirichlet_mask = false(size(head_m)); end
if nargin < 8 || isempty(dirichlet_head_m), dirichlet_head_m = head_m; end
if nargin < 9 || isempty(courant), courant = 0.25; end
if nargin < 10 || isempty(max_substeps), max_substeps = 100000; end
[ny, nx] = size(head_m);
valid = reshape(column.valid, ny, nx);
fixed = logical(dirichlet_mask) & valid;
free = valid & ~fixed;
cell_area = dx_m*dy_m;
min_spacing = min(dx_m, dy_m);
head = head_m;
head(fixed) = dirichlet_head_m(fixed);
water = water_mm;
exfil = zeros(ny, nx, 'like', head_m);
unmet = exfil;
boundary = 0;
t = 0;
n = 0;
while t < duration_s*(1 - 1e-12) && n < max_substeps
    T = reshape(hp2d_layered_transmissivity(column, head(:)), ny, nx);
    sy = reshape(max(hp2d_effective_specific_yield(column, head(:)), 1e-3), ny, nx);
    ratio = inf(ny, nx, 'like', head_m);
    ok = valid & T > 0;
    ratio(ok) = sy(ok)./max(T(ok), 1e-300);
    dt = min(courant*min_spacing^2*min(ratio(:))/4, duration_s - t);
    both_x = valid(:, 1:end-1) & valid(:, 2:end);
    both_y = valid(1:end-1, :) & valid(2:end, :);
    tx = 0.5*(T(:, 1:end-1) + T(:, 2:end));
    ty = 0.5*(T(1:end-1, :) + T(2:end, :));
    qx = -tx.*(head(:, 2:end) - head(:, 1:end-1))/dx_m*dy_m;
    qy = -ty.*(head(2:end, :) - head(1:end-1, :))/dy_m*dx_m;
    qx(~both_x) = 0;
    qy(~both_y) = 0;
    out_rate = zeros(ny, nx, 'like', head_m);
    out_rate(:, 1:end-1) = out_rate(:, 1:end-1) + max(qx, 0);
    out_rate(:, 2:end) = out_rate(:, 2:end) + max(-qx, 0);
    out_rate(1:end-1, :) = out_rate(1:end-1, :) + max(qy, 0);
    out_rate(2:end, :) = out_rate(2:end, :) + max(-qy, 0);
    available = reshape(hp2d_drainable_storage_mm(column, head(:)), ny, nx)/1000*cell_area;
    factor = ones(ny, nx, 'like', head_m);
    limited = free & out_rate*dt > available;
    factor(limited) = available(limited)./max(out_rate(limited)*dt, 1e-300);
    fx_left = factor(:, 1:end-1); fx_right = factor(:, 2:end);
    qx(qx > 0) = qx(qx > 0).*fx_left(qx > 0);
    qx(qx < 0) = qx(qx < 0).*fx_right(qx < 0);
    fy_up = factor(1:end-1, :); fy_down = factor(2:end, :);
    qy(qy > 0) = qy(qy > 0).*fy_up(qy > 0);
    qy(qy < 0) = qy(qy < 0).*fy_down(qy < 0);
    net = zeros(ny, nx, 'like', head_m);
    net(:, 1:end-1) = net(:, 1:end-1) - qx;
    net(:, 2:end) = net(:, 2:end) + qx;
    net(1:end-1, :) = net(1:end-1, :) - qy;
    net(2:end, :) = net(2:end, :) + qy;
    net_m3 = net*dt;
    boundary = boundary - sum(net_m3(fixed));
    net_mm = net_m3/cell_area*1000;
    net_mm(~free) = 0;
    [new_head, water, exf, unm] = hp2d_water_table_update(column, head(:), water, net_mm(:));
    head = reshape(new_head, ny, nx);
    head(fixed) = dirichlet_head_m(fixed);
    exfil = exfil + reshape(exf, ny, nx);
    unmet = unmet + reshape(unm, ny, nx);
    t = t + dt;
    n = n + 1;
end
R.head_m = head;
R.water_mm = water;
R.exfiltration_mm = exfil;
R.unmet_removal_mm = unmet;
R.boundary_inflow_m3 = boundary;
R.substeps = n;
end
