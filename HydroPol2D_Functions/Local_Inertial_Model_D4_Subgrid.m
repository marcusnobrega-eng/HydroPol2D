function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a,eta_t,V_t] = ...
    Local_Inertial_Model_D4_Subgrid( ...
    flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2, ...
    flag_reservoir,z,d_tot,d_p,roughness_cell,roughness_squared,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet, ...
    row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,nc,nf,River_Width,River_Depth, ...
    Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev,Subgrid_Properties,flag_inflow,SubgridTables)
%LOCAL_INERTIAL_MODEL_D4_SUBGRID SFINCS-style lookup-subgrid local inertial solver.
%
% This production branch follows van Ormondt et al. (2025): continuity is
% updated in wet volume, water levels are recovered from cell volume tables,
% and momentum is solved with grid-average unit discharge q_G using
% velocity-point lookup tables for H_G, n_rep, and wet fraction.

unused_inputs = {flag_numerical_scheme, z, d_p, roughness_cell, roughness_squared, ...
    d_tolerance, flag_critical, nc, nf, River_Width, River_Depth, Qc_prev, ...
    Qf_prev, Qci_prev, Qfi_prev, C_a_prev, Subgrid_Properties, flag_inflow}; %#ok<NASGU>

[ny, nx] = size(d_tot);
dt_total = time_step * 60;
g = 9.81;

if isempty(idx_nan)
    inactive = false(ny, nx);
else
    inactive = logical(idx_nan);
end

zmin_cell = SubgridTables.z_zmin;
eta_n = zmin_cell + max(d_tot ./ 1000, 0);
V_n = hp2d_sfincs_cell_volume_from_zs(SubgridTables, eta_n);
eta_t = eta_n;
V_t = V_n;
C_a = hp2d_sfincs_cell_wet_area_from_zs(SubgridTables, eta_t);

% HydroPol stores previous face memory in equivalent areal flux units for
% compatibility. In the SFINCS branch:
%   1:2 encode q_G [m2/s] after unit conversion;
%   4:5 encode kfuv wet-state memory;
%   6   encodes previous cell water level zs0 [m];
%   7   encodes zsderv [m] for wiggle suppression.
%   8   encodes SFINCS boundary uvmean memory [m2/s].
%   9   encodes SFINCS cell storage volume z_volume [m3].
if isempty(outflow)
    qG_prev = zeros(ny, nx, 2, 'like', eta_n);
    kfuv_prev = false(ny, nx, 2);
    zs_prev = eta_n;
    zsderv_prev = zeros(ny, nx, 'like', eta_n);
    uvmean_boundary_prev = zeros(ny, nx, 'like', eta_n);
    V_t = V_n;
else
    qG_prev = zeros(ny, nx, 2, 'like', eta_n);
    nchan = min(size(outflow, 3), 2);
    qG_prev(:, :, 1:nchan) = outflow(:, :, 1:nchan) ./ 1000 ./ 3600 .* Resolution;
    kfuv_prev = false(ny, nx, 2);
    if size(outflow, 3) >= 5
        kfuv_prev(:, :, 1) = outflow(:, :, 4) > 0.5;
        kfuv_prev(:, :, 2) = outflow(:, :, 5) > 0.5;
    end
    if size(outflow, 3) >= 6 && any(isfinite(outflow(:, :, 6)), 'all') && max(abs(outflow(:, :, 6)), [], 'all') > 0
        zs_prev = outflow(:, :, 6);
    else
        zs_prev = eta_n;
    end
    if size(outflow, 3) >= 7
        zsderv_prev = outflow(:, :, 7);
        zsderv_prev(~isfinite(zsderv_prev)) = 0;
    else
        zsderv_prev = zeros(ny, nx, 'like', eta_n);
    end
    if size(outflow, 3) >= 8
        uvmean_boundary_prev = outflow(:, :, 8);
        uvmean_boundary_prev(~isfinite(uvmean_boundary_prev)) = 0;
    else
        uvmean_boundary_prev = zeros(ny, nx, 'like', eta_n);
    end
    if size(outflow, 3) >= 9
        V_prev_state = outflow(:, :, 9);
        bad_volume = ~isfinite(V_prev_state);
        V_prev_state(bad_volume) = V_n(bad_volume);
        % HydroPol public depth represents max(z_volume,0). Preserve any
        % SFINCS negative storage memory and apply hydrologic source/sink
        % changes inferred from the public depth state.
        V_t = V_prev_state + (V_n - max(V_prev_state, 0));
    else
        V_t = V_n;
    end
end

outflow = zeros(ny, nx, 9, 'like', eta_n);
Hf_x = zeros(ny, nx, 'like', eta_n);
Hf_y = zeros(ny, nx, 'like', eta_n);
phi_x = zeros(ny, nx, 'like', eta_n); %#ok<NASGU>
phi_y = zeros(ny, nx, 'like', eta_n); %#ok<NASGU>
dt = dt_total;
[qG_candidate, kfuv_new, Hf_x, Hf_y, phi_x, phi_y] = ...
    sfincs_flux_update(eta_t, V_t, qG_prev, kfuv_prev, zsderv_prev, ...
    Resolution, dt, g, inactive, SubgridTables);

Q_face = qG_candidate .* Resolution;
Qwest = [zeros(ny, 1, 'like', eta_n), Q_face(:, 1:(nx-1), 1)];
Qeast = Q_face(:, :, 1);
Qnorth = [zeros(1, nx, 'like', eta_n); Q_face(1:(ny-1), :, 2)];
Qsouth = Q_face(:, :, 2);

Vol_Flux = dt .* (Qwest - Qeast + Qnorth - Qsouth);

if flag_reservoir == 1
    Vol_Flux = apply_reservoir_volume_exchange(Vol_Flux, eta_t, zmin_cell, ...
        reservoir_x, reservoir_y, k1, h1, k2, k3, h2, k4, yds1, xds1, ...
        yds2, xds2, dt / 60, cell_area);
end

eta_before_continuity = eta_t;
V_t = V_t + Vol_Flux;
eta_t = hp2d_sfincs_zs_from_cell_volume(SubgridTables, V_t);
zsderv_new = eta_t - 2 .* eta_before_continuity + zs_prev;
C_a = hp2d_sfincs_cell_wet_area_from_zs(SubgridTables, eta_t);

outlet_flow_sub = zeros(ny, nx, 'like', eta_n);
Hf = zeros(ny, nx, 3, 'like', eta_n);
Hf(:, :, 1) = Hf_x;
Hf(:, :, 2) = Hf_y;
if ~isempty(row_outlet)
    [outlet_flow_sub, ~, Hf, V_t, eta_t, C_a, uvmean_boundary_new] = apply_sfincs_outlet( ...
        outlet_flow_sub, outflow, Hf, V_t, eta_t, C_a, SubgridTables, ...
        row_outlet(:), col_outlet(:), outlet_type, slope_outlet, ...
        Resolution, dt, g, ny, nx, uvmean_boundary_prev);
end

qG_prev = qG_candidate;
outflow(:, :, 1:2) = qG_prev(:, :, 1:2) ./ Resolution .* 1000 .* 3600;
outlet_flow = outlet_flow_sub;
outflow(:, :, 3) = outlet_flow;
outflow(:, :, 4) = double(kfuv_new(:, :, 1));
outflow(:, :, 5) = double(kfuv_new(:, :, 2));
outflow(:, :, 6) = eta_before_continuity;
outflow(:, :, 7) = zsderv_new;
if exist('uvmean_boundary_new', 'var')
    outflow(:, :, 8) = uvmean_boundary_new;
else
    outflow(:, :, 8) = uvmean_boundary_prev;
end
outflow(:, :, 9) = V_t;
matrix_store = outflow(:, :, 1:2);

d_t = max(eta_t - zmin_cell, 0) .* 1000;
I_tot_end_cell = abs(sum(outflow(:, :, 1:3), 3)) .* dt_total ./ 1000 ./ 3600 .* Resolution^2;

qout_left = -[zeros(ny, 1, 'like', matrix_store(:, :, 1)), matrix_store(:, 1:end-1, 1)];
qout_right = matrix_store(:, :, 1);
qout_up = matrix_store(:, :, 2);
qout_down = -[matrix_store(2:end, :, 2); zeros(1, nx, 'like', matrix_store(:, :, 2))];

Qc = 0;
Qf = 0;
Qci = 0;
Qfi = 0;
end

function [qG_face, kfuv_new, Hf_x, Hf_y, phi_x, phi_y] = sfincs_flux_update( ...
    eta, V, qG_prev, kfuv_prev, zsderv_prev, dx, dt, g, inactive, S)
[ny, nx] = size(eta);
qG_face = zeros(ny, nx, 2, 'like', eta);
kfuv_new = false(ny, nx, 2);
Hf_x = zeros(ny, nx, 'like', eta);
Hf_y = zeros(ny, nx, 'like', eta);
phi_x = zeros(ny, nx, 'like', eta);
phi_y = zeros(ny, nx, 'like', eta);

if nx > 1
    etaL = eta(:, 1:nx-1);
    etaR = eta(:, 2:nx);
    zu = max(etaL, etaR);
    face = hp2d_sfincs_velocity_state(S, zu, 'u');
    slope = (etaR - etaL) ./ dx;
    active = ~(inactive(:, 1:nx-1) | inactive(:, 2:nx)) & zu > face.zmin & ...
        face.HG > 0 & face.phi > 0 & isfinite(face.gnavg2) & face.gnavg2 > 0 & ...
        isfinite(slope);
    qold = qG_prev(:, 1:nx-1, 1);
    kold = kfuv_prev(:, 1:nx-1, 1);
    qnew = sfincs_lie_step(qold, kold, face.HG, face.gnavg2, slope, dt, g, active);
    qnew = sfincs_wiggle_limiter(qnew, zsderv_prev(:, 1:nx-1), ...
        zsderv_prev(:, 2:nx), active);
    qnew(V(:, 1:nx-1) <= 0) = min(qnew(V(:, 1:nx-1) <= 0), 0);
    qnew(V(:, 2:nx) <= 0) = max(qnew(V(:, 2:nx) <= 0), 0);
    qnew = min(max(qnew, -face.HG .* 10), face.HG .* 10);
    qG_face(:, 1:nx-1, 1) = qnew;
    kfuv_new(:, 1:nx-1, 1) = active;
    Hf_x(:, 1:nx-1) = face.HG;
    phi_x(:, 1:nx-1) = face.phi;
end

if ny > 1
    etaN = eta(1:ny-1, :);
    etaS = eta(2:ny, :);
    zu = max(etaN, etaS);
    face = hp2d_sfincs_velocity_state(S, zu, 'v');
    slope = (etaS - etaN) ./ dx;
    active = ~(inactive(1:ny-1, :) | inactive(2:ny, :)) & zu > face.zmin & ...
        face.HG > 0 & face.phi > 0 & isfinite(face.gnavg2) & face.gnavg2 > 0 & ...
        isfinite(slope);
    qold = qG_prev(1:ny-1, :, 2);
    kold = kfuv_prev(1:ny-1, :, 2);
    qnew = sfincs_lie_step(qold, kold, face.HG, face.gnavg2, slope, dt, g, active);
    qnew = sfincs_wiggle_limiter(qnew, zsderv_prev(1:ny-1, :), ...
        zsderv_prev(2:ny, :), active);
    qnew(V(1:ny-1, :) <= 0) = min(qnew(V(1:ny-1, :) <= 0), 0);
    qnew(V(2:ny, :) <= 0) = max(qnew(V(2:ny, :) <= 0), 0);
    qnew = min(max(qnew, -face.HG .* 10), face.HG .* 10);
    qG_face(1:ny-1, :, 2) = qnew;
    kfuv_new(1:ny-1, :, 2) = active;
    Hf_y(1:ny-1, :) = face.HG;
    phi_y(1:ny-1, :) = face.phi;
end
end

function qnew = sfincs_lie_step(qold, kfuv_prev, HG, gnavg2, slope, dt, g, active)
qnew = zeros(size(qold), 'like', qold);
qfr = abs(qold);
first_wet = active & ~kfuv_prev;
qfr(first_wet) = sqrt(abs(slope(first_wet)) ./ ...
    (max(gnavg2(first_wet), 1.0e-5) ./ 10)) .* HG(first_wet).^(5/3);
den = 1 + dt .* gnavg2 .* qfr ./ max(HG.^(7/3), eps_like(qold));
rhs = qold - g .* dt .* HG .* slope;
qnew(active) = rhs(active) ./ den(active);
qnew(~active) = 0;
qnew(~isfinite(qnew)) = 0;
end

function qnew = sfincs_wiggle_limiter(qnew, zderv_a, zderv_b, active)
wiggle_threshold = 0.1;
wiggle_factor = 0.1;
mdrv = abs(zderv_a - zderv_b) - wiggle_threshold;
mask = active & mdrv > 0;
qnew(mask) = qnew(mask) .* wiggle_threshold ./ ...
    (wiggle_factor .* mdrv(mask) + wiggle_threshold);
end

function Vol_Flux = apply_reservoir_volume_exchange(Vol_Flux, eta, zmin, ...
    reservoir_x, reservoir_y, k1, h1, k2, k3, h2, k4, yds1, xds1, yds2, xds2, ...
    time_step, cell_area)
for ii = 1:length(reservoir_y)
    dtsup = max(eta(reservoir_y(ii), reservoir_x(ii)) - zmin(reservoir_y(ii), reservoir_x(ii)), 0);
    dt_h = time_step / 60;
    if ~isnan(yds1(ii))
        dh = k1(ii) * max(dtsup - h1(ii), 0)^k2(ii) / cell_area * 1000 * 3600 * dt_h;
        Vol_Flux(reservoir_y(ii), reservoir_x(ii)) = ...
            Vol_Flux(reservoir_y(ii), reservoir_x(ii)) - dh / 1000 * cell_area;
        Vol_Flux(yds1(ii), xds1(ii)) = Vol_Flux(yds1(ii), xds1(ii)) + dh / 1000 * cell_area;
    end
    if ~isnan(yds2(ii))
        dh = k3(ii) * max(dtsup - h2(ii), 0)^k4(ii) / cell_area * 1000 * 3600 * dt_h;
        Vol_Flux(reservoir_y(ii), reservoir_x(ii)) = ...
            Vol_Flux(reservoir_y(ii), reservoir_x(ii)) - dh / 1000 * cell_area;
        Vol_Flux(yds2(ii), xds2(ii)) = Vol_Flux(yds2(ii), xds2(ii)) + dh / 1000 * cell_area;
    end
end
end

function [outlet_flow, outflow, Hf, V, eta, C_a, uvmean_boundary] = apply_sfincs_outlet( ...
    outlet_flow, outflow, Hf, V, eta, C_a, S, row_outlet, col_outlet, ...
    outlet_type, slope_outlet, dx, dt, g, ny, nx, uvmean_boundary_prev)
outlet_sub = sub2ind(size(eta), row_outlet, col_outlet);
side_out = outlet_boundary_side(row_outlet, col_outlet, ny, nx);
if isempty(uvmean_boundary_prev)
    uvmean_boundary = zeros(size(eta), 'like', eta);
else
    uvmean_boundary = uvmean_boundary_prev;
    uvmean_boundary(~isfinite(uvmean_boundary)) = 0;
end
Q_out = zeros(size(outlet_sub), 'like', eta);
H_out = zeros(size(outlet_sub), 'like', eta);

btfilter = 60.0;
btrelax = 3600.0;
factime = min(dt / btfilter, 1.0);
one_minus_factime = 1.0 - factime;
facrel = 1.0 - min(dt / btrelax, 1.0);

for ii = 1:numel(outlet_sub)
    side = side_out{ii};
    idx = outlet_sub(ii);
    zsnmi = eta(idx);
    if outlet_type == 1
        if isscalar(slope_outlet)
            Sout = abs(slope_outlet);
        else
            Sout = abs(slope_outlet(idx));
        end
        % SFINCS downstream-river boundary (kcs=5): boundary water level is
        % taken from the inside model and adjusted by the imposed downstream
        % slope over one cell.
        zsnmb = max(zsnmi - Sout * dx, S.z_zmin(idx));
    else
        % SFINCS weakly reflective open boundary with a dry exterior stage.
        zsnmb = S.z_zmin(idx);
    end
    zs0nmb = zsnmb;
    zu = max(zsnmb, zsnmi);
    face = hp2d_sfincs_boundary_velocity_state(S, zu .* ones(size(eta), 'like', eta), side);
    HG = face.HG(idx);
    if HG > 1.0e-6 && isfinite(HG)
        ibuvdir = outlet_ibuvdir(side);
        ui = sqrt(g / HG) * (zsnmb - zs0nmb);
        ub = ibuvdir * (2 * ui - sqrt(g / HG) * (zsnmi - zs0nmb));
        q = ub * HG + uvmean_boundary(idx);

        % SFINCS wet/dry sign constraints at an open boundary. The exterior
        % boundary is dry for HydroPol free-outlet cells, so only outward
        % flux is retained.
        if V(idx) <= 0
            if ibuvdir == 1
                q = max(q, 0);
            else
                q = min(q, 0);
            end
        end
        if zsnmb - S.z_zmin(idx) < S.huthresh
            if ibuvdir == 1
                q = min(q, 0);
            else
                q = max(q, 0);
            end
        end

        qout = max(outlet_outward_sign(side) * q, 0);
        q = outlet_outward_sign(side) * qout;
        Q_out(ii) = qout * dx;
        H_out(ii) = HG;
        uvmean_boundary(idx) = factime * q + facrel * one_minus_factime * uvmean_boundary(idx);
    else
        uvmean_boundary(idx) = 0;
    end
end

Q_out = max(Q_out, 0);
q_equiv = Q_out ./ (dx^2) .* 1000 .* 3600;
outlet_flow(outlet_sub) = q_equiv;
outflow(:, :, 3) = outlet_flow;
V(outlet_sub) = V(outlet_sub) - Q_out .* dt;
eta = hp2d_sfincs_zs_from_cell_volume(S, V);
C_a = hp2d_sfincs_cell_wet_area_from_zs(S, eta);
H3 = Hf(:, :, 3);
H3(outlet_sub) = H_out;
Hf(:, :, 3) = H3;
end

function s = outlet_ibuvdir(side)
switch lower(char(side))
    case {'west', 'north'}
        s = 1;
    case {'east', 'south'}
        s = -1;
    otherwise
        error('Local_Inertial_Model_D4_Subgrid:badOutletSide', ...
            'Unsupported outlet side %s.', side);
end
end

function s = outlet_outward_sign(side)
switch lower(char(side))
    case {'east', 'south'}
        s = 1;
    case {'west', 'north'}
        s = -1;
    otherwise
        error('Local_Inertial_Model_D4_Subgrid:badOutletSide', ...
            'Unsupported outlet side %s.', side);
end
end

function side = outlet_boundary_side(row_outlet, col_outlet, ny, nx)
side = cell(numel(row_outlet), 1);
for ii = 1:numel(row_outlet)
    r = row_outlet(ii);
    c = col_outlet(ii);
    if r <= 1
        side{ii} = 'north';
    elseif r >= ny
        side{ii} = 'south';
    elseif c <= 1
        side{ii} = 'west';
    elseif c >= nx
        side{ii} = 'east';
    else
        dist = [r - 1, ny - r, c - 1, nx - c];
        labels = {'north', 'south', 'west', 'east'};
        [~, idx] = min(dist);
        side{ii} = labels{idx};
        warning('Local_Inertial_Model_D4_Subgrid:interiorOutlet', ...
            ['Outlet cell (%d,%d) is not on a model boundary. ', ...
             'Using nearest boundary side %s.'], r, c, side{ii});
    end
end
end

function e = eps_like(x)
if isa(x, 'gpuArray')
    e = eps(classUnderlying(x));
else
    e = eps(class(x));
end
end
