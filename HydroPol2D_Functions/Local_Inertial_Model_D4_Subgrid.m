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
% compatibility. In the SFINCS branch, channels 1 and 2 encode q_G.
if isempty(outflow)
    qG_prev = zeros(ny, nx, 2, 'like', eta_n);
else
    qG_prev = zeros(ny, nx, 2, 'like', eta_n);
    nchan = min(size(outflow, 3), 2);
    qG_prev(:, :, 1:nchan) = outflow(:, :, 1:nchan) ./ 1000 ./ 3600 .* Resolution;
end

outflow = zeros(ny, nx, 3, 'like', eta_n);
Hf_x = zeros(ny, nx, 'like', eta_n);
Hf_y = zeros(ny, nx, 'like', eta_n);
phi_x = zeros(ny, nx, 'like', eta_n); %#ok<NASGU>
phi_y = zeros(ny, nx, 'like', eta_n); %#ok<NASGU>
outlet_volume_total = zeros(ny, nx, 'like', eta_n);

dt_remaining = dt_total;
dt_try = dt_total;
min_sub_dt = max(dt_total / 64, 1e-4);
negative_tol = 1e-9;

while dt_remaining > 1e-12
    dt_sub = min(dt_try, dt_remaining);
    [qG_candidate, Hf_x_candidate, Hf_y_candidate, phi_x_candidate, phi_y_candidate] = ...
        sfincs_flux_update(eta_t, qG_prev, Resolution, dt_sub, g, inactive, SubgridTables);

    Q_face = qG_candidate .* Resolution;
    Q_face = limit_face_discharge_by_available_volume(Q_face, V_t, dt_sub);
    qG_candidate = Q_face ./ Resolution;
    Qwest = [zeros(ny, 1, 'like', eta_n), Q_face(:, 1:(nx-1), 1)];
    Qeast = Q_face(:, :, 1);
    Qnorth = [zeros(1, nx, 'like', eta_n); Q_face(1:(ny-1), :, 2)];
    Qsouth = Q_face(:, :, 2);

    Vol_Flux = dt_sub .* (Qwest - Qeast + Qnorth - Qsouth);

    if flag_reservoir == 1
        Vol_Flux = apply_reservoir_volume_exchange(Vol_Flux, eta_t, zmin_cell, ...
            reservoir_x, reservoir_y, k1, h1, k2, k3, h2, k4, yds1, xds1, ...
            yds2, xds2, dt_sub / 60, cell_area);
    end

    V_candidate = V_t + Vol_Flux;
    if any(V_candidate(:) < -negative_tol) && dt_sub > min_sub_dt
        dt_try = max(0.5 * dt_sub, min_sub_dt);
        continue;
    end

    V_t = max(V_candidate, 0);
    eta_t = hp2d_sfincs_zs_from_cell_volume(SubgridTables, V_t);
    C_a = hp2d_sfincs_cell_wet_area_from_zs(SubgridTables, eta_t);

    outlet_flow_sub = zeros(ny, nx, 'like', eta_n);
    Hf_sub = zeros(ny, nx, 3, 'like', eta_n);
    Hf_sub(:, :, 1) = Hf_x_candidate;
    Hf_sub(:, :, 2) = Hf_y_candidate;
    if ~isempty(row_outlet)
        [outlet_flow_sub, ~, Hf_sub, V_t, eta_t, C_a] = apply_sfincs_outlet( ...
            outlet_flow_sub, outflow, Hf_sub, V_t, eta_t, C_a, SubgridTables, ...
            row_outlet(:), col_outlet(:), outlet_type, slope_outlet, ...
            Resolution, dt_sub, g, ny, nx);
    end

    outlet_volume_total = outlet_volume_total + ...
        outlet_flow_sub ./ 1000 ./ 3600 .* Resolution^2 .* dt_sub;

    qG_prev = qG_candidate;
    Hf_x = Hf_x_candidate;
    Hf_y = Hf_y_candidate;
    phi_x = phi_x_candidate; %#ok<NASGU>
    phi_y = phi_y_candidate; %#ok<NASGU>
    Hf = Hf_sub;

    dt_remaining = dt_remaining - dt_sub;
    dt_try = min(dt_try * 1.25, max(dt_remaining, min_sub_dt));
end

outflow(:, :, 1:2) = qG_prev(:, :, 1:2) ./ Resolution .* 1000 .* 3600;
outlet_flow = outlet_volume_total ./ Resolution^2 .* 1000 .* 3600 ./ dt_total;
outflow(:, :, 3) = outlet_flow;
matrix_store = outflow(:, :, 1:2);

d_t = max(eta_t - zmin_cell, 0) .* 1000;
I_tot_end_cell = abs(sum(outflow, 3)) .* dt_total ./ 1000 ./ 3600 .* Resolution^2;

qout_left = -[zeros(ny, 1, 'like', matrix_store(:, :, 1)), matrix_store(:, 1:end-1, 1)];
qout_right = matrix_store(:, :, 1);
qout_up = matrix_store(:, :, 2);
qout_down = -[matrix_store(2:end, :, 2); zeros(1, nx, 'like', matrix_store(:, :, 2))];

Qc = 0;
Qf = 0;
Qci = 0;
Qfi = 0;
end

function Q_face = limit_face_discharge_by_available_volume(Q_face, V, dt)
[ny, nx, ~] = size(Q_face);
Qwest = [zeros(ny, 1, 'like', V), Q_face(:, 1:(nx-1), 1)];
Qeast = Q_face(:, :, 1);
Qnorth = [zeros(1, nx, 'like', V); Q_face(1:(ny-1), :, 2)];
Qsouth = Q_face(:, :, 2);

out_rate = max(Qeast, 0) + max(-Qwest, 0) + max(Qsouth, 0) + max(-Qnorth, 0);
scale = ones(ny, nx, 'like', V);
needs_limit = out_rate .* dt > max(V, 0) & out_rate > 0;
scale(needs_limit) = max(V(needs_limit), 0) ./ (out_rate(needs_limit) .* dt);
scale = max(min(scale, 1), 0);

if nx > 1
    Qx = Q_face(:, 1:(nx-1), 1);
    donor_scale = ones(size(Qx), 'like', Qx);
    pos = Qx >= 0;
    left_scale = scale(:, 1:(nx-1));
    right_scale = scale(:, 2:nx);
    donor_scale(pos) = left_scale(pos);
    donor_scale(~pos) = right_scale(~pos);
    Q_face(:, 1:(nx-1), 1) = Qx .* donor_scale;
end

if ny > 1
    Qy = Q_face(1:(ny-1), :, 2);
    donor_scale = ones(size(Qy), 'like', Qy);
    pos = Qy >= 0;
    north_scale = scale(1:(ny-1), :);
    south_scale = scale(2:ny, :);
    donor_scale(pos) = north_scale(pos);
    donor_scale(~pos) = south_scale(~pos);
    Q_face(1:(ny-1), :, 2) = Qy .* donor_scale;
end
end

function [qG_face, Hf_x, Hf_y, phi_x, phi_y] = sfincs_flux_update( ...
    eta, qG_prev, dx, dt, g, inactive, S)
[ny, nx] = size(eta);
qG_face = zeros(ny, nx, 2, 'like', eta);
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
    active = ~(inactive(:, 1:nx-1) | inactive(:, 2:nx)) & ...
        face.HG > 0 & face.phi > 0 & isfinite(face.n) & face.n > 0 & ...
        isfinite(slope);
    qold = qG_prev(:, 1:nx-1, 1);
    qnew = sfincs_lie_step(qold, face.HG, face.n, face.phi, slope, dt, g, active);
    qG_face(:, 1:nx-1, 1) = qnew;
    Hf_x(:, 1:nx-1) = face.HG;
    phi_x(:, 1:nx-1) = face.phi;
end

if ny > 1
    etaN = eta(1:ny-1, :);
    etaS = eta(2:ny, :);
    zu = max(etaN, etaS);
    face = hp2d_sfincs_velocity_state(S, zu, 'v');
    slope = (etaS - etaN) ./ dx;
    active = ~(inactive(1:ny-1, :) | inactive(2:ny, :)) & ...
        face.HG > 0 & face.phi > 0 & isfinite(face.n) & face.n > 0 & ...
        isfinite(slope);
    qold = qG_prev(1:ny-1, :, 2);
    qnew = sfincs_lie_step(qold, face.HG, face.n, face.phi, slope, dt, g, active);
    qG_face(1:ny-1, :, 2) = qnew;
    Hf_y(1:ny-1, :) = face.HG;
    phi_y(1:ny-1, :) = face.phi;
end
end

function qnew = sfincs_lie_step(qold, HG, nrep, phi, slope, dt, g, active)
qnew = zeros(size(qold), 'like', qold);
den = 1 + g .* dt .* nrep.^2 .* abs(qold) ./ max(HG.^(7/3), eps_like(qold));
rhs = qold - g .* dt .* HG .* slope;
% External forcing F is not used in HydroPol's current local-inertial
% subgrid branch. The SFINCS equation term phi*F*dt is therefore zero.
qnew(active) = rhs(active) ./ den(active);
qnew(~active) = 0;
qnew(~isfinite(qnew)) = 0;
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

function [outlet_flow, outflow, Hf, V, eta, C_a] = apply_sfincs_outlet( ...
    outlet_flow, outflow, Hf, V, eta, C_a, S, row_outlet, col_outlet, ...
    outlet_type, slope_outlet, dx, dt, g, ny, nx)
outlet_sub = sub2ind(size(eta), row_outlet, col_outlet);
side_out = outlet_boundary_side(row_outlet, col_outlet, ny, nx);
Q_out = zeros(size(outlet_sub), 'like', eta);
H_out = zeros(size(outlet_sub), 'like', eta);

for ii = 1:numel(outlet_sub)
    side = side_out{ii};
    face = hp2d_sfincs_boundary_velocity_state(S, eta, side);
    idx = outlet_sub(ii);
    HG = face.HG(idx);
    nrep = face.n(idx);
    active = HG > 0 && isfinite(nrep) && nrep > 0;
    if active
        if outlet_type == 1
            if isscalar(slope_outlet)
                Sout = abs(slope_outlet);
            else
                Sout = abs(slope_outlet(idx));
            end
            qG = (1 / nrep) * HG^(5/3) * sqrt(max(Sout, 0));
        else
            qG = HG * sqrt(g * max(HG, 0));
        end
        Q_out(ii) = qG * dx;
        H_out(ii) = HG;
    end
end

available_Q = max(V(outlet_sub), 0) ./ dt;
Q_out = min(max(Q_out, 0), available_Q);
q_equiv = Q_out ./ (dx^2) .* 1000 .* 3600;
outlet_flow(outlet_sub) = q_equiv;
outflow(:, :, 3) = outlet_flow;
V(outlet_sub) = max(V(outlet_sub) - Q_out .* dt, 0);
eta = hp2d_sfincs_zs_from_cell_volume(S, V);
C_a = hp2d_sfincs_cell_wet_area_from_zs(S, eta);
H3 = Hf(:, :, 3);
H3(outlet_sub) = H_out;
Hf(:, :, 3) = H3;
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
