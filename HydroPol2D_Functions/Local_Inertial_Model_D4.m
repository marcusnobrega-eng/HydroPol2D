function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
    Local_Inertial_Model_D4(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2, ...
    flag_reservoir,z,d_tot,d_p,roughness_cell,roughness_squared,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet, ...
    row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width,River_Depth, ...
    Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev,Subgrid_Properties,flag_overbanks,flag_inflow,SubgridTables)

%% Domain Dimensions
if isgpuarray(cell_area)
    d_t_min = gpuArray(d_tolerance/1000); %#ok<NASGU>
    ny = size(z,1);
    nx = size(z,2);
else
    d_t_min = d_tolerance/1000; %#ok<NASGU>
    ny = size(z,1);
    nx = size(z,2);
end

River_Width(logical(outlet_index)) = 0;
River_Depth(logical(outlet_index)) = 0;

idx_rivers = River_Width > 0;

%% Minimum operative depth
if flag_inflow ~= 1
    h_min = 0;
else
    h_min = 0;
end

%% Cell depth and WSE
depth_cell = max(d_tot/1000,0); % [m]

% Shared-face subgrid only in this specific case
use_sharedface_subgrid = (flag_subgrid == 1 && flag_overbanks ~= 1);
if use_sharedface_subgrid && ~isempty(SubgridTables)
    z_dem = SubgridTables.invert_el;
else
    z_dem = z;
end

% Absolute water surface elevation. In lookup-subgrid mode, d_tot is the
% representative depth above the subgrid cell invert.
y = z_dem + depth_cell;

%% Current stored volume [m3]
if flag_subgrid == 1 && flag_overbanks == 1

    current_volume = nansum(nansum((Resolution - River_Width).*Resolution.*max((depth_cell - River_Depth),0))) + ...
                     nansum(nansum(Resolution.*River_Width.*depth_cell));

elseif use_sharedface_subgrid && ~isempty(SubgridTables)

    V_now = hp2d_lookup_eta_uniform( ...
        SubgridTables.volume_cell, ...
        y, ...
        SubgridTables.invert_el, ...
        SubgridTables.dz, ...
        SubgridTables.maxDepth);

    current_volume = nansum(V_now(:)); %#ok<NASGU>

else

    if isscalar(C_a_prev)
        current_volume = nansum(nansum(C_a_prev * depth_cell)); %#ok<NASGU>
    else
        current_volume = nansum(nansum(C_a_prev .* depth_cell)); %#ok<NASGU>
    end

end

%% ---------------- Faces: hydrostatic reconstruction, slopes, Hf, and Wf ----------------
Sx   = zeros(size(y), 'like', y);
Sy   = zeros(size(y), 'like', y);
Hf_x = zeros(size(y), 'like', y);
Hf_y = zeros(size(y), 'like', y);

if use_sharedface_subgrid
    Wf_x = zeros(size(y), 'like', y);
    Wf_y = zeros(size(y), 'like', y);
else
    Wf_x = Resolution * ones(size(y), 'like', y);
    Wf_y = Resolution * ones(size(y), 'like', y);
end

%% ============================================================
% SLOPES: WITH / WITHOUT HYDROSTATIC RECONSTRUCTION
% ============================================================
flag_HC = 0;
% These are needed in BOTH schemes because Hf uses them
zf_x = max(z_dem(:,1:end-1), z_dem(:,2:end));
zf_y = max(z_dem(1:end-1,:), z_dem(2:end,:));

if flag_HC == 1

    % --------------------------------------------------------
    % HYDROSTATIC RECONSTRUCTION
    % Robust near wet/dry fronts, but more diffusive
    % --------------------------------------------------------

    % X faces
    etaL_x = max(y(:,1:end-1), zf_x);
    etaR_x = max(y(:,2:end),   zf_x);

    Sx(:,1:end-1) = (etaR_x - etaL_x) ./ Resolution;
    Sx(:,end) = 0;

    % Y faces
    etaS_y = max(y(1:end-1,:), zf_y);
    etaN_y = max(y(2:end,:),   zf_y);

    Sy(2:end,:) = (etaS_y - etaN_y) ./ Resolution;
    Sy(1,:) = 0;

else

    % --------------------------------------------------------
    % ORIGINAL LOCAL INERTIAL
    % Best for smooth benchmark / rainfall-runoff waves
    % --------------------------------------------------------

    % X faces
    Sx(:,1:end-1) = (y(:,2:end) - y(:,1:end-1)) ./ Resolution;
    Sx(:,end) = 0;

    % Y faces
    Sy(2:end,:) = (y(1:end-1,:) - y(2:end,:)) ./ Resolution;
    Sy(1,:) = 0;

end

if ~use_sharedface_subgrid
    Hf_x(:,1:end-1) = max(max(y(:,1:end-1), y(:,2:end)) - zf_x, 0);
    Hf_y(2:end,:)   = max(max(y(1:end-1,:), y(2:end,:)) - zf_y, 0);

    Wf_x(:,end) = 0;
    Wf_y(1,:)   = 0;
end

% Boundary faces always closed
Hf_x(:,end) = 0;
Hf_y(1,:)   = 0;

% Store for solver
matrix_store = 0*outflow;
matrix_store(:,:,1) = Sx;
matrix_store(:,:,2) = Sy;

Hf = 0*outflow;

if ~use_sharedface_subgrid
    Hf(:,:,1) = Hf_x;
    Hf(:,:,2) = Hf_y;
end

%% Minimum operative depth logic
if h_min > 0
    shallow_cell = (depth_cell < h_min);

    shallow_x = false(size(shallow_cell), 'like', idx_nan);
    shallow_y = false(size(shallow_cell), 'like', idx_nan);

    shallow_x(:,1:end-1) = shallow_cell(:,1:end-1) | shallow_cell(:,2:end);
    shallow_x(:,end)     = shallow_cell(:,end);

    shallow_y(2:end,:) = shallow_cell(1:end-1,:) | shallow_cell(2:end,:);
    shallow_y(1,:)     = shallow_cell(1,:);

    Sx(shallow_x) = 0;
    Sy(shallow_y) = 0;

    matrix_store(:,:,1) = Sx;
    matrix_store(:,:,2) = Sy;

    tmp = Hf(:,:,1);
    tmp(shallow_x) = 0;
    Hf(:,:,1) = tmp;

    tmp = Hf(:,:,2);
    tmp(shallow_y) = 0;
    Hf(:,:,2) = tmp;
end

Hf(Hf < 0) = 0;

%% Local-Inertial Formulation
outflow_mmh_prev = outflow;
if use_sharedface_subgrid
    Q_prev = zeros(ny,nx,2,'like',depth_cell);
    if isscalar(cell_area)
        Q_prev(:,:,1:2) = outflow_mmh_prev(:,:,1:2) ./ 1000 ./ 3600 .* cell_area;
    else
        Q_prev(:,:,1) = outflow_mmh_prev(:,:,1) ./ 1000 ./ 3600 .* cell_area;
        Q_prev(:,:,2) = outflow_mmh_prev(:,:,2) ./ 1000 ./ 3600 .* cell_area;
    end
    outflow_prev = [];
else
    outflow = outflow/1000/3600*Resolution^2/Resolution; % [m2/s]
    outflow_prev = outflow;
end
dt = time_step*60; % [s]

%% Sub-grid channels / functions
if flag_subgrid == 1 && flag_overbanks == 1

    % Original overbank subgrid channel mode
    [Q,Qc,Qf,Qci,Qfi,C_a] = subgrid_channel(depth_cell, River_Width, z_dem, z_dem - River_Depth , ...
                                            Resolution, nc, nf, Qc_prev, Qf_prev, Qci_prev, Qfi_prev, ...
                                            g, dt, idx_rivers, outlet_index);

    cell_width = Resolution;
    outflow = Q./cell_width; % [m2/s]
    outflow(~isfinite(outflow)) = 0;

elseif flag_subgrid == 1

    % Shared-face lookup-table subgrid mode
    if isempty(SubgridTables)
        error('flag_subgrid = 1 requires SubgridTables for shared-face subgrid mode.');
    end

    Qc  = 0;
    Qf  = 0;
    Qci = 0;
    Qfi = 0;
    Q   = 0;

    [Q_face, Hf_x, Hf_y, Wf_x, Wf_y] = hp2d_subgrid_local_inertial_flux( ...
        y,Q_prev,nc,Resolution,dt,h_min,idx_nan,SubgridTables,flag_numerical_scheme);

    Hf(:,:,1) = Hf_x;
    Hf(:,:,2) = Hf_y;

    outflow = zeros(ny,nx,2,'like',depth_cell);
    outflow(:,:,1) = Q_face(:,:,1) ./ max(Wf_x, eps_like(Wf_x));
    outflow(:,:,2) = Q_face(:,:,2) ./ max(Wf_y, eps_like(Wf_y));
    outflow(:,:,1) = zero_dry_width(outflow(:,:,1), Wf_x);
    outflow(:,:,2) = zero_dry_width(outflow(:,:,2), Wf_y);
    outflow(~isfinite(outflow)) = 0;

    % Current wetted cell area for outlet and diagnostic use
    C_a = hp2d_lookup_eta_uniform( ...
        SubgridTables.area_cell, ...
        y, ...
        SubgridTables.invert_el, ...
        SubgridTables.dz, ...
        SubgridTables.maxDepth);

    C_a(~isfinite(C_a)) = 0;

    cell_width = C_a ./ Resolution;
    cell_width(~isfinite(cell_width)) = 0;

else

    % Original coarse-grid inertial solver
    Qc  = 0;
    Qf  = 0;
    Qci = 0;
    Qfi = 0;
    Q   = 0;

    outflow = Inertial_Solver(flag_numerical_scheme, outflow_prev, dt, Hf, matrix_store, roughness_squared, Resolution, idx_nan);
    outflow(~isfinite(outflow)) = 0;

    C_a = Resolution^2;
    cell_width = Resolution;

end

%% Limiting outflow to critical velocity
if flag_critical == 1
    critical_velocity = Hf(:,:,1:size(outflow,3)) .* sqrt(9.81 * Hf(:,:,1:size(outflow,3)));
    outflow = min(max(outflow, -critical_velocity), critical_velocity);
end

%% Convert q -> discharge rate Q
q = outflow; % [m2/s]

if use_sharedface_subgrid
    outflow_rate = zeros(size(q), 'like', q);
    outflow_rate(:,:,1) = q(:,:,1) .* Wf_x;
    outflow_rate(:,:,2) = q(:,:,2) .* Wf_y;
else
    outflow_rate = outflow .* cell_width;
end

outflow = outflow_rate./cell_area*1000*3600; % [mm/h]
matrix_store = outflow;

%% Limiting outflow to maximum velocity
max_velocity = 10; % [m/s]
threshold_velocity = Hf(:,:,1:size(outflow,3)) * max_velocity;
q = min(max(q, -threshold_velocity), threshold_velocity);

% Recompute discharge after clipping
if use_sharedface_subgrid
    outflow_rate = zeros(size(q), 'like', q);
    outflow_rate(:,:,1) = q(:,:,1) .* Wf_x;
    outflow_rate(:,:,2) = q(:,:,2) .* Wf_y;
    outflow_rate = apply_subgrid_draining_limiter( ...
        outflow_rate, y, dt, SubgridTables, Resolution);
	else
	    outflow_rate = q .* cell_width;
	    outflow_rate = apply_coarse_draining_limiter( ...
	        outflow_rate, d_tot, dt, cell_area);
	end

outflow = outflow_rate./cell_area*1000*3600;
matrix_store = outflow;

%% Intercell Volume
Vol_Flux = dt * ( ...
    [zeros(ny,1,'like',outflow_rate(:,:,1)) , outflow_rate(:,1:(nx-1),1)] - outflow_rate(:,:,1) ...
    - outflow_rate(:,:,2) + [outflow_rate(2:end,:,2); zeros(1,nx,'like',outflow_rate(:,:,2))] );

%% Reservoir Boundary Condition
if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = d_tot(reservoir_y(ii),reservoir_x(ii))./1000;
            dt_h = time_step/60;

            available_volume = 1000*(max(dtsup - h1(ii),0))/dt_h;
            dh = min(k1(ii)*(max(dtsup - h1(ii),0))^k2(ii)/cell_area*1000*3600,available_volume)*dt_h;

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
            dtsup = dtsup - dh/1000;
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
        else
            dh = 0;
        end

        if ~isnan(yds2(ii))
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;

            available_volume = 1000*(max(dtsup - h2(ii),0))/dt_h;
            dh = min(k3(ii)*(max(dtsup - h2(ii),0))^k4(ii)/cell_area*1000*3600,available_volume)*dt_h;

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
            d_tot(yds2(ii),xds2(ii)) = d_tot(yds2(ii),xds2(ii)) + dh;
        end
    end
end

%% Output Fluxes
qout_left  = -[zeros(ny,1,'like',matrix_store(:,:,1)), matrix_store(:,1:end-1,1)];
qout_right =  matrix_store(:,:,1);
qout_up    =  matrix_store(:,:,2);
qout_down  = -[matrix_store(2:end,:,2); zeros(1,nx,'like',matrix_store(:,:,2))];

%% Final depth at the cell
if use_sharedface_subgrid && ~isempty(SubgridTables)

    V_now = hp2d_lookup_eta_uniform( ...
        SubgridTables.volume_cell, ...
        y, ...
        SubgridTables.invert_el, ...
        SubgridTables.dz, ...
        SubgridTables.maxDepth);

    V_new = V_now + Vol_Flux;
    V_new(~isfinite(V_new)) = 0;
    V_new = max(V_new, 0);

    drep_new = hp2d_inverse_volume_uniform( ...
        SubgridTables.volume_cell, ...
        V_new, ...
        SubgridTables.dz, ...
        SubgridTables.maxDepth, ...
        Resolution);

    eta_new = SubgridTables.invert_el + drep_new;

    d_t = 1000 .* max(eta_new - z_dem, 0);
    d_t(idx_nan) = NaN;

    C_a = hp2d_lookup_eta_uniform( ...
        SubgridTables.area_cell, ...
        eta_new, ...
        SubgridTables.invert_el, ...
        SubgridTables.dz, ...
        SubgridTables.maxDepth);

    C_a(~isfinite(C_a)) = 0;

else

    if isscalar(C_a)
        d_t = d_tot + 1000 * Vol_Flux / C_a;
    else
        d_t = d_tot + 1000 * Vol_Flux ./ C_a;
    end

end

if flag_subgrid == 1 && flag_overbanks
    idx = d_t/1000 > River_Depth & d_p/1000 <= River_Depth & idx_rivers;
    if sum(sum(idx)) > 0
        d_t(idx) = 1000*(d_t(idx)/1000 - (d_t(idx)/1000 - z_dem(idx) + (z_dem(idx) - River_Depth(idx))).*(1 - River_Width(idx)/Resolution));
        d_t(idx) = inbank_to_overbank(River_Width(idx),River_Depth(idx),d_t(idx)/1000,Resolution);
        C_a(idx) = Resolution^2;
    end

    idx = d_t/1000 <= River_Depth & d_p/1000 >= River_Depth & idx_rivers;
    if sum(sum(idx)) > 0
        factor = d_t(idx)/1000 - z_dem(idx) + (z_dem(idx) - River_Depth(idx));
        d_t(idx) = 1000*(d_t(idx)/1000 + (Resolution*(factor))./River_Width(idx) ...
            - (d_t(idx)/1000 - z_dem(idx) + (z_dem(idx) - River_Depth(idx))));
        d_t(idx) = overbank_to_inbank(River_Width(idx),River_Depth(idx),d_t(idx)/1000,Resolution);
        C_a(idx) = River_Width(idx)*Resolution;
    end
end

d_t(d_t<0) = 0;

%% Outlet Flow
outlet_flow = zeros(size(d_t), 'like', d_t);
outflow(:,:,3) = 0;
Hf(:,:,3) = 0;

if ~isempty(row_outlet)
    outlet_sub = sub2ind(size(d_t), row_outlet(:), col_outlet(:));

    h_out = max(d_t(outlet_sub), 0) / 1000;

    if use_sharedface_subgrid
        [q_out, h_face_out, area_out] = subgrid_outlet_flow_mmh( ...
            d_t, outlet_sub, row_outlet(:), col_outlet(:), slope_outlet, ...
            outlet_type, roughness_squared, roughness_cell, time_step, ...
            h_min, Resolution, SubgridTables, cell_area);
        Hf3 = Hf(:,:,3);
        Hf3(outlet_sub) = h_face_out;
        Hf(:,:,3) = Hf3;

        outlet_flow(outlet_sub) = q_out;
        outflow(:,:,3) = outlet_flow;

        V_before = hp2d_subgrid_lookup_depth( ...
            SubgridTables.volume_cell, max(d_t ./ 1000, 0), ...
            SubgridTables.dz, SubgridTables.maxDepth);
        outlet_volume = outlet_flow(outlet_sub) .* (time_step/60) ./ 1000 .* area_out(outlet_sub);
        V_after = max(V_before, 0);
        V_after(outlet_sub) = max(V_before(outlet_sub) - outlet_volume, 0);
        drep_after = hp2d_subgrid_inverse_volume( ...
            SubgridTables.volume_cell, V_after, ...
            SubgridTables.dz, SubgridTables.maxDepth, Resolution);
        d_t = 1000 .* drep_after;
        d_t(idx_nan) = NaN;
        C_a = hp2d_subgrid_lookup_depth( ...
            SubgridTables.area_cell, max(d_t ./ 1000, 0), ...
            SubgridTables.dz, SubgridTables.maxDepth);
        C_a(~isfinite(C_a)) = 0;
    else
    if isscalar(cell_width)
        width_out = cell_width * ones(size(h_out), 'like', h_out);
    else
        width_out = cell_width(outlet_sub);
    end

    if isscalar(cell_area)
        area_out = cell_area * ones(size(h_out), 'like', h_out);
    else
        area_out = cell_area(outlet_sub);
    end

    if outlet_type == 1
        if isscalar(slope_outlet)
            sqrtS_out = sqrt(abs(slope_outlet)) * ones(size(h_out), 'like', h_out);
        else
            sqrtS_out = sqrt(abs(slope_outlet(outlet_sub)));
        end
    else
        h_safe = max(h_out, 1e-12);
        sqrtS_out = sqrt(9.81 .* roughness_squared(outlet_sub) .* h_safe.^(-1/3));
        sqrtS_out(~isfinite(sqrtS_out)) = 0;
    end

    Hf3 = Hf(:,:,3);
    Hf3(outlet_sub) = h_out;
    Hf(:,:,3) = Hf3;

    q_out = (1 ./ roughness_cell(outlet_sub)) .* ...
            width_out .* ...
            h_out.^(5/3) .* ...
            sqrtS_out ./ ...
            area_out * 1000 * 3600;

    qmax_out = max((d_t(outlet_sub) - 1000*h_min), 0) / (time_step/60);
    q_out = min(q_out, qmax_out);
    q_out(~isfinite(q_out)) = 0;

    outlet_flow(outlet_sub) = q_out;
    outflow(:,:,3) = outlet_flow;

    d_t(outlet_sub) = d_t(outlet_sub) - q_out * (time_step/60);
    end
end

%% Mass Balance Check
if flag_subgrid == 1 && flag_overbanks == 1
    final_volume = nansum(nansum((Resolution - River_Width).*Resolution.*max((d_t/1000 - River_Depth),0))) + ...
                   nansum(nansum(Resolution.*River_Width.*d_t/1000)); %#ok<NASGU>
else
    if isscalar(C_a)
        final_volume = nansum(nansum(C_a * d_t / 1000)); %#ok<NASGU>
    else
        final_volume = nansum(nansum(C_a .* d_t / 1000)); %#ok<NASGU>
    end
end

%% Total Flow that Leaves the Cell
mask = outflow;
I_tot_end_cell = abs(sum(mask,3))*dt/1000*1/3600*Resolution^2; % [m3]

end

function h = inbank_to_overbank(B,H,h_p,R)
    h = 1000*(B.*H + (R - B).*h_p) ./ R;
end

function h = overbank_to_inbank(B,H,h_p,R)
    h = 1000*(B .* h_p + (R - B) .* (H - h_p)) ./ B;
end

function Vq = hp2d_lookup_eta_uniform(Vtab, eta, invert_el, dz, maxDepth)

[ny,nx,nz] = size(Vtab);

Vq = NaN(size(eta),'like',Vtab);

d = eta - invert_el;

valid = isfinite(d) & isfinite(invert_el);
dry   = valid & (d <= 0);

d_clamped = d;
d_clamped(valid) = max(d_clamped(valid), 0);
d_clamped(valid) = min(d_clamped(valid), maxDepth);

idxL = floor(d_clamped ./ dz) + 1;
idxL = max(idxL, 1);
idxL = min(idxL, nz-1);
idxU = idxL + 1;

dL = (idxL - 1) .* dz;
w  = (d_clamped - dL) ./ dz;

atTop = valid & (d_clamped >= maxDepth);
idxL(atTop) = nz - 1;
idxU(atTop) = nz;
w(atTop)    = 1;

[I,J] = ndgrid(1:ny, 1:nx);

indL = sub2ind([ny,nx,nz], I, J, idxL);
indU = sub2ind([ny,nx,nz], I, J, idxU);

vL = Vtab(indL);
vU = Vtab(indU);

ok = valid & isfinite(vL) & isfinite(vU) & isfinite(w);

Vq(ok) = vL(ok) + w(ok) .* (vU(ok) - vL(ok));
Vq(dry) = 0;

end

function drep = hp2d_inverse_volume_uniform(Vtab, Vq, dz, maxDepth, Resolution)

[ny,nx,nz] = size(Vtab);

drep = NaN(ny,nx,'like',Vq);

depth_axis = reshape((0:nz-1) * dz, 1, 1, nz);

finite_count = sum(isfinite(Vtab), 3);
Vtop = Vtab(:,:,end);

valid_cell = (finite_count >= 2) & isfinite(Vtop) & isfinite(Vq);

if ~any(valid_cell(:))
    return;
end

Vq_work = Vq;
Vq_work(~valid_cell) = NaN;
Vq_work(valid_cell) = max(Vq_work(valid_cell), 0);

Vq_clamped = Vq_work;
Vq_clamped(valid_cell) = min(Vq_work(valid_cell), Vtop(valid_cell));

cmp = (Vtab <= Vq_clamped);
cmp(~isfinite(Vtab)) = false;

idxL = sum(cmp, 3);
idxL = max(idxL, 1);
idxL = min(idxL, nz-1);
idxU = idxL + 1;

[I,J] = ndgrid(1:ny, 1:nx);

indL = sub2ind([ny,nx,nz], I, J, idxL);
indU = sub2ind([ny,nx,nz], I, J, idxU);

V1 = Vtab(indL);
V2 = Vtab(indU);

D = repmat(depth_axis, ny, nx, 1);
d1 = D(indL);
d2 = D(indU);

V1 = reshape(V1, ny, nx);
V2 = reshape(V2, ny, nx);
d1 = reshape(d1, ny, nx);
d2 = reshape(d2, ny, nx);

interp_ok = valid_cell & isfinite(V1) & isfinite(V2) & isfinite(d1) & isfinite(d2);

denom = V2 - V1;
bad = interp_ok & (abs(denom) <= eps(class_underlying_like(Vtab)));

denom_safe = denom;
denom_safe(bad) = 1;

a = zeros(ny,nx,'like',Vq);
a(interp_ok) = (Vq_clamped(interp_ok) - V1(interp_ok)) ./ denom_safe(interp_ok);
a(interp_ok) = max(min(a(interp_ok), 1), 0);
a(bad) = 0;

drep(interp_ok) = d1(interp_ok) + a(interp_ok) .* (d2(interp_ok) - d1(interp_ok));

above_top = valid_cell & (Vq_work > Vtop);

if any(above_top(:))
    drep(above_top) = maxDepth + ...
        (Vq_work(above_top) - Vtop(above_top)) ./ (Resolution^2);
end

zero_or_neg = valid_cell & (Vq <= 0);
drep(zero_or_neg) = 0;

bad_final = valid_cell & ~isfinite(drep);
drep(bad_final) = 0;

end

function c = class_underlying_like(x)

if isa(x,'gpuArray')
    c = classUnderlying(x);
else
    c = class(x);
end

end

function q = zero_dry_width(q, W)
q(W <= 0 | ~isfinite(W)) = 0;
end

function e = eps_like(x)
if isa(x, 'gpuArray')
    e = eps(classUnderlying(x));
else
    e = eps(class(x));
end
end

function [q_out,h_face_out,area_map] = subgrid_outlet_flow_mmh( ...
    d_t, outlet_sub, row_outlet, col_outlet, slope_outlet, outlet_type, ...
    roughness_squared, roughness_cell, time_step, h_min, Resolution, ...
    SubgridTables, cell_area)

ny = size(d_t,1);
nx = size(d_t,2);
dt_h = time_step / 60;
g = 9.81;

drep = max(d_t ./ 1000, 0);
cell_state = hp2d_subgrid_cell_state(SubgridTables, drep, 'depth', Resolution);
area_map = cell_state.area;
area_map(~isfinite(area_map) | area_map <= 0) = cell_area;
volume_map = cell_state.volume;

q_out = zeros(numel(outlet_sub), 1, 'like', d_t);
h_face_out = zeros(numel(outlet_sub), 1, 'like', d_t);

for kk = 1:numel(outlet_sub)
    r = row_outlet(kk);
    c = col_outlet(kk);
    eta = SubgridTables.invert_el(r,c) + drep(r,c);

    if nx > 1 && c >= nx
        [A,Rh,n_face,K] = lookup_outlet_face(SubgridTables, 'x', r, nx-1, eta);
    elseif nx > 1 && c <= 1
        [A,Rh,n_face,K] = lookup_outlet_face(SubgridTables, 'x', r, 1, eta);
    elseif ny > 1 && r <= 1
        [A,Rh,n_face,K] = lookup_outlet_face(SubgridTables, 'y', 1, c, eta);
    elseif ny > 1 && r >= ny
        [A,Rh,n_face,K] = lookup_outlet_face(SubgridTables, 'y', ny-1, c, eta);
    elseif nx > 1
        [A,Rh,n_face,K] = lookup_outlet_face(SubgridTables, 'x', r, min(c, nx-1), eta);
    else
        A = area_map(r,c);
        Rh = max(drep(r,c), 0);
        n_face = roughness_cell(r,c);
        K = (1 / max(n_face, eps_like(n_face))) * A * Rh^(2/3);
    end

    if ~isfinite(A) || A <= 0 || ~isfinite(Rh) || Rh <= 0
        continue;
    end

    if outlet_type == 1
        if isscalar(slope_outlet)
            sqrtS = sqrt(abs(slope_outlet));
        elseif isequal(size(slope_outlet), size(d_t))
            sqrtS = sqrt(abs(slope_outlet(r,c)));
        else
            sqrtS = sqrt(abs(slope_outlet(1)));
        end
        Q = K * sqrtS;
    else
        n2 = roughness_squared(r,c);
        if ~isfinite(n2) || n2 <= 0
            n2 = max(n_face, eps_like(n_face))^2;
        end
        h_safe = max(A / Resolution, 1e-12);
        sqrtS = sqrt(g .* n2 .* h_safe.^(-1/3));
        Q = (1 / max(n_face, eps_like(n_face))) * A * Rh^(2/3) * sqrtS;
    end

    A_cell = max(area_map(r,c), eps_like(area_map));
    q_candidate = Q / A_cell * 1000 * 3600;
    qmax = max(volume_map(r,c) - h_min * A_cell, 0) / A_cell / dt_h * 1000;
    q_out(kk) = min(max(q_candidate, 0), qmax);
    h_face_out(kk) = A / max(Resolution, eps_like(A));
end
end

function Q = apply_subgrid_draining_limiter(Q, eta, dt, SubgridTables, Resolution)
%APPLY_SUBGRID_DRAINING_LIMITER Prevent intercell fluxes from overdrawing storage.
%
% Q(:,:,1) is positive from (i,j) to (i,j+1). Q(:,:,2) follows the HydroPol2D
% continuity convention and is positive from (i,j) to (i-1,j). The limiter
% scales all outgoing fluxes from a donor cell when their explicit
% time-integrated volume exceeds available stored water in the subgrid cell
% volume table.

[ny,nx,~] = size(Q);
cell_state = hp2d_subgrid_cell_state(SubgridTables, eta, 'eta', Resolution);
available_rate = max(cell_state.volume, 0) ./ max(dt, eps_like(eta));

out_rate = zeros(ny,nx,'like',eta);
Qx = Q(:,:,1);
Qy = Q(:,:,2);

out_rate = out_rate + max(Qx, 0);
out_rate(:,2:end) = out_rate(:,2:end) + max(-Qx(:,1:end-1), 0);
out_rate = out_rate + max(Qy, 0);
out_rate(1:end-1,:) = out_rate(1:end-1,:) + max(-Qy(2:end,:), 0);

scale = ones(ny,nx,'like',eta);
needs_limit = out_rate > available_rate & out_rate > 0;
scale(needs_limit) = available_rate(needs_limit) ./ out_rate(needs_limit);
scale(~isfinite(scale)) = 0;
scale = max(min(scale, 1), 0);

if nx > 1
    donor_scale_pos = scale(:,1:nx-1);
    donor_scale_neg = scale(:,2:nx);
    pos = Qx(:,1:nx-1) > 0;
    neg = Qx(:,1:nx-1) < 0;
    Qx_sub = Qx(:,1:nx-1);
    Qx_sub(pos) = Qx_sub(pos) .* donor_scale_pos(pos);
    Qx_sub(neg) = Qx_sub(neg) .* donor_scale_neg(neg);
    Qx(:,1:nx-1) = Qx_sub;
    Qx(:,nx) = 0;
end

if ny > 1
    donor_scale_pos = scale(2:ny,:);
    donor_scale_neg = scale(1:ny-1,:);
    pos = Qy(2:ny,:) > 0;
    neg = Qy(2:ny,:) < 0;
    Qy_sub = Qy(2:ny,:);
    Qy_sub(pos) = Qy_sub(pos) .* donor_scale_pos(pos);
    Qy_sub(neg) = Qy_sub(neg) .* donor_scale_neg(neg);
    Qy(2:ny,:) = Qy_sub;
    Qy(1,:) = 0;
end

Q(:,:,1) = Qx;
Q(:,:,2) = Qy;
Q(~isfinite(Q)) = 0;
end

function Q = apply_coarse_draining_limiter(Q, d_tot_mm, dt, cell_area)
%APPLY_COARSE_DRAINING_LIMITER Prevent non-subgrid LI fluxes from overdrawing cells.
%
% Q(:,:,1) is positive from (i,j) to (i,j+1). Q(:,:,2) follows the HydroPol2D
% continuity convention and is positive from (i,j) to (i-1,j). Q is discharge
% [m3/s]. The limiter scales outgoing fluxes from each donor cell if their
% explicit volume exceeds current coarse-cell storage.

[ny,nx,~] = size(Q);
if isscalar(cell_area)
    area = cell_area .* ones(ny,nx,'like',d_tot_mm);
else
    area = cell_area;
end

available_rate = max(d_tot_mm ./ 1000, 0) .* area ./ max(dt, eps_like(d_tot_mm));

out_rate = zeros(ny,nx,'like',d_tot_mm);
Qx = Q(:,:,1);
Qy = Q(:,:,2);

out_rate = out_rate + max(Qx, 0);
out_rate(:,2:end) = out_rate(:,2:end) + max(-Qx(:,1:end-1), 0);
out_rate = out_rate + max(Qy, 0);
out_rate(1:end-1,:) = out_rate(1:end-1,:) + max(-Qy(2:end,:), 0);

scale = ones(ny,nx,'like',d_tot_mm);
needs_limit = out_rate > available_rate & out_rate > 0;
scale(needs_limit) = available_rate(needs_limit) ./ out_rate(needs_limit);
scale(~isfinite(scale)) = 0;
scale = max(min(scale, 1), 0);

if nx > 1
    donor_scale_pos = scale(:,1:nx-1);
    donor_scale_neg = scale(:,2:nx);
    pos = Qx(:,1:nx-1) > 0;
    neg = Qx(:,1:nx-1) < 0;
    Qx_sub = Qx(:,1:nx-1);
    Qx_sub(pos) = Qx_sub(pos) .* donor_scale_pos(pos);
    Qx_sub(neg) = Qx_sub(neg) .* donor_scale_neg(neg);
    Qx(:,1:nx-1) = Qx_sub;
    Qx(:,nx) = 0;
end

if ny > 1
    donor_scale_pos = scale(2:ny,:);
    donor_scale_neg = scale(1:ny-1,:);
    pos = Qy(2:ny,:) > 0;
    neg = Qy(2:ny,:) < 0;
    Qy_sub = Qy(2:ny,:);
    Qy_sub(pos) = Qy_sub(pos) .* donor_scale_pos(pos);
    Qy_sub(neg) = Qy_sub(neg) .* donor_scale_neg(neg);
    Qy(2:ny,:) = Qy_sub;
    Qy(1,:) = 0;
end

Q(:,:,1) = Qx;
Q(:,:,2) = Qy;
Q(~isfinite(Q)) = 0;
end

function [A,Rh,n_face,K] = lookup_outlet_face(T, dirn, r, c, eta)
if dirn == 'x'
    invert = T.invert_x(r,c);
    d = max(eta - invert, 0);
    A = lookup_scalar_depth(T.area_x, r, c, d, T.dz, T.maxDepth);
    Rh = lookup_scalar_depth(T.Rh_x, r, c, d, T.dz, T.maxDepth);
    n_face = lookup_scalar_depth(T.n_x, r, c, d, T.dz, T.maxDepth);
    K = lookup_scalar_depth(T.K_x, r, c, d, T.dz, T.maxDepth);
else
    invert = T.invert_y(r,c);
    d = max(eta - invert, 0);
    A = lookup_scalar_depth(T.area_y, r, c, d, T.dz, T.maxDepth);
    Rh = lookup_scalar_depth(T.Rh_y, r, c, d, T.dz, T.maxDepth);
    n_face = lookup_scalar_depth(T.n_y, r, c, d, T.dz, T.maxDepth);
    K = lookup_scalar_depth(T.K_y, r, c, d, T.dz, T.maxDepth);
end
end

function v = lookup_scalar_depth(tab, r, c, d, dz, maxDepth)
nz = size(tab,3);
d = min(max(d, 0), maxDepth);
idxL = floor(d / dz) + 1;
idxL = max(1, min(idxL, nz-1));
idxU = idxL + 1;
w = (d - (idxL - 1) * dz) / dz;
if d >= maxDepth
    idxL = nz - 1;
    idxU = nz;
    w = 1;
end
vL = tab(r,c,idxL);
vU = tab(r,c,idxU);
if isfinite(vL) && isfinite(vU)
    v = vL + w * (vU - vL);
else
    v = NaN;
end
end
