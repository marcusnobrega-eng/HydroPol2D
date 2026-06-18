function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = Routing_Wave_Model_D4(wave_type,flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2,flag_reservoir,z,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width,River_Depth,Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev) %#ok<INUSD>
%ROUTING_WAVE_MODEL_D4 Shared conservative D4 kinematic/diffusive routing.
%
% Sign convention:
%   outflow(:,:,1) is positive from a cell to its right/east neighbor.
%   outflow(:,:,2) is positive from a cell to its up/north neighbor.
%
% Flux outputs follow HydroPol2D's historical convention and are returned
% in mm/h. The internal update uses conservative m3 volume balances.

[ny,nx] = size(z);
dt_total = time_step * 60;
dx = Resolution;

if isempty(outflow) || size(outflow,1) ~= ny || size(outflow,2) ~= nx || size(outflow,3) < 3
    outflow = zeros(ny,nx,3,'like',z);
end

if isempty(idx_nan)
    idx_nan = false(ny,nx);
else
    idx_nan = logical(idx_nan);
end

rough = roughness_cell;
rough(~isfinite(rough) | rough <= 0) = 1e-6;

C_a = resolve_cell_area(C_a_prev,cell_area,dx,ny,nx,z);
C_a(idx_nan) = 0;

Qc = 0;
Qf = 0;
Qci = 0;
Qfi = 0;

d_work = d_tot;
d_work(idx_nan) = NaN;
d_tolerance_m = max(d_tolerance,0) / 1000;
n_sub = estimate_substeps(wave_type,z,d_work,rough,dt_total,dx,idx_nan);
dt_sub = dt_total / n_sub;
time_step_sub = time_step / n_sub;

outlet_volume_total = zeros(ny,nx,'like',z);
last_outflow_rate = zeros(ny,nx,3,'like',z);
last_qx = zeros(ny,nx,'like',z);
last_qy = zeros(ny,nx,'like',z);
Hf = zeros(ny,nx,3,'like',z);

for isub = 1:n_sub
    h = max(d_work,0) / 1000;
    h(idx_nan) = 0;
    eta = z + h;

    [qx,qy,Hf] = compute_face_fluxes(wave_type,z,eta,h,rough,dx,idx_nan,outflow);

    if flag_critical == 1
        qcrit_x = Hf(:,:,1) .* sqrt(9.81 * Hf(:,:,1));
        qcrit_y = Hf(:,:,2) .* sqrt(9.81 * Hf(:,:,2));
        qx = min(max(qx,-qcrit_x),qcrit_x);
        qy = min(max(qy,-qcrit_y),qcrit_y);
    end

    [qx,qy] = limit_to_available_storage(qx,qy,h,C_a,dt_sub,dx,idx_nan);

    outflow_rate = zeros(ny,nx,3,'like',z);
    outflow_rate(:,:,1) = qx .* dx;
    outflow_rate(:,:,2) = qy .* dx;

    Vol_Flux = dt_sub * ( ...
        [zeros(ny,1,'like',z), outflow_rate(:,1:end-1,1)] - outflow_rate(:,:,1) ...
        - outflow_rate(:,:,2) + [outflow_rate(2:end,:,2); zeros(1,nx,'like',z)] );

    if flag_reservoir == 1
        [Vol_Flux,d_work] = apply_reservoir_fluxes(Vol_Flux,d_work,C_a,time_step_sub,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2);
    end

    d_work = d_work + 1000 .* safe_divide(Vol_Flux,C_a);
    d_work(~isfinite(d_work) & ~idx_nan) = 0;
    d_work(d_work < 0) = 0;
    d_work(idx_nan) = NaN;

    [outlet_flow_sub,d_work,Hf(:,:,3)] = apply_outlet_flux(d_work,C_a,rough,dx,time_step_sub,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance_m,idx_nan);
    outlet_volume_total = outlet_volume_total + outlet_flow_sub ./ 1000 ./ 3600 .* C_a .* dt_sub;

    last_outflow_rate = outflow_rate;
    last_qx = qx;
    last_qy = qy;
end

d_t = d_work;
outlet_flow = safe_divide(outlet_volume_total,C_a) * 1000 * 3600 / dt_total;
outlet_flow(~isfinite(outlet_flow)) = 0;
outlet_flow(idx_nan) = 0;

outflow(:,:,1) = safe_divide(last_outflow_rate(:,:,1),C_a) * 1000 * 3600;
outflow(:,:,2) = safe_divide(last_outflow_rate(:,:,2),C_a) * 1000 * 3600;
outflow(:,:,3) = outlet_flow;
outflow(~isfinite(outflow)) = 0;
outflow(repmat(idx_nan,1,1,3)) = 0;

matrix_store = outflow;
qout_left = -[zeros(ny,1,'like',z), matrix_store(:,1:end-1,1)];
qout_right = matrix_store(:,:,1);
qout_up = matrix_store(:,:,2);
qout_down = -[matrix_store(2:end,:,2); zeros(1,nx,'like',z)];

I_tot_end_cell = compute_outgoing_volume(last_qx,last_qy,dx,dt_total,ny,nx,z);
I_tot_end_cell(idx_nan) = 0;
end

function [qx,qy,Hf] = compute_face_fluxes(wave_type,z,eta,h,rough,dx,idx_nan,outflow_template)
[ny,nx] = size(z);
qx = zeros(ny,nx,'like',z);
qy = zeros(ny,nx,'like',z);
Hf = zeros(size(outflow_template),'like',z);

if nx > 1
    zf_x = max(z(:,1:end-1),z(:,2:end));
    eta_left = eta(:,1:end-1);
    eta_right = eta(:,2:end);
    h_left = h(:,1:end-1);
    h_right = h(:,2:end);
    rough_left = rough(:,1:end-1);
    rough_right = rough(:,2:end);

    switch lower(wave_type)
        case 'diffusive'
            drive_x = (eta_left - eta_right) ./ dx;
            hf_x = max(max(eta_left,eta_right) - zf_x,0);
            rough_x = 0.5 * (rough_left + rough_right);
        case 'kinematic'
            drive_x = (z(:,1:end-1) - z(:,2:end)) ./ dx;
            donor_left = drive_x >= 0;
            hf_x = zeros(ny,nx-1,'like',z);
            rough_x = zeros(ny,nx-1,'like',z);
            hf_x(donor_left) = h_left(donor_left);
            hf_x(~donor_left) = h_right(~donor_left);
            rough_x(donor_left) = rough_left(donor_left);
            rough_x(~donor_left) = rough_right(~donor_left);
        otherwise
            error('Unknown routing wave type "%s".',wave_type);
    end

    qx(:,1:end-1) = manning_unit_flux(hf_x,rough_x,drive_x);
    Hf(:,1:end-1,1) = hf_x;
end

if ny > 1
    zf_y = max(z(1:end-1,:),z(2:end,:));
    eta_north = eta(1:end-1,:);
    eta_south = eta(2:end,:);
    h_north = h(1:end-1,:);
    h_south = h(2:end,:);
    rough_north = rough(1:end-1,:);
    rough_south = rough(2:end,:);

    switch lower(wave_type)
        case 'diffusive'
            drive_y = (eta_south - eta_north) ./ dx;
            hf_y = max(max(eta_north,eta_south) - zf_y,0);
            rough_y = 0.5 * (rough_north + rough_south);
        case 'kinematic'
            drive_y = (z(2:end,:) - z(1:end-1,:)) ./ dx;
            donor_south = drive_y >= 0;
            hf_y = zeros(ny-1,nx,'like',z);
            rough_y = zeros(ny-1,nx,'like',z);
            hf_y(donor_south) = h_south(donor_south);
            hf_y(~donor_south) = h_north(~donor_south);
            rough_y(donor_south) = rough_south(donor_south);
            rough_y(~donor_south) = rough_north(~donor_south);
        otherwise
            error('Unknown routing wave type "%s".',wave_type);
    end

    qy(2:end,:) = manning_unit_flux(hf_y,rough_y,drive_y);
    Hf(2:end,:,2) = hf_y;
end

qx = close_inactive_faces_x(qx,idx_nan);
qy = close_inactive_faces_y(qy,idx_nan);
end

function n_sub = estimate_substeps(wave_type,z,d_work,rough,dt,dx,idx_nan)
h = max(d_work,0) / 1000;
h(idx_nan) = 0;
active = h > 0;
if ~any(active,'all')
    n_sub = 1;
    return
end

eta = z + h;
slopes = [];
if size(z,2) > 1
    if strcmpi(wave_type,'diffusive')
        slopes = [slopes; reshape(abs((eta(:,1:end-1) - eta(:,2:end)) ./ dx),[],1)]; %#ok<AGROW>
    else
        slopes = [slopes; reshape(abs((z(:,1:end-1) - z(:,2:end)) ./ dx),[],1)]; %#ok<AGROW>
    end
end
if size(z,1) > 1
    if strcmpi(wave_type,'diffusive')
        slopes = [slopes; reshape(abs((eta(2:end,:) - eta(1:end-1,:)) ./ dx),[],1)]; %#ok<AGROW>
    else
        slopes = [slopes; reshape(abs((z(2:end,:) - z(1:end-1,:)) ./ dx),[],1)]; %#ok<AGROW>
    end
end

slope_max = max(slopes(:),[],'omitnan');
if ~isfinite(slope_max) || slope_max <= 0
    n_sub = 1;
    return
end

h_max = max(h(active),[],'omitnan');
n_min = min(rough(active),[],'omitnan');
n_min = max(n_min,1e-6);
velocity_est = h_max^(2/3) / n_min * sqrt(slope_max);
courant_target = 0.20;
n_sub = max(1,ceil(velocity_est * dt / dx / courant_target));
n_sub = min(n_sub,200);
end

function q = manning_unit_flux(hf,rough,drive_slope)
hf = max(hf,0);
drive_slope(~isfinite(drive_slope)) = 0;
rough(~isfinite(rough) | rough <= 0) = 1e-6;
q = sign(drive_slope) .* hf.^(5/3) ./ rough .* sqrt(abs(drive_slope));
q(hf <= 0 | abs(drive_slope) <= 0) = 0;
q(~isfinite(q)) = 0;
end

function C_a = resolve_cell_area(C_a_prev,cell_area,dx,ny,nx,prototype)
if isempty(C_a_prev)
    area_in = cell_area;
else
    area_in = C_a_prev;
end

if isempty(area_in)
    C_a = dx^2 * ones(ny,nx,'like',prototype);
elseif isscalar(area_in)
    C_a = area_in * ones(ny,nx,'like',prototype);
else
    C_a = area_in;
end

C_a(~isfinite(C_a) | C_a <= 0) = dx^2;
end

function qx = close_inactive_faces_x(qx,idx_nan)
[~,nx] = size(idx_nan);
if nx > 1
    closed = idx_nan(:,1:end-1) | idx_nan(:,2:end);
    tmp = qx(:,1:end-1);
    tmp(closed) = 0;
    qx(:,1:end-1) = tmp;
end
qx(:,end) = 0;
end

function qy = close_inactive_faces_y(qy,idx_nan)
[ny,~] = size(idx_nan);
if ny > 1
    closed = idx_nan(1:end-1,:) | idx_nan(2:end,:);
    tmp = qy(2:end,:);
    tmp(closed) = 0;
    qy(2:end,:) = tmp;
end
qy(1,:) = 0;
end

function [qx,qy] = limit_to_available_storage(qx,qy,h,C_a,dt,dx,idx_nan)
[ny,nx] = size(h);
available = max(h,0) .* C_a;
available(idx_nan) = 0;
outgoing = compute_outgoing_volume(qx,qy,dx,dt,ny,nx,h);
scale = ones(ny,nx,'like',h);
active = outgoing > available & outgoing > 0;
scale(active) = available(active) ./ outgoing(active);
scale(idx_nan) = 0;

if nx > 1
    face = qx(:,1:end-1);
    src_scale = scale(:,1:end-1);
    right_scale = scale(:,2:end);
    src_scale(face < 0) = right_scale(face < 0);
    qx(:,1:end-1) = face .* src_scale;
end

if ny > 1
    face = qy(2:end,:);
    src_scale = scale(2:end,:);
    north_scale = scale(1:end-1,:);
    src_scale(face < 0) = north_scale(face < 0);
    qy(2:end,:) = face .* src_scale;
end
end

function outgoing = compute_outgoing_volume(qx,qy,dx,dt,ny,nx,prototype)
outgoing = zeros(ny,nx,'like',prototype);
if nx > 1
    Qx = qx(:,1:end-1) .* dx;
    outgoing(:,1:end-1) = outgoing(:,1:end-1) + max(Qx,0) .* dt;
    outgoing(:,2:end) = outgoing(:,2:end) + max(-Qx,0) .* dt;
end
if ny > 1
    Qy = qy(2:end,:) .* dx;
    outgoing(2:end,:) = outgoing(2:end,:) + max(Qy,0) .* dt;
    outgoing(1:end-1,:) = outgoing(1:end-1,:) + max(-Qy,0) .* dt;
end
end

function [Vol_Flux,d_tot] = apply_reservoir_fluxes(Vol_Flux,d_tot,C_a,time_step,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2)
dt_h = time_step / 60;
for ii = 1:length(reservoir_y)
    if ~isnan(yds1(ii))
        area_res = C_a(reservoir_y(ii),reservoir_x(ii));
        dtsup = d_tot(reservoir_y(ii),reservoir_x(ii)) / 1000;
        available_mm_h = 1000 * max(dtsup - h1(ii),0) / max(dt_h,eps);
        q_mm_h = k1(ii) * max(dtsup - h1(ii),0)^k2(ii) / area_res * 1000 * 3600;
        dh = min(q_mm_h,available_mm_h) * dt_h;
        Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh / 1000 * area_res;
        dtsup = dtsup - dh / 1000;
        d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
    else
        dh = 0; %#ok<NASGU>
        dtsup = d_tot(reservoir_y(ii),reservoir_x(ii)) / 1000;
    end

    if ~isnan(yds2(ii))
        area_res = C_a(reservoir_y(ii),reservoir_x(ii));
        available_mm_h = 1000 * max(dtsup - h2(ii),0) / max(dt_h,eps);
        q_mm_h = k3(ii) * max(dtsup - h2(ii),0)^k4(ii) / area_res * 1000 * 3600;
        dh = min(q_mm_h,available_mm_h) * dt_h;
        Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh / 1000 * area_res;
        d_tot(yds2(ii),xds2(ii)) = d_tot(yds2(ii),xds2(ii)) + dh;
    end
end
end

function [outlet_flow,d_t,Hf_outlet] = apply_outlet_flux(d_t,C_a,rough,dx,time_step,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance_m,idx_nan)
[ny,nx] = size(d_t);
dt_h = time_step / 60;
outlet_flow = zeros(ny,nx,'like',d_t);
Hf_outlet = zeros(ny,nx,'like',d_t);
if isempty(row_outlet)
    return
end

mask_outlet = false(ny,nx);
mask_outlet(sub2ind([ny,nx],row_outlet(:),col_outlet(:))) = true;
mask_outlet = mask_outlet & logical(outlet_index) & ~idx_nan;
if ~any(mask_outlet,'all')
    return
end

if outlet_type == 1
    if isscalar(slope_outlet)
        S0 = zeros(ny,nx,'like',d_t);
        S0(mask_outlet) = slope_outlet;
    else
        S0 = slope_outlet;
        S0(~mask_outlet) = 0;
    end
else
    h_out_tmp = max(d_t,0) / 1000;
    S0 = zeros(ny,nx,'like',d_t);
    S0(mask_outlet) = h_out_tmp(mask_outlet).^(-1/6) .* sqrt(9.81) .* rough(mask_outlet);
    S0(~isfinite(S0)) = 0;
end

h_out = max(d_t,0) / 1000;
Hf_outlet(mask_outlet) = h_out(mask_outlet);
Q_out = zeros(ny,nx,'like',d_t);
Q_out(mask_outlet) = 1 ./ rough(mask_outlet) .* dx .* h_out(mask_outlet).^(5/3) .* sqrt(abs(S0(mask_outlet)));
outlet_flow(mask_outlet) = safe_divide(Q_out(mask_outlet),C_a(mask_outlet)) * 1000 * 3600;
available_mm_h = max(d_t - 1000 * d_tolerance_m,0) / max(dt_h,eps);
outlet_flow = min(outlet_flow,available_mm_h);
outlet_flow(~isfinite(outlet_flow)) = 0;
d_t = d_t - outlet_flow * dt_h;
d_t(d_t < 0) = 0;
d_t(idx_nan) = NaN;
end

function y = safe_divide(a,b)
y = zeros(size(a),'like',a);
valid = b ~= 0 & isfinite(b);
y(valid) = a(valid) ./ b(valid);
end
