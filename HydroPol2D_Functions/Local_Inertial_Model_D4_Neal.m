function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
    Local_Inertial_Model_D4_Neal(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2, ...
    flag_reservoir,z,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet, ...
    row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,nc,nf,River_Width,River_Depth, ...
    Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev) %#ok<INUSD>
%LOCAL_INERTIAL_MODEL_D4_NEAL Neal et al. (2012) simple channel-subgrid step.

ny = size(z,1);
nx = size(z,2);
g = 9.81;
dt = time_step * 60;

hp2d_neal_validate_geometry(River_Width, River_Depth, Resolution);

depth_cell = max(d_tot ./ 1000, 0);   % representative depth above channel bed [m]
River_Width = max(River_Width, 0);
River_Depth = max(River_Depth, 0);

z_channel = z - River_Depth;
idx_rivers = River_Width > 0 & River_Depth > 0;

V_before = hp2d_neal_cell_volume(depth_cell, River_Width, River_Depth, Resolution);

[Q,Qc,Qf,Qci,Qfi,C_a] = subgrid_channel( ...
    depth_cell, River_Width, z, z_channel, Resolution, nc, nf, ...
    Qc_prev, Qf_prev, Qci_prev, Qfi_prev, g, dt, idx_rivers, outlet_index);

outflow_rate = Q; % [m3/s]
outflow = outflow_rate ./ cell_area * 1000 * 3600; % [mm/h]
outflow(~isfinite(outflow)) = 0;

qout_left  = -[zeros(ny,1,'like',outflow(:,:,1)), outflow(:,1:end-1,1)];
qout_right =  outflow(:,:,1);
qout_up    =  outflow(:,:,2);
qout_down  = -[outflow(2:end,:,2); zeros(1,nx,'like',outflow(:,:,2))];

Vol_Flux = dt * ( ...
    [zeros(ny,1,'like',outflow_rate(:,:,1)) , outflow_rate(:,1:nx-1,1)] - outflow_rate(:,:,1) ...
    - outflow_rate(:,:,2) + [outflow_rate(2:end,:,2); zeros(1,nx,'like',outflow_rate(:,:,2))] );

if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = depth_cell(reservoir_y(ii),reservoir_x(ii));
            dt_h = time_step / 60;

            available_volume = 1000 * max(dtsup - h1(ii), 0) / dt_h;
            dh = min(k1(ii) * (max(dtsup - h1(ii),0))^k2(ii) / cell_area * 1000 * 3600, available_volume) * dt_h;

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh / 1000 * cell_area;
            dtsup = dtsup - dh / 1000;
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
        else
            dh = 0;
        end

        if ~isnan(yds2(ii))
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
            available_volume = 1000 * max(dtsup - h2(ii), 0) / dt_h;
            dh = min(k3(ii) * (max(dtsup - h2(ii),0))^k4(ii) / cell_area * 1000 * 3600, available_volume) * dt_h;

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh / 1000 * cell_area;
            d_tot(yds2(ii),xds2(ii)) = d_tot(yds2(ii),xds2(ii)) + dh;
        end
    end
end

V_after_faces = max(V_before + Vol_Flux, 0);
h_after_faces = hp2d_neal_depth_from_volume(V_after_faces, River_Width, River_Depth, Resolution);
C_a = hp2d_neal_cell_area(h_after_faces, River_Width, River_Depth, Resolution);

outlet_flow = zeros(size(d_tot), 'like', d_tot);
outflow(:,:,3) = 0;
Hf = zeros(ny,nx,3,'like',d_tot);

if ~isempty(row_outlet)
    outlet_sub = sub2ind(size(d_tot), row_outlet(:), col_outlet(:));

    h_out = h_after_faces(outlet_sub);
    w_out = River_Width(outlet_sub);
    H_out = River_Depth(outlet_sub);
    n_channel = nc(outlet_sub);
    n_flood = nf(outlet_sub);
    qc_old_out = zeros(size(h_out), 'like', h_out);
    qf_old_out = zeros(size(h_out), 'like', h_out);

    if ~isempty(Qci_prev) && ~isscalar(Qci_prev)
        qc_old_out = Qci_prev(outlet_sub);
    end
    if ~isempty(Qfi_prev) && ~isscalar(Qfi_prev)
        qf_old_out = Qfi_prev(outlet_sub);
    end

    if isscalar(slope_outlet)
        slope_face = slope_outlet * ones(size(h_out), 'like', h_out);
    else
        slope_face = slope_outlet(outlet_sub);
    end
    slope_face = max(abs(slope_face), 0);

    [Q_out,qc_old_new,qf_old_new] = hp2d_neal_outlet_free_flux( ...
        h_out, w_out, H_out, slope_face, outlet_type, dt, Resolution, g, ...
        n_channel, n_flood, qc_old_out, qf_old_out, d_tolerance / 1000);

    Q_available = V_after_faces(outlet_sub) ./ max(dt, eps);
    Q_uncapped = Q_out;
    Q_out = min(Q_uncapped, Q_available);
    cap_scale = min(Q_out ./ max(Q_uncapped, eps), 1);
    qc_old_new = qc_old_new .* cap_scale;
    qf_old_new = qf_old_new .* cap_scale;
    Q_out(~isfinite(Q_out)) = 0;
    Q_out = max(Q_out, 0);

    outlet_flow(outlet_sub) = Q_out ./ cell_area * 1000 * 3600;
    outflow(:,:,3) = outlet_flow;
    Qci(outlet_sub) = qc_old_new;
    Qfi(outlet_sub) = qf_old_new;

    V_after_faces(outlet_sub) = max(V_after_faces(outlet_sub) - Q_out .* dt, 0);
    h_after_faces = hp2d_neal_depth_from_volume(V_after_faces, River_Width, River_Depth, Resolution);
    C_a = hp2d_neal_cell_area(h_after_faces, River_Width, River_Depth, Resolution);
end

d_t = 1000 .* h_after_faces;
d_t(idx_nan) = NaN;

I_tot_end_cell = abs(sum(outflow,3)) * dt / 1000 / 3600 * Resolution^2;

if flag_critical == 1
    outflow = max(outflow, 0);
end

end

function [Q_out,qc_old_new,qf_old_new] = hp2d_neal_outlet_free_flux( ...
    h_out, w_out, H_out, slope_face, outlet_type, dt, cell_width, g, ...
    n_channel, n_flood, qc_old_prev, qf_old_prev, depth_tolerance_m)
%HP2D_NEAL_OUTLET_FREE_FLUX Neal/LISFLOOD-style free-boundary update.
%   Channel memory is carried in m3/s, floodplain memory in m2/s.

Q_out = zeros(size(h_out), 'like', h_out);
qc_old_new = zeros(size(h_out), 'like', h_out);
qf_old_new = zeros(size(h_out), 'like', h_out);

wet = h_out > depth_tolerance_m;
if ~any(wet)
    return;
end

% A full-width channel is an explicit wide-flow cell. Use the ordinary
% HydroPol2D outlet relation so the limiting case is numerically identical.
has_channel = wet & (w_out > 0) & (H_out > 0);
full_width = has_channel & w_out >= cell_width .* (1 - 10 .* eps(cell_width));
if any(full_width)
    h_full = h_out(full_width);
    n_full = n_channel(full_width);
    if outlet_type == 1
        sqrt_slope = sqrt(slope_face(full_width));
    else
        sqrt_slope = sqrt(g .* n_full.^2 .* max(h_full,eps).^(-1/3));
    end
    q_full = (1 ./ n_full) .* w_out(full_width) .* ...
        h_full.^(5/3) .* sqrt_slope;
    q_full(~isfinite(q_full)) = 0;
    qc_old_new(full_width) = q_full;
    Q_out(full_width) = q_full;
end

% Embedded channel component: same inertial form as the interior solve.
has_channel = has_channel & ~full_width;
if any(has_channel)
    w_channel = w_out(has_channel);
    h_channel = h_out(has_channel);
    A = w_channel .* h_channel;
    R = A ./ max(w_channel + 2 .* h_channel, eps);
    qc_prev = abs(qc_old_prev(has_channel));
    num = qc_prev + abs(g .* dt .* A .* slope_face(has_channel));
    den = 1 + g .* dt .* (n_channel(has_channel) .^ 2) .* qc_prev ./ ...
        max((R .^ (4/3)) .* A, eps);
    qc_now = num ./ den;
    qc_now(~isfinite(qc_now)) = 0;
    qc_old_new(has_channel) = qc_now;
    Q_out(has_channel) = Q_out(has_channel) + qc_now;
end

% Floodplain component: unit discharge with width applied afterward.
hflood = max(h_out - H_out, 0);
active_fp = wet & (hflood > depth_tolerance_m) & (w_out < cell_width);
if any(active_fp)
    qf_prev = abs(qf_old_prev(active_fp));
    num = qf_prev + abs(g .* dt .* hflood(active_fp) .* slope_face(active_fp));
    den = 1 + g .* dt .* (n_flood(active_fp) .^ 2) .* qf_prev ./ ...
        max(hflood(active_fp) .^ (7/3), eps);
    qf_now = num ./ den;
    qf_now(~isfinite(qf_now)) = 0;
    qf_old_new(active_fp) = qf_now;
    Q_out(active_fp) = Q_out(active_fp) + qf_now .* max(cell_width - w_out(active_fp), 0);
end

Q_out(~isfinite(Q_out)) = 0;
Q_out = max(Q_out, 0);
end
