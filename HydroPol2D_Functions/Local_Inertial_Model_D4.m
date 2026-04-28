function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
    Local_Inertial_Model_D4(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2, ...
    flag_reservoir,z,d_tot,d_p,roughness_cell,roughness_squared,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet, ...
    row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width,River_Depth, ...
    Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev,Subgrid_Properties,flag_overbanks,flag_inflow,SubgridTables)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Last Updated: 5 August, 2024                    %
%                                                                 %
%  Updated here: March 2026                                       %
%  - preserves original behavior when subgrid is OFF              %
%  - avoids unnecessary Hf/Wf computations in shared-face mode    %
%  - does NOT compute C_a inside the face block                   %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

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
    h_min = 1e-6;
else
    h_min = 0;
end

%% Cell depth and WSE
depth_cell = d_tot;
depth_cell = max(depth_cell/1000,0); % m

% Current Stored Volume
if flag_subgrid == 1 && flag_overbanks == 1
    current_volume = nansum(nansum((Resolution - River_Width).*Resolution.*max((depth_cell - River_Depth),0))) + ...
                     nansum(nansum(Resolution.*River_Width.*depth_cell)); % m3
else
    if isscalar(C_a_prev)
        current_volume = nansum(nansum(C_a_prev * depth_cell)); % m3
    else
        current_volume = nansum(nansum(C_a_prev .* depth_cell)); % m3
    end
end

% Keep real DEM untouched
z_dem = z;

% Absolute water surface elevation
y = z_dem + depth_cell;

% Shared-face subgrid only in this specific case
use_sharedface_subgrid = (flag_subgrid == 1 && flag_overbanks ~= 1);

%% ---------------- Faces: hydrostatic reconstruction, slopes, Hf, and Wf ----------------
Sx   = zeros(size(y), 'like', y);
Sy   = zeros(size(y), 'like', y);
Hf_x = zeros(size(y), 'like', y);
Hf_y = zeros(size(y), 'like', y);

if use_sharedface_subgrid
    Wf_x = zeros(size(y), 'like', y);
    Wf_y = zeros(size(y), 'like', y);
else
    % Preserve original coarse behavior when not using shared-face subgrid
    Wf_x = Resolution * ones(size(y), 'like', y);
    Wf_y = Resolution * ones(size(y), 'like', y);
end

% ===== X faces =====
zf_x = max(z_dem(:,1:end-1), z_dem(:,2:end));
etaL_x = max(y(:,1:end-1), zf_x);
etaR_x = max(y(:,2:end),   zf_x);

Sx(:,1:end-1) = (etaR_x - etaL_x) ./ Resolution;
Sx(:,end) = 0;

% ===== Y faces =====
zf_y = max(z_dem(1:end-1,:), z_dem(2:end,:));
etaS_y = max(y(1:end-1,:), zf_y);
etaN_y = max(y(2:end,:),   zf_y);

Sy(2:end,:) = (etaS_y - etaN_y) ./ Resolution;
Sy(1,:) = 0;

if ~use_sharedface_subgrid
    % ----- Original coarse-grid Hf only -----
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

% %% ---------------- Wet/dry and NaN hygiene ----------------
% dry_x = (Hf(:,:,1) <= 0);
% dry_y = (Hf(:,:,2) <= 0);
% 
% nan_x_face = false(size(y), 'like', idx_nan);
% nan_y_face = false(size(y), 'like', idx_nan);
% 
% % x-faces: columns 1:end-1 valid
% nan_x_face(:,1:end-1) = ...
%       isnan(z_dem(:,1:end-1)) | isnan(z_dem(:,2:end)) ...
%     | isnan(y(:,1:end-1))     | isnan(y(:,2:end));
% nan_x_face(:,end) = true;
% 
% % y-faces: rows 2:end valid
% nan_y_face(2:end,:) = ...
%       isnan(z_dem(1:end-1,:)) | isnan(z_dem(2:end,:)) ...
%     | isnan(y(1:end-1,:))     | isnan(y(2:end,:));
% nan_y_face(1,:) = true;
% 
% mask_x = dry_x | nan_x_face;
% mask_y = dry_y | nan_y_face;
% 
% Sx(mask_x) = 0;
% Sy(mask_y) = 0;
% 
% matrix_store(:,:,1) = Sx;
% matrix_store(:,:,2) = Sy;
% 
% tmp = Hf(:,:,1);
% tmp(mask_x) = 0;
% Hf(:,:,1) = tmp;
% 
% tmp = Hf(:,:,2);
% tmp(mask_y) = 0;
% Hf(:,:,2) = tmp;
% 
% if use_sharedface_subgrid
%     tmp = Wf_x;
%     tmp(mask_x) = 0;
%     Wf_x = tmp;
% 
%     tmp = Wf_y;
%     tmp(mask_y) = 0;
%     Wf_y = tmp;
% end

%% ---------------- Minimum operative depth logic (h_min) ----------------
if h_min > 0
    shallow_cell = (depth_cell < h_min);

    shallow_x = false(size(shallow_cell), 'like', idx_nan);
    shallow_y = false(size(shallow_cell), 'like', idx_nan);

    % x-face touches cells j and j+1
    shallow_x(:,1:end-1) = shallow_cell(:,1:end-1) | shallow_cell(:,2:end);
    shallow_x(:,end)     = shallow_cell(:,end);

    % y-face touches cells i-1 and i
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

    if use_sharedface_subgrid
        tmp = Wf_x;
        tmp(shallow_x) = 0;
        Wf_x = tmp;

        tmp = Wf_y;
        tmp(shallow_y) = 0;
        Wf_y = tmp;
    end
end

Hf(Hf < 0) = 0;

%% Local-Inertial Formulation
outflow = outflow/1000/3600*Resolution^2/Resolution; % m2/s
dt = time_step*60; % s
g = 9.81;
outflow_prev = outflow;

%% Sub-grid channels / functions
if flag_subgrid == 1 && flag_overbanks == 1
    % Original overbank subgrid channel mode (unchanged)
    [Q,Qc,Qf,Qci,Qfi,C_a] = subgrid_channel(depth_cell, River_Width, z_dem, z_dem - River_Depth , ...
                                            Resolution, nc, nf, Qc_prev, Qf_prev, Qci_prev, Qfi_prev, ...
                                            g, dt, idx_rivers, outlet_index);

    cell_width = Resolution;
    outflow = Q./(cell_width); % m2/s

    bad = ~isfinite(outflow);
    outflow(bad) = 0;

elseif flag_subgrid == 1
    % Shared-face subgrid functions mode
    Qc = 0;
    Qf = 0;
    Qci = 0;
    Qfi = 0;
    Q = 0;

    % Previous discharge [m3/s] on x and y faces
    Qc_prev = outflow_prev(:,:,1:2) * cell_area / 1000 / 3600;
    Qci_prev = 0*Qc_prev;

    [outflow, C_a, Hf_x, Hf_y, Wf_x, Wf_y] = subgrid_channel_functions( ...
        depth_cell, River_Width, z_dem - River_Depth, Resolution, nc, ...
        Qc_prev, Qci_prev, g, dt, idx_rivers, Subgrid_Properties, SubgridTables);

    Hf(:,:,1) = Hf_x;
    Hf(:,:,2) = Hf_y;

    bad = ~isfinite(outflow);
    if any(bad(:))
        outflow(bad) = 0;
    end

    % kept for outlet logic below
    if isscalar(C_a)
        cell_width = C_a / Resolution;
    else
        cell_width = C_a ./ Resolution;
    end

else
    % Original coarse-grid inertial solver (UNCHANGED)
    Qc = 0;
    Qf = 0;
    Qci = 0;
    Qfi = 0;
    Q = 0;

    outflow = Inertial_Solver(flag_numerical_scheme, outflow_prev, dt, Hf, matrix_store, roughness_squared, Resolution, idx_nan);
    
    bad = ~isfinite(outflow);
    if any(bad(:))
        outflow(bad) = 0;
    end

    C_a = Resolution^2;      % scalar
    cell_width = Resolution; % scalar
end

%% Limiting outflow to critical velocity
if flag_critical == 1
    critical_velocity = Hf(:,:,1:size(outflow,3)) .* sqrt(9.81 * Hf(:,:,1:size(outflow,3)));
    outflow = min(max(outflow, -critical_velocity), critical_velocity);
end

%% Convert q -> discharge rate Q
q = outflow; % m2/s

if use_sharedface_subgrid
    outflow_rate = zeros(size(q), 'like', q);
    outflow_rate(:,:,1) = q(:,:,1) .* Wf_x; % x-face discharge [m3/s]
    outflow_rate(:,:,2) = q(:,:,2) .* Wf_y; % y-face discharge [m3/s]
else
    outflow_rate = outflow .* (cell_width); % original behavior
end

outflow = outflow_rate./(cell_area)*1000*3600; % mm/h
matrix_store = outflow; % mm/h

%% Limiting outflow to maximum velocity
max_velocity = 10; % m/s
threshold_velocity = Hf(:,:,1:size(outflow,3)) * max_velocity;
q = min(max(q, -threshold_velocity), threshold_velocity);

% Recompute discharge AFTER clipping
if use_sharedface_subgrid
    outflow_rate = zeros(size(q), 'like', q);
    outflow_rate(:,:,1) = q(:,:,1) .* Wf_x;
    outflow_rate(:,:,2) = q(:,:,2) .* Wf_y;
else
    outflow_rate = q .* (cell_width);
end

outflow = outflow_rate./(cell_area)*1000*3600;
matrix_store = outflow;

%% Intercell Volume
Vol_Flux = dt * ( ...
    [zeros(ny,1,'like',outflow_rate(:,:,1)) , outflow_rate(:,1:(nx-1),1)] - outflow_rate(:,:,1) ...
    - outflow_rate(:,:,2) + [outflow_rate(2:end,:,2); zeros(1,nx,'like',outflow_rate(:,:,2))] );

if flag_subgrid == 1
    % reserved for future optional subgrid diagonal/additional fluxes
end

%% Reservoir Boundary Condition
if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = d_tot(reservoir_y(ii),reservoir_x(ii))./1000; % m
            dt_h = (time_step)/60; % h

            available_volume = 1000*(max(dtsup - h1(ii),0))/dt_h; % mm/h
            dh = min(k1(ii)*(max(dtsup - h1(ii),0))^k2(ii)/cell_area*1000*3600,available_volume)*dt_h; % mm

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
            dtsup = dtsup - dh/1000;
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
        else
            dh = 0;
        end

        if ~isnan(yds2(ii))
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;

            available_volume = 1000*(max(dtsup - h2(ii),0))/dt_h; % mm/h
            dh = min(k3(ii)*(max(dtsup - h2(ii),0))^k4(ii)/cell_area*1000*3600,available_volume)*dt_h; % mm

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
if isscalar(C_a)
    d_t = d_tot + 1000 * Vol_Flux / C_a;
else
    d_t = d_tot + 1000 * Vol_Flux ./ C_a;
end

if flag_subgrid == 1 && flag_overbanks
    % Inbank -> overbank
    idx = d_t/1000 > River_Depth & d_p/1000 <= River_Depth & idx_rivers;
    if sum(sum(idx)) > 0
        d_t(idx) = 1000*(d_t(idx)/1000 - (d_t(idx)/1000 - z_dem(idx) + (z_dem(idx) - River_Depth(idx))).*(1 - River_Width(idx)/Resolution));
        [d_t(idx)] = inbank_to_overbank(River_Width(idx),River_Depth(idx),d_t(idx)/1000,Resolution);
        C_a(idx) = Resolution^2;
    end

    % Overbank -> inbank
    idx = d_t/1000 <= River_Depth & d_p/1000 >= River_Depth & idx_rivers;
    if sum(sum(idx)) > 0
        factor = d_t(idx)/1000 - z_dem(idx) + (z_dem(idx) - River_Depth(idx));
        d_t(idx) = 1000*(d_t(idx)/1000 + (Resolution*(factor))./River_Width(idx) ...
            - (d_t(idx)/1000 - z_dem(idx) + (z_dem(idx) - River_Depth(idx))));
        [d_t(idx)] = overbank_to_inbank(River_Width(idx),River_Depth(idx),d_t(idx)/1000,Resolution);
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

    % Local water depth at outlet cells [m]
    h_out = max(d_t(outlet_sub), 0) / 1000;

    % Local width at outlet cells
    if isscalar(cell_width)
        width_out = cell_width * ones(size(h_out), 'like', h_out);
    else
        width_out = cell_width(outlet_sub);
    end

    % Local cell area at outlet cells
    if isscalar(cell_area)
        area_out = cell_area * ones(size(h_out), 'like', h_out);
    else
        area_out = cell_area(outlet_sub);
    end

    % Local slope term
    if outlet_type == 1
        % Normal-flow slope prescribed at outlet
        if isscalar(slope_outlet)
            sqrtS_out = sqrt(abs(slope_outlet)) * ones(size(h_out), 'like', h_out);
        else
            sqrtS_out = sqrt(abs(slope_outlet(outlet_sub)));
        end
    else
        % Critical slope only at outlet cells
        h_safe = max(h_out, 1e-12);
        sqrtS_out = sqrt(9.81 .* roughness_squared(outlet_sub) .* h_safe.^(-1/3));
        sqrtS_out(~isfinite(sqrtS_out)) = 0;
    end

    % Keep Hf(:,:,3) only at outlet cells, if you still want it stored
    Hf3 = Hf(:,:,3);
    Hf3(outlet_sub) = h_out;
    Hf(:,:,3) = Hf3;

    % Outlet discharge [mm/h]
    q_out = (1 ./ roughness_cell(outlet_sub)) .* ...
            width_out .* ...
            h_out.^(5/3) .* ...
            sqrtS_out ./ ...
            area_out * 1000 * 3600;

    % Limit by available water
    qmax_out = max((d_t(outlet_sub) - 1000*h_min), 0) / (time_step/60);
    q_out = min(q_out, qmax_out);
    q_out(~isfinite(q_out)) = 0;

    % Write back only at outlet cells
    outlet_flow(outlet_sub) = q_out;
    outflow(:,:,3) = outlet_flow;

    % Update depth only at outlet cells
    d_t(outlet_sub) = d_t(outlet_sub) - q_out * (time_step/60);
end

%% Mass Balance Check
if flag_subgrid == 1 && flag_overbanks == 1
    final_volume = nansum(nansum((Resolution - River_Width).*Resolution.*max((d_t/1000 - River_Depth),0))) + ...
                   nansum(nansum(Resolution.*River_Width.*d_t/1000));
else
    if isscalar(C_a)
        final_volume = nansum(nansum(C_a * d_t / 1000));
    else
        final_volume = nansum(nansum(C_a .* d_t / 1000));
    end
end

% error_vol = (final_volume - current_volume) + nansum(nansum(Resolution^2.*1/1000*time_step/60.*outlet_flow));
% if abs(error_vol) > 10
%     ttt = 1; %#ok<NASGU>
% end

%% Total Flow that Leaves the Cell
mask = outflow;
I_tot_end_cell = abs(sum(mask,3))*dt/1000*1/3600*Resolution^2; % m3

end

function [h] = inbank_to_overbank(B,H,h_p,R)
    h = 1000*(B.*H + (R - B).*h_p) ./ R;
end

function [h] = overbank_to_inbank(B,H,h_p,R)
    h = 1000*(B .* h_p + (R - B) .* (H - h_p)) ./ B;
end

function Vq = hp2d_lookup_uniform(Vtab, q, dz)
% Fast uniform-depth lookup for GPU/CPU
% Vtab: [ny x nx x nz]
% q   : [ny x nx]
% dz  : scalar depth spacing

    [ny, nx, nz] = size(Vtab);
    N = ny * nx;
    maxDepth = (nz - 1) * dz;

    qv = reshape(q, N, 1);
    qv(~isfinite(qv)) = 0;
    qv = max(qv, 0);
    qv = min(qv, maxDepth);

    idxL = floor(qv ./ dz) + 1;
    idxL = min(max(idxL, 1), nz-1);

    q0 = (idxL - 1) .* dz;
    a  = (qv - q0) ./ dz;

    atTop = (qv >= maxDepth);
    idxL(atTop) = nz-1;
    a(atTop)    = 1;

    V = reshape(Vtab, N, nz);

    base = (0:N-1)' * nz;
    ind1 = base + idxL;
    ind2 = ind1 + 1;

    v1 = V(ind1);
    v2 = V(ind2);

    Vqi = v1 + a .* (v2 - v1);
    Vqi(~isfinite(Vqi)) = 0;
    Vqi = max(Vqi, 0);

    Vq = reshape(Vqi, ny, nx);
end