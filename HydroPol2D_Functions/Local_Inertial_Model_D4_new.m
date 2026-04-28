function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a,LI_Work] = ...
    Local_Inertial_Model_D4(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2, ...
    flag_reservoir,z,Precomputed,LI_Work,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type, ...
    slope_outlet,row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf, ...
    River_Width,River_Depth,Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev,Subgrid_Properties,flag_overbanks, ...
    flag_inflow,SubgridTables)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Updated for precomputed geometry: March 2026    %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% =========================================================================
% Expected precomputed fields (to be built outside this function later)
% -------------------------------------------------------------------------
% Precomputed.ny                  : number of rows
% Precomputed.nx                  : number of columns
% Precomputed.zf_x                : max(z(:,1:end-1), z(:,2:end))   [ny x (nx-1)]
% Precomputed.zf_y                : max(z(1:end-1,:), z(2:end,:))   [(ny-1) x nx]
% Precomputed.valid_x             : valid x-face mask               [ny x (nx-1)]
% Precomputed.valid_y             : valid y-face mask               [(ny-1) x nx]
% Precomputed.active_mask         : active cell mask                [ny x nx]   (optional)
% Precomputed.outlet_linear_index : linear indices of outlet cells  (optional)
% Precomputed.compute_error_vol   : true/false for debug MB check   (optional)
% =========================================================================

%% -------------------- Basic dimensions / constants -------------------- %%
if isfield(Precomputed,'ny') && ~isempty(Precomputed.ny)
    ny = Precomputed.ny;
else
    ny = size(z,1);
end

if isfield(Precomputed,'nx') && ~isempty(Precomputed.nx)
    nx = Precomputed.nx;
else
    nx = size(z,2);
end

if ~isfield(Precomputed,'compute_error_vol') || isempty(Precomputed.compute_error_vol)
    Precomputed.compute_error_vol = false;
end

% For the coarse-grid / non-subgrid case we are optimizing here
C_a_const  = Resolution^2;   % [m2]
cell_width = Resolution;     % [m]

% Minimum depth logic used to suppress spurious wet-dry fluxes
if flag_inflow ~= 1
    h_min = 1e-6;
else
    h_min = 0;
end

%% -------------------- Current water depth / storage ------------------- %%
depth_cell = max(d_tot / 1000, 0);    % [m]
y          = z + depth_cell;          % water surface elevation [m]

if Precomputed.compute_error_vol
    if isfield(Precomputed,'active_mask') && ~isempty(Precomputed.active_mask)
        current_volume = sum((depth_cell(Precomputed.active_mask) .* cell_area), 'all');
    else
        current_volume = nansum(nansum(depth_cell .* cell_area));
    end
else
    current_volume = 0;
end

%% ---------------- Face reconstruction using precomputed bed faces ----- %%
% X faces: (i, j+1/2)
zf_x   = Precomputed.zf_x;                              % [ny x (nx-1)]
etaL_x = max(y(:,1:end-1), zf_x);
etaR_x = max(y(:,2:end  ), zf_x);
Sx_int = (etaR_x - etaL_x) ./ Resolution;              % positive to the east
Hf_x_int = max(max(y(:,1:end-1), y(:,2:end)) - zf_x, 0);

% Y faces: (i+1/2, j)
zf_y   = Precomputed.zf_y;                              % [(ny-1) x nx]
etaS_y = max(y(1:end-1,:), zf_y);
etaN_y = max(y(2:end  ,:), zf_y);
Sy_int = (etaS_y - etaN_y) ./ Resolution;              % positive to the north
Hf_y_int = max(max(y(1:end-1,:), y(2:end,:)) - zf_y, 0);

% Preallocate full-face arrays without concatenation
Sx = Precomputed.Zeros;
Sy = Precomputed.Zeros;
Hf = Precomputed.Zeros_3D;

Sx(:,1:nx-1)   = Sx_int;
Sy(2:ny,:)     = Sy_int;
Hf(:,1:nx-1,1) = Hf_x_int;
Hf(2:ny,:,2)   = Hf_y_int;

%% ---------------- Wet/dry and invalid-face masking -------------------- %%
% Dynamic finite masks from current free surface
finite_x = isfinite(y(:,1:end-1)) & isfinite(y(:,2:end));
finite_y = isfinite(y(1:end-1,:)) & isfinite(y(2:end,:));

% Build full-size dynamic finite masks
finite_x_full = false(ny,nx,'like',idx_nan);
finite_y_full = false(ny,nx,'like',idx_nan);

finite_x_full(:,1:nx-1) = finite_x;
finite_y_full(2:ny,:)   = finite_y;

% Full-size dry masks
dry_x_full = false(ny,nx,'like',idx_nan);
dry_y_full = false(ny,nx,'like',idx_nan);

dry_x_full(:,1:nx-1) = (Hf_x_int <= 0);
dry_y_full(2:ny,:)   = (Hf_y_int <= 0);

% Static + dynamic face masks
mask_x_full = ~(Precomputed.valid_x_full & finite_x_full) | dry_x_full;
mask_y_full = ~(Precomputed.valid_y_full & finite_y_full) | dry_y_full;

% Minimum operative depth around shallow cells
if h_min > 0
    shallow_cell = (depth_cell < h_min) | ~isfinite(depth_cell);

    shallow_x_full = false(ny,nx,'like',idx_nan);
    shallow_y_full = false(ny,nx,'like',idx_nan);

    shallow_x_full(:,1:nx-1) = shallow_cell(:,1:end-1) | shallow_cell(:,2:end);
    shallow_y_full(2:ny,:)   = shallow_cell(1:end-1,:) | shallow_cell(2:end,:);

    mask_x_full = mask_x_full | shallow_x_full;
    mask_y_full = mask_y_full | shallow_y_full;
end

% Apply masks directly on full-size face arrays
Sx(mask_x_full) = 0;
Sy(mask_y_full) = 0;

tmp = Hf(:,:,1);
tmp(mask_x_full) = 0;
Hf(:,:,1) = tmp;

tmp = Hf(:,:,2);
tmp(mask_y_full) = 0;
Hf(:,:,2) = tmp;

% Final guard
Hf(Hf < 0) = 0;

% matrix_store carries slopes for the inertial solver
matrix_store = zeros(ny,nx,3,'like',y);
matrix_store(:,:,1) = Sx;
matrix_store(:,:,2) = Sy;

%% -------------------- Local inertial formulation ---------------------- %%
% Input outflow arrives in mm/h; convert to unit discharge [m2/s]
outflow     = outflow / 1000 / 3600 * cell_width;
dt          = time_step * 60;  % [s]
g           = 9.81;
outflow_prev = outflow;

% Subgrid is intentionally ignored here for optimization of the standard case
Qc  = 0;
Qf  = 0;
Qci = 0;
Qfi = 0;
C_a = Precomputed.C_a_full;

[outflow_solver] = Inertial_Solver(flag_numerical_scheme,outflow_prev(:,:,1:2),dt,Hf(:,:,1:2),matrix_store(:,:,1:2),Precomputed.n_face_sq,Resolution,idx_nan);

% Rebuild full outflow container:
%   (:,:,1) -> x-face flux
%   (:,:,2) -> y-face flux
%   (:,:,3) -> outlet flux (filled later)
outflow = zeros(ny,nx,3,'like',Hf);
outflow(:,:,1:2) = outflow_solver;

% Hygiene
outflow(isnan(outflow)) = 0;

%% -------------------- Optional critical-flow limiter ------------------ %%
if flag_critical == 1
    critical_velocity = Hf(:,:,1:size(outflow,3)) .* sqrt(g * Hf(:,:,1:size(outflow,3))); % [m2/s]
    outflow(outflow >  critical_velocity) =  critical_velocity(outflow >  critical_velocity);
    outflow(outflow < -critical_velocity) = -critical_velocity(outflow < -critical_velocity);
end

%% -------------------- Absolute max-velocity limiter ------------------- %%
max_velocity = 10;
threshold_velocity = Hf(:,:,1:size(outflow,3)) * max_velocity; % [m2/s]

outflow(outflow >  threshold_velocity) =  threshold_velocity(outflow >  threshold_velocity);
outflow(outflow < -threshold_velocity) = -threshold_velocity(outflow < -threshold_velocity);

% Convert to volumetric flux [m3/s]
outflow_rate = outflow * cell_width;

% Keep a compact container for output-face bookkeeping
matrix_store = outflow;

%% -------------------- Intercell volume balance ------------------------ %%
% X-direction uses face 1, Y-direction uses face 2
Vol_Flux = LI_Work.Vol_Flux;
Vol_Flux(:) = 0;

% Start with outflows leaving each cell
Vol_Flux = Vol_Flux - outflow_rate(:,:,1);   % flow leaving to the right
Vol_Flux = Vol_Flux - outflow_rate(:,:,2);   % flow leaving upward

% Add inflow from left neighbor (their right-face outflow enters current cell)
Vol_Flux(:,2:nx) = Vol_Flux(:,2:nx) + outflow_rate(:,1:nx-1,1);

% Add inflow from lower neighbor in the y-bookkeeping used here
Vol_Flux(1:ny-1,:) = Vol_Flux(1:ny-1,:) + outflow_rate(2:ny,:,2);

% Convert from m3/s to m3 over the time step
Vol_Flux = dt * Vol_Flux;   % [m3]

%% -------------------- Reservoir boundary conditions ------------------- %%
if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = d_tot(reservoir_y(ii),reservoir_x(ii)) / 1000;  % [m]
            dt_h  = time_step / 60;                                 % [h]

            available_volume = 1000 * max(dtsup - h1(ii), 0) / dt_h; % [mm/h]
            dh = min(k1(ii) * (max(dtsup - h1(ii),0))^k2(ii) / cell_area * 1000 * 3600, available_volume) * dt_h; % [mm]

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000 * cell_area;
            dtsup = dtsup - dh/1000;

            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
        else
            dh = 0;
        end

        if ~isnan(yds2(ii))
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;

            available_volume = 1000 * max(dtsup - h2(ii), 0) / dt_h; % [mm/h]
            dh = min(k3(ii) * (max(dtsup - h2(ii),0))^k4(ii) / cell_area * 1000 * 3600, available_volume) * dt_h; % [mm]

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000 * cell_area;
            d_tot(yds2(ii),xds2(ii)) = d_tot(yds2(ii),xds2(ii)) + dh;
        end
    end
end

%% -------------------- Output fluxes by direction ---------------------- %%
qout_left  = LI_Work.qout_left;  qout_left(:)  = 0;
qout_right = LI_Work.qout_right; qout_right(:) = 0;
qout_up    = LI_Work.qout_up;    qout_up(:)    = 0;
qout_down  = LI_Work.qout_down;  qout_down(:)  = 0;

% Right-face outflow from the left neighbor enters as "left-directed" bookkeeping
qout_left(:,2:nx) = -matrix_store(:,1:nx-1,1);

% Rightward outflow from current cell
qout_right(:,:) = matrix_store(:,:,1);

% Upward outflow from current cell
qout_up(:,:) = matrix_store(:,:,2);

% Downward bookkeeping from the cell above/below arrangement used here
qout_down(1:ny-1,:) = -matrix_store(2:ny,:,2);

%% -------------------- Cell depth update before outlet ----------------- %%
d_t = d_tot + 1000 * Vol_Flux / C_a_const;   % [mm]

lost_mass = 1/1000 * abs(sum(sum(d_t(d_t<0) * C_a_const))); %#ok<NASGU>
%% -------------------- Outlet flow ------------------------------------- %%
if outlet_type == 1
    S_0 = slope_outlet .* outlet_index;
else
    S_0 = zeros(size(z),'like',z);

    if isfield(Precomputed,'outlet_linear_index') && ~isempty(Precomputed.outlet_linear_index)
        h_out = max(d_t(Precomputed.outlet_linear_index)/1000, 0);
        S_0(Precomputed.outlet_linear_index) = (roughness_cell(Precomputed.outlet_linear_index).^2) .* g .* h_out.^(-1/3);
    else
        for i = 1:length(row_outlet)
            row = row_outlet(i);
            col = col_outlet(i);
            h   = max(d_t(row,col)/1000, 0);
            S_0(row,col) = (roughness_cell(row,col)^2) * g * h^(-1/3);
        end
    end

    S_0(isinf(S_0)) = 0;
end

mask_outlet = Precomputed.mask_outlet_cells & (S_0 ~= 0);

Hf_outlet = zeros(ny,nx,'like',z);
Hf_outlet(mask_outlet) = max(d_t(mask_outlet),0)/1000;
Hf(:,:,3) = Hf_outlet;

outlet_flow = (1 ./ roughness_cell) .* cell_width .* Hf(:,:,3).^(5/3) .* abs(S_0).^(0.5) ./ C_a_const * 1000 * 3600; % [mm/h]
outlet_flow = min(outlet_flow, max((d_t - 1000*h_min),0) / (time_step/60));
outlet_flow(~mask_outlet) = 0;

outflow(:,:,3) = outlet_flow;

% Final depth after outlet release
d_t = d_t - outlet_flow * (time_step/60);

d_t(d_t < 0) = 0;
d_t(idx_nan) = nan;

%% -------------------- Optional debug mass-balance check --------------- %%
if Precomputed.compute_error_vol
    if isfield(Precomputed,'active_mask') && ~isempty(Precomputed.active_mask)
        final_volume = sum((d_t(Precomputed.active_mask)/1000 .* cell_area), 'all');
    else
        final_volume = nansum(nansum((d_t/1000) .* cell_area));
    end

    error_vol = (final_volume - current_volume) + nansum(nansum(C_a_const .* 1/1000 .* time_step/60 .* outlet_flow)); %#ok<NASGU>

    if abs(error_vol) > 10
        ttt = 1; %#ok<NASGU>
    end
end

%% -------------------- Total flow leaving each cell -------------------- %%
I_tot_end_cell = abs(sum(outflow,3)) * dt / 1000 / 3600 .* cell_area; % [m3]

if max(max(d_t)) > 20000
    ttt = 1;
end
end

function h = inbank_to_overbank(B,H,h_p,R)
h = 1000 * (B .* H + (R - B) .* h_p) ./ R;
end

function h = overbank_to_inbank(B,H,h_p,R)
h = 1000 * (B .* h_p + (R - B) .* (H - h_p)) ./ B;
end