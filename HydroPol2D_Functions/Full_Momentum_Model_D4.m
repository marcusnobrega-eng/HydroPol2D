function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
    Full_Momentum_Model_D4(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2, ...
    flag_reservoir,z,d_tot,d_p,roughness_cell,roughness_squared,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet, ...
    row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width,River_Depth, ...
    Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev,Subgrid_Properties,flag_overbanks,flag_inflow,SubgridTables,Outlet_Properties)
%FULL_MOMENTUM_MODEL_D4
% -------------------------------------------------------------------------
% Conservative 2-D full-momentum shallow-water solver for HydroPol2D-style
% depth updates.
%
% This function keeps the same call signature and main outputs used by your
% Local_Inertial_Model_D4, but replaces the local-inertial flux computation
% with a robust finite-volume full-momentum solver:
%
%   U = [h, hu, hv]^T
%
%   dU/dt + dF(U)/dx + dG(U)/dy = S_bed + S_friction
%
% where:
%   h  = water depth [m]
%   hu = unit discharge in x direction [m^2/s], positive east/right
%   hv = unit discharge in y direction [m^2/s], positive north/up
%
% Robustness features:
%   1) Hydrostatic reconstruction for bed-slope source terms.
%   2) HLL Riemann flux for full momentum.
%   3) Positivity-preserving draining limiter for wetting/drying fronts.
%   4) Semi-implicit Manning friction update.
%   5) DEM NaNs and idx_nan cells are treated as inactive solid walls.
%   6) Closed reflective/no-normal-flow boundaries everywhere except the
%      outlet, which is handled as an explicit sink using normal/critical
%      flow logic equivalent to the local-inertial version.
%
% Important storage convention:
%   outflow(:,:,1) = x-face mass flux converted to HydroPol2D mm/h units
%   outflow(:,:,2) = y-face mass flux converted to HydroPol2D mm/h units
%   outflow(:,:,3) = outlet sink [mm/h]
%   outflow(:,:,4) = updated cell-centered hu [m^2/s]
%   outflow(:,:,5) = updated cell-centered hv [m^2/s]
%
% You must keep passing the returned outflow back into the next call. Pages
% 4 and 5 are what make this a true full-momentum time-marching solver.
%
% Notes:
%   - Lookup-table subgrid mode uses representative depth above subgrid
%     invert for momentum, but applies continuity in volume space through
%     SubgridTables.volume_cell.
%   - flag_numerical_scheme is accepted for interface compatibility. The
%     robust HLL + hydrostatic-reconstruction scheme is used by default.
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Compatibility / unused-input guards
% -------------------------------------------------------------------------
% These inputs are kept because the function is called with the same
% HydroPol2D routing interface as Local_Inertial_Model_D4.
% The current coarse-grid full-momentum version does not use the subgrid
% structures, but keeping them explicit prevents input-list ambiguity.
% -------------------------------------------------------------------------

if nargin < 45 || nargin > 46
    error('Full_Momentum_Model_D4 expects 45 or 46 inputs, but received %d. Check the call in HydroPol2D_Main_While.', nargin);
end

if nargin < 46
    Outlet_Properties = [];
end

if isempty(Subgrid_Properties)
    Subgrid_Properties = [];
end

if isempty(SubgridTables)
    SubgridTables = [];
end

% These are intentionally accepted for interface compatibility.
% Do not delete them from the function signature.
unused_flag_numerical_scheme = flag_numerical_scheme; %#ok<NASGU>
unused_nc                    = nc;                    %#ok<NASGU>
unused_nf                    = nf;                    %#ok<NASGU>
unused_River_Width           = River_Width;           %#ok<NASGU>
unused_River_Depth           = River_Depth;           %#ok<NASGU>
unused_Qc_prev               = Qc_prev;               %#ok<NASGU>
unused_Qf_prev               = Qf_prev;               %#ok<NASGU>
unused_Qci_prev              = Qci_prev;              %#ok<NASGU>
unused_Qfi_prev              = Qfi_prev;              %#ok<NASGU>
unused_C_a_prev              = C_a_prev;              %#ok<NASGU>
unused_Subgrid_Properties    = Subgrid_Properties;    %#ok<NASGU>
unused_flag_inflow           = flag_inflow;           %#ok<NASGU>


%% ------------------------------------------------------------------------
% 0. Basic setup
% -------------------------------------------------------------------------
ny = size(z,1);
nx = size(z,2);

g  = 9.81;
dx = Resolution;
dt = time_step * 60;          % [s]

% Dry tolerance used by the full-momentum solver [m].
% Do not use this to limit the outlet sink; the outlet can drain to h_min.
h_dry = max(d_tolerance/1000, 1.0e-12);
h_min = 0.0;

%% ------------------------------------------------------------------------
% Source/friction treatment switches
% -------------------------------------------------------------------------
% flag_HR_full_momentum:
%   1 = hydrostatic reconstruction for complex terrain
%   0 = raw HLL fluxes + explicit/forecast bed-slope source
%
% flag_source_friction_predictor:
%   1 = apply bed-slope acceleration + exact implicit Manning friction
%       BEFORE computing HLL fluxes.
%       This is the recommended test for normal roughness.
%
% flag_post_flux_friction_corrector:
%   0 = no second friction damping after the flux update
%   1 = optional light post-flux correction
% -------------------------------------------------------------------------

flag_HR_full_momentum = 0;

flag_source_friction_predictor = 1;
flag_post_flux_friction_corrector = 0;
post_flux_friction_fraction = 0.50;

% Optional override: set Outlet_Properties.flag_HR_full_momentum = 1
% to activate hydrostatic reconstruction for complex terrain runs.
if exist('Outlet_Properties','var') && isstruct(Outlet_Properties) && ...
        isfield(Outlet_Properties,'flag_HR_full_momentum') && ...
        ~isempty(Outlet_Properties.flag_HR_full_momentum)
    flag_HR_full_momentum = Outlet_Properties.flag_HR_full_momentum;
end

% Safety velocity limiter. This is a numerical protection, not a boundary
% condition. You can increase this after verification if needed.
max_velocity = 10.0;          % [m/s]

% Coarse-grid full-momentum version.
Qc  = 0; Qf  = 0; Qci = 0; Qfi = 0;

use_lookup_subgrid = flag_subgrid == 1 && flag_overbanks ~= 1 && ~isempty(SubgridTables);
if flag_subgrid == 1 && flag_overbanks == 1
    warning('HydroPol2D:FullMomentumLegacySubgrid', ...
        ['Full_Momentum_Model_D4 received flag_subgrid=1 and flag_overbanks=1. ', ...
         'This legacy overbank pathway is diagnostic only; lookup-table subgrid validation uses flag_overbanks=0.']);
end

if use_lookup_subgrid
    hp2d_validate_subgrid_tables(SubgridTables);
    if flag_reservoir == 1
        error('HydroPol2D:FullMomentumSubgridReservoirUnsupported', ...
            'Full-momentum lookup subgrid storage is not yet coupled to reservoir structures.');
    end
end

% Cell area matrix.
if isscalar(cell_area)
    C_a = cell_area * ones(size(z), 'like', z);
else
    C_a = cell_area;
end
C_a(~isfinite(C_a) | C_a <= 0) = dx^2;

% Active domain. DEM NaNs are solid inactive cells.
if use_lookup_subgrid
    z_dem = SubgridTables.invert_el;
    active = isfinite(z_dem) & ~idx_nan;
else
    z_dem = z;
    active = isfinite(z) & ~idx_nan;
end

if ~isempty(outlet_index)
    % Outlet cells must remain active if DEM is valid; outlet is handled as
    % a sink, not as a dry/closed wall.
    active(logical(outlet_index) & isfinite(z_dem)) = true;
end

h     = max(d_tot/1000, 0);   % [m]
h(~active) = 0;

V_subgrid_old = [];
if use_lookup_subgrid
    V_subgrid_old = hp2d_subgrid_lookup_depth( ...
        SubgridTables.volume_cell, h, SubgridTables.dz, SubgridTables.maxDepth);
    C_a = hp2d_subgrid_lookup_depth( ...
        SubgridTables.area_cell, h, SubgridTables.dz, SubgridTables.maxDepth);
    C_a(~isfinite(C_a) | C_a <= 0) = dx^2;
end

eta = z_dem + h;
eta(~active) = NaN;

%% ------------------------------------------------------------------------
% Face-based reflective wall masks
% -------------------------------------------------------------------------
% Current flux arrays only represent INTERNAL active-active faces.
% For full momentum, inactive/outside faces that are not outlets must still
% exert hydrostatic wall pressure on momentum. Otherwise a lake-at-rest or
% dam-break reservoir next to a wall is not balanced.
%
% These masks identify exterior faces of active cells:
%   wall_left  : west/left face is a reflective wall
%   wall_right : east/right face is a reflective wall
%   wall_up    : north/up face is a reflective wall
%   wall_down  : south/down face is a reflective wall
%
% Outlet faces from Outlet_Properties are excluded from the reflective wall
% masks; they are handled later by the outlet sink logic.
% -------------------------------------------------------------------------
[wall_left, wall_right, wall_up, wall_down] = ...
    build_reflective_wall_face_masks(active, Outlet_Properties);

%% ------------------------------------------------------------------------
% Manning roughness squared
% -------------------------------------------------------------------------
% Needed before the flux calculation because the source/friction predictor
% uses n^2 to relax the predicted momentum before HLL fluxes are computed.
% -------------------------------------------------------------------------

if isempty(roughness_squared)
    n2 = roughness_cell.^2;
else
    n2 = roughness_squared;
end

if isscalar(n2)
    n2 = n2 * ones(size(h), 'like', h);
end

n2(~isfinite(n2) | n2 < 0) = 0;

%% ------------------------------------------------------------------------
% 1. Recover full-momentum state hu/hv from outflow memory
% -------------------------------------------------------------------------
[hu,hv] = recover_momentum_from_outflow(outflow, C_a, dx, h, h_dry, active);

% Remove momentum from dry/inactive cells and prevent unrealistic velocities.
[hu,hv] = clean_and_limit_momentum(hu, hv, h, h_dry, active, max_velocity);

%% ------------------------------------------------------------------------
% 1B. Bed-slope + exact implicit Manning friction predictor
% -------------------------------------------------------------------------
% Current problem:
%   If HLL fluxes are computed using only old hu/hv, the rainfall-added
%   water has to wait for source terms to create momentum. This delays the
%   rising limb for normal Manning n.
%
% New predictor:
%   1) accelerate momentum with bed slope:
%        hu_star = hu - dt*g*h*dz/dx
%        hv_star = hv - dt*g*h*dz/dy
%
%   2) relax hu_star/hv_star with exact pointwise implicit Manning friction.
%
%   3) use hu_flux/hv_flux to compute HLL fluxes in this same time step.
% -------------------------------------------------------------------------

hu_flux = hu;
hv_flux = hv;

use_source_friction_predictor = false;

if flag_source_friction_predictor == 1

    [dzdx,dzdy] = compute_centered_bed_gradients(z_dem, active, dx);

    hu_star = hu - dt .* g .* h .* dzdx;
    hv_star = hv - dt .* g .* h .* dzdy;

    hu_star(~active) = 0;
    hv_star(~active) = 0;

    [hu_flux,hv_flux] = apply_manning_friction_exact_implicit( ...
        hu_star, hv_star, h, n2, dt, g, h_dry, active);

    [hu_flux,hv_flux] = clean_and_limit_momentum( ...
        hu_flux, hv_flux, h, h_dry, active, max_velocity);

    use_source_friction_predictor = true;

end

%% ------------------------------------------------------------------------
% 2. HLL full-momentum fluxes
% -------------------------------------------------------------------------
% flag_HR_full_momentum = 1:
%   hydrostatic reconstruction + well-balanced topographic source treatment
%
% flag_HR_full_momentum = 0:
%   raw HLL fluxes + explicit bed-slope source term added later
% -------------------------------------------------------------------------

if flag_HR_full_momentum == 1

    [Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
     Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N,Hf_x,Hf_y] = ...
        compute_hr_hll_fluxes(h, hu_flux, hv_flux, z_dem, active, g, h_dry);

    apply_explicit_bed_source = false;

else

    [Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
     Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N,Hf_x,Hf_y] = ...
        compute_raw_hll_fluxes_no_hr(h, hu_flux, hv_flux, active, g, h_dry);

    % If the source/friction predictor is active, do NOT add the bed-slope
    % source again after the flux update.
    apply_explicit_bed_source = ~use_source_friction_predictor;

end

%% ------------------------------------------------------------------------
% 3. Positivity-preserving draining limiter
% -------------------------------------------------------------------------
[Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
 Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N] = ...
    apply_draining_limiter(Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
                           Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N, ...
                           h, C_a, dx, dt, h_dry, active);

%% ------------------------------------------------------------------------
% 4. Conservative depth update from face mass fluxes
% -------------------------------------------------------------------------
% Fx_h(:,j) is positive from cell (i,j) to (i,j+1).
% Fy_h(i,:) is positive from cell (i,j) to (i-1,j), i.e., upward/northward.
Vol_Flux = dt * dx * ( ...
    [zeros(ny,1,'like',h), Fx_h(:,1:end-1)] - Fx_h ...
    - Fy_h + [Fy_h(2:end,:); zeros(1,nx,'like',h)] );

%% ------------------------------------------------------------------------
% 5. Conservative momentum update from corrected momentum fluxes
% -------------------------------------------------------------------------
% If the pre-flux source/friction predictor is active, the conservative
% state used for this step starts from hu_flux/hv_flux, not old hu/hv.
% -------------------------------------------------------------------------

if use_source_friction_predictor
    hu_new = hu_flux;
    hv_new = hv_flux;
else
    hu_new = hu;
    hv_new = hv;
end

% X faces, j = 1:nx-1
if nx > 1
    areaL = C_a(:,1:end-1);
    areaR = C_a(:,2:end);

    hu_new(:,1:end-1) = hu_new(:,1:end-1) - dt*dx .* Fx_hu_L(:,1:end-1) ./ areaL;
    hu_new(:,2:end)   = hu_new(:,2:end)   + dt*dx .* Fx_hu_R(:,1:end-1) ./ areaR;

    hv_new(:,1:end-1) = hv_new(:,1:end-1) - dt*dx .* Fx_hv_L(:,1:end-1) ./ areaL;
    hv_new(:,2:end)   = hv_new(:,2:end)   + dt*dx .* Fx_hv_R(:,1:end-1) ./ areaR;
end

% Y faces, i = 2:ny. S = southern/current cell, N = northern/upstream row.
if ny > 1
    areaS = C_a(2:end,:);
    areaN = C_a(1:end-1,:);

    hu_new(2:end,:)   = hu_new(2:end,:)   - dt*dx .* Fy_hu_S(2:end,:) ./ areaS;
    hu_new(1:end-1,:) = hu_new(1:end-1,:) + dt*dx .* Fy_hu_N(2:end,:) ./ areaN;

    hv_new(2:end,:)   = hv_new(2:end,:)   - dt*dx .* Fy_hv_S(2:end,:) ./ areaS;
    hv_new(1:end-1,:) = hv_new(1:end-1,:) + dt*dx .* Fy_hv_N(2:end,:) ./ areaN;
end

%% ------------------------------------------------------------------------
% 5B. Reflective wall pressure fluxes for exterior/inactive faces
% -------------------------------------------------------------------------
% Internal HLL fluxes are only computed where both neighboring cells are
% active. At a closed wall, mass flux is zero, but the hydrostatic pressure
% flux is NOT zero. This correction supplies the missing wall pressure term:
%
%   left wall  -> pushes hu positive/right
%   right wall -> pushes hu negative/left
%   up wall    -> pushes hv negative/down
%   down wall  -> pushes hv positive/up
%
% This is essential for lake-at-rest balance and Ritter/dam-break tests.
% -------------------------------------------------------------------------
[hu_new,hv_new] = apply_reflective_wall_pressure_fluxes( ...
    hu_new, hv_new, h, C_a, dt, dx, g, h_dry, active, ...
    wall_left, wall_right, wall_up, wall_down);

% Depth after conservative intercell flux update.
if use_lookup_subgrid
    V_subgrid_after_flux = max(V_subgrid_old + Vol_Flux, 0);
    h_new = hp2d_subgrid_inverse_volume( ...
        SubgridTables.volume_cell, V_subgrid_after_flux, ...
        SubgridTables.dz, SubgridTables.maxDepth, dx);
    h_new(~active) = 0;
    h_new = max(h_new, 0);
else
    h_new = h + Vol_Flux ./ C_a;
    h_new(~active) = 0;
    h_new = max(h_new, 0);
end

% If a cell was drained strongly, scale momentum consistently.
scale_h = ones(size(h), 'like', h);
wet_old = h > h_dry;
scale_h(wet_old) = min(1, h_new(wet_old) ./ max(h(wet_old), h_dry));
hu_new = hu_new .* scale_h;
hv_new = hv_new .* scale_h;

%% ------------------------------------------------------------------------
% Explicit bed-slope source for non-hydrostatic-reconstruction test
% -------------------------------------------------------------------------
% Momentum source:
%   dhu/dt = -g h dz/dx
%   dhv/dt = -g h dz/dy
%
% Here hv is positive north/up, so dz/dy is computed with y positive upward.
% -------------------------------------------------------------------------

if apply_explicit_bed_source

    [dzdx,dzdy] = compute_centered_bed_gradients(z_dem, active, dx);

    hu_new = hu_new - dt .* g .* h_new .* dzdx;
    hv_new = hv_new - dt .* g .* h_new .* dzdy;

    hu_new(~active) = 0;
    hv_new(~active) = 0;

end

%% ------------------------------------------------------------------------
% 6. Optional post-flux exact implicit Manning friction corrector
% -------------------------------------------------------------------------
% The main friction treatment is now done BEFORE the HLL flux computation.
% Therefore, by default we do not apply a second full friction damping here,
% otherwise the rising limb can become over-damped.
% -------------------------------------------------------------------------

if flag_post_flux_friction_corrector == 1

    dt_corr = post_flux_friction_fraction * dt;

    [hu_new,hv_new] = apply_manning_friction_exact_implicit( ...
        hu_new, hv_new, h_new, n2, dt_corr, g, h_dry, active);

end

[hu_new,hv_new] = clean_and_limit_momentum(hu_new, hv_new, h_new, h_dry, active, max_velocity);
%% ------------------------------------------------------------------------
% 7. Reservoir boundary condition copied in structure from local inertial
% -------------------------------------------------------------------------
% Convert preliminary depth to mm using intercell volume first.
if use_lookup_subgrid
    d_t = 1000 * h_new;
else
    d_t = d_tot + 1000 * Vol_Flux ./ C_a;
end
d_t(~active) = NaN;
d_t(active) = max(d_t(active), 0);

if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = d_tot(reservoir_y(ii),reservoir_x(ii))./1000;
            dt_h = time_step/60;

            available_volume = 1000*(max(dtsup - h1(ii),0))/dt_h;
            area_res = C_a(reservoir_y(ii), reservoir_x(ii));
            dh = min(k1(ii)*(max(dtsup - h1(ii),0))^k2(ii)/area_res*1000*3600,available_volume)*dt_h;

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000*area_res;
            dtsup = dtsup - dh/1000;
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
        else
            dh = 0;
        end

        if ~isnan(yds2(ii))
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;

            available_volume = 1000*(max(dtsup - h2(ii),0))/dt_h;
            area_res = C_a(reservoir_y(ii), reservoir_x(ii));
            dh = min(k3(ii)*(max(dtsup - h2(ii),0))^k4(ii)/area_res*1000*3600,available_volume)*dt_h;

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000*area_res;
            d_tot(yds2(ii),xds2(ii)) = d_tot(yds2(ii),xds2(ii)) + dh;
        end
    end

    d_t = d_tot + 1000 * Vol_Flux ./ C_a;
    d_t(~active) = NaN;
    d_t(active) = max(d_t(active), 0);
    h_new = max(d_t/1000, 0);
    h_new(~active) = 0;
end

V_subgrid_before_outlet = [];
if use_lookup_subgrid
    V_subgrid_before_outlet = V_subgrid_after_flux;
    C_a = hp2d_subgrid_lookup_depth( ...
        SubgridTables.area_cell, h_new, SubgridTables.dz, SubgridTables.maxDepth);
    C_a(~isfinite(C_a) | C_a <= 0) = dx^2;
end

%% ------------------------------------------------------------------------
% 8. Convert face unit fluxes to HydroPol2D outflow storage [mm/h]
% -------------------------------------------------------------------------
outflow_rate_x = Fx_h .* dx;      % [m^3/s] through vertical face
outflow_rate_y = Fy_h .* dx;      % [m^3/s] through horizontal face

outflow_mmh_x = outflow_rate_x ./ C_a * 1000 * 3600;
outflow_mmh_y = outflow_rate_y ./ C_a * 1000 * 3600;

outflow_mmh_x(~isfinite(outflow_mmh_x)) = 0;
outflow_mmh_y(~isfinite(outflow_mmh_y)) = 0;
outflow_mmh_x(~active) = 0;
outflow_mmh_y(~active) = 0;

%% ------------------------------------------------------------------------
% 9. Face-filtered outlet sink: normal-flow or critical-flow
% -------------------------------------------------------------------------
% Outlet_Properties converts outlet cells into outlet FACES.
% Only selected faces are allowed to drain.
%
% This avoids the corner problem:
%   - a corner cell can touch two exterior faces
%   - but only the face selected by steepest outward descent is opened
% -------------------------------------------------------------------------

outlet_flow = zeros(size(d_t), 'like', d_t);  % [mm/h]

Hf = zeros(ny,nx,3,'like',h);
Hf(:,:,1) = Hf_x;
Hf(:,:,2) = Hf_y;

use_face_based_outlet = false;

if ~isempty(Outlet_Properties) && isstruct(Outlet_Properties)

    required_fields = {'face_right','face_left','face_up','face_down'};

    has_all_fields = true;
    for ff = 1:numel(required_fields)
        has_all_fields = has_all_fields && isfield(Outlet_Properties, required_fields{ff});
    end

    if has_all_fields
        face_right = logical(Outlet_Properties.face_right);
        face_left  = logical(Outlet_Properties.face_left);
        face_up    = logical(Outlet_Properties.face_up);
        face_down  = logical(Outlet_Properties.face_down);

        % GPU compatibility
        if isa(h,'gpuArray')
            face_right = gpuArray(face_right);
            face_left  = gpuArray(face_left);
            face_up    = gpuArray(face_up);
            face_down  = gpuArray(face_down);
        end

        face_right = face_right & active;
        face_left  = face_left  & active;
        face_up    = face_up    & active;
        face_down  = face_down  & active;

        use_face_based_outlet = any(face_right(:)) || any(face_left(:)) || ...
                                any(face_up(:))    || any(face_down(:));
    end
end

if use_face_based_outlet

    % Right/east outlet faces
    [d_t,outlet_flow,Hf,hu_new,hv_new] = apply_outlet_face_sink( ...
        face_right, d_t, outlet_flow, Hf, hu_new, hv_new, ...
        C_a, n2, roughness_cell, outlet_type, slope_outlet, ...
        time_step, dx, g, h_min);

    % Left/west outlet faces
    [d_t,outlet_flow,Hf,hu_new,hv_new] = apply_outlet_face_sink( ...
        face_left, d_t, outlet_flow, Hf, hu_new, hv_new, ...
        C_a, n2, roughness_cell, outlet_type, slope_outlet, ...
        time_step, dx, g, h_min);

    % Up/north outlet faces
    [d_t,outlet_flow,Hf,hu_new,hv_new] = apply_outlet_face_sink( ...
        face_up, d_t, outlet_flow, Hf, hu_new, hv_new, ...
        C_a, n2, roughness_cell, outlet_type, slope_outlet, ...
        time_step, dx, g, h_min);

    % Down/south outlet faces
    [d_t,outlet_flow,Hf,hu_new,hv_new] = apply_outlet_face_sink( ...
        face_down, d_t, outlet_flow, Hf, hu_new, hv_new, ...
        C_a, n2, roughness_cell, outlet_type, slope_outlet, ...
        time_step, dx, g, h_min);

else

    % ---------------------------------------------------------------------
    % Fallback: old cell-based outlet behavior
    % ---------------------------------------------------------------------
    % This keeps the function backward-compatible if Outlet_Properties is
    % empty or was not passed.
    % ---------------------------------------------------------------------
    if ~isempty(row_outlet)
        outlet_sub = sub2ind(size(d_t), row_outlet(:), col_outlet(:));

        h_out = max(d_t(outlet_sub), 0) / 1000;  % [m]
        h_out(~isfinite(h_out)) = 0;

        width_out = dx * ones(size(h_out), 'like', h_out);

        if isscalar(C_a)
            area_out = C_a * ones(size(h_out), 'like', h_out);
        else
            area_out = C_a(outlet_sub);
        end
        area_out(~isfinite(area_out) | area_out <= 0) = dx^2;

        if outlet_type == 1
            if isscalar(slope_outlet)
                sqrtS_out = sqrt(abs(slope_outlet)) * ones(size(h_out), 'like', h_out);
            elseif isequal(size(slope_outlet), size(d_t))
                sqrtS_out = sqrt(abs(slope_outlet(outlet_sub)));
            else
                sqrtS_out = sqrt(abs(slope_outlet(1))) * ones(size(h_out), 'like', h_out);
            end
        else
            h_safe = max(h_out, 1e-12);
            sqrtS_out = sqrt(g .* n2(outlet_sub) .* h_safe.^(-1/3));
            sqrtS_out(~isfinite(sqrtS_out)) = 0;
        end

        Hf3 = Hf(:,:,3);
        Hf3(outlet_sub) = h_out;
        Hf(:,:,3) = Hf3;

        if isscalar(roughness_cell)
            n_out = roughness_cell * ones(size(h_out), 'like', h_out);
        else
            n_out = roughness_cell(outlet_sub);
        end

        n_fallback = sqrt(max(n2(outlet_sub), 1e-12));
        bad_n = ~isfinite(n_out) | n_out <= 0;
        n_out(bad_n) = n_fallback(bad_n);

        q_out = (1 ./ n_out) .* ...
                width_out .* ...
                h_out.^(5/3) .* ...
                sqrtS_out ./ ...
                area_out * 1000 * 3600;   % [mm/h]

        qmax_out = max((d_t(outlet_sub) - 1000*h_min), 0) / (time_step/60);
        q_out = min(q_out, qmax_out);
        q_out(~isfinite(q_out)) = 0;

        outlet_flow(outlet_sub) = q_out;

        h_before_out = max(d_t(outlet_sub)/1000, 0);
        d_t(outlet_sub) = d_t(outlet_sub) - q_out * (time_step/60);
        d_t(outlet_sub) = max(d_t(outlet_sub), 1000*h_min);
        h_after_out = max(d_t(outlet_sub)/1000, 0);

        mom_scale = ones(size(h_before_out), 'like', h_before_out);
        idx_scale = h_before_out > 1e-12;
        mom_scale(idx_scale) = h_after_out(idx_scale) ./ h_before_out(idx_scale);

        hu_new(outlet_sub) = hu_new(outlet_sub) .* mom_scale;
        hv_new(outlet_sub) = hv_new(outlet_sub) .* mom_scale;
    end

end

if use_lookup_subgrid
    outlet_volume = outlet_flow .* (time_step/60) ./ 1000 .* C_a;
    outlet_volume(~isfinite(outlet_volume)) = 0;
    V_subgrid_after_outlet = max(V_subgrid_before_outlet - outlet_volume, 0);
    d_t = 1000 * hp2d_subgrid_inverse_volume( ...
        SubgridTables.volume_cell, V_subgrid_after_outlet, ...
        SubgridTables.dz, SubgridTables.maxDepth, dx);
    d_t(~active) = NaN;
    d_t(active) = max(d_t(active), 0);
end

% Final dry cleanup after outlet.
h_final = max(d_t/1000, 0);
h_final(~active) = 0;
[hu_new,hv_new] = clean_and_limit_momentum(hu_new, hv_new, h_final, h_dry, active, max_velocity);

%% ------------------------------------------------------------------------
% 10. Assemble outputs in the same style as Local_Inertial_Model_D4
% -------------------------------------------------------------------------
outflow = zeros(ny,nx,5,'like',h);
outflow(:,:,1) = outflow_mmh_x;
outflow(:,:,2) = outflow_mmh_y;
outflow(:,:,3) = outlet_flow;
outflow(:,:,4) = hu_new;        % [m^2/s], full-momentum memory
outflow(:,:,5) = hv_new;        % [m^2/s], full-momentum memory

matrix_store = outflow(:,:,1:3);

qout_left  = -[zeros(ny,1,'like',h), matrix_store(:,1:end-1,1)];
qout_right =  matrix_store(:,:,1);
qout_up    =  matrix_store(:,:,2);
qout_down  = -[matrix_store(2:end,:,2); zeros(1,nx,'like',h)];

% Total absolute flux leaving/entering each cell, same diagnostic structure
% as the local-inertial model. Uses the first three mm/h flux pages.
mask = matrix_store;
I_tot_end_cell = abs(sum(mask,3)) * dt/1000 * 1/3600 .* C_a; % [m3]
I_tot_end_cell(~active) = NaN;

end

%% ========================================================================
% Helper functions
% ========================================================================
function [hu,hv] = recover_momentum_from_outflow(outflow, C_a, dx, h, h_dry, active)

hu = zeros(size(h), 'like', h);
hv = zeros(size(h), 'like', h);

if isempty(outflow)
    return;
end

if ndims(outflow) >= 3 && size(outflow,3) >= 5
    % Preferred full-momentum memory from previous call.
    hu = outflow(:,:,4);
    hv = outflow(:,:,5);
else
    % First call fallback: infer a crude cell-centered momentum from the
    % previous HydroPol2D face-flux pages [mm/h]. This is only for startup.
    if ndims(outflow) >= 3 && size(outflow,3) >= 1
        hu = outflow(:,:,1) ./ 1000 ./ 3600 .* C_a ./ dx;
    end
    if ndims(outflow) >= 3 && size(outflow,3) >= 2
        hv = outflow(:,:,2) ./ 1000 ./ 3600 .* C_a ./ dx;
    end
end

hu(~isfinite(hu)) = 0;
hv(~isfinite(hv)) = 0;
hu(~active | h <= h_dry) = 0;
hv(~active | h <= h_dry) = 0;

end

function [hu,hv] = clean_and_limit_momentum(hu, hv, h, h_dry, active, max_velocity)

hu(~isfinite(hu)) = 0;
hv(~isfinite(hv)) = 0;

iswet = active & h > h_dry;
hu(~iswet) = 0;
hv(~iswet) = 0;

speed = zeros(size(h), 'like', h);
speed(iswet) = sqrt(hu(iswet).^2 + hv(iswet).^2) ./ max(h(iswet), h_dry);

idx = iswet & speed > max_velocity;
if any(idx(:))
    fac = max_velocity ./ max(speed(idx), 1e-12);
    hu(idx) = hu(idx) .* fac;
    hv(idx) = hv(idx) .* fac;
end

end

function [Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
          Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N,Hf_x,Hf_y] = ...
    compute_hr_hll_fluxes(h, hu, hv, z, active, g, h_dry)

[ny,nx] = size(h);

Fx_h    = zeros(ny,nx,'like',h);
Fx_hu_L = zeros(ny,nx,'like',h);
Fx_hu_R = zeros(ny,nx,'like',h);
Fx_hv_L = zeros(ny,nx,'like',h);
Fx_hv_R = zeros(ny,nx,'like',h);

Fy_h    = zeros(ny,nx,'like',h);
Fy_hu_S = zeros(ny,nx,'like',h);
Fy_hu_N = zeros(ny,nx,'like',h);
Fy_hv_S = zeros(ny,nx,'like',h);
Fy_hv_N = zeros(ny,nx,'like',h);

Hf_x = zeros(ny,nx,'like',h);
Hf_y = zeros(ny,nx,'like',h);

eta = z + h;

% ----------------------- X faces: left -> right -------------------------
if nx > 1
    face_active = active(:,1:end-1) & active(:,2:end);

    hL = h(:,1:end-1);    hR = h(:,2:end);
    zL = z(:,1:end-1);    zR = z(:,2:end);
    etaL = eta(:,1:end-1); etaR = eta(:,2:end);

    zf = max(zL, zR);
    hLs = max(etaL - zf, 0);
    hRs = max(etaR - zf, 0);

    uL = zeros(size(hL), 'like', h); vL = uL;
    uR = zeros(size(hR), 'like', h); vR = uR;
    wetL = hL > h_dry;
    wetR = hR > h_dry;
    huL = hu(:,1:end-1); hvL = hv(:,1:end-1);
    huR = hu(:,2:end);   hvR = hv(:,2:end);
    uL(wetL) = huL(wetL) ./ hL(wetL);
    vL(wetL) = hvL(wetL) ./ hL(wetL);
    uR(wetR) = huR(wetR) ./ hR(wetR);
    vR(wetR) = hvR(wetR) ./ hR(wetR);

    qnL = hLs .* uL; qtL = hLs .* vL;
    qnR = hRs .* uR; qtR = hRs .* vR;

    [Fh,Fqn,Fqt] = hll_normal_flux(hLs, qnL, qtL, hRs, qnR, qtR, g, h_dry);

    corrL = 0.5*g*(hL.^2 - hLs.^2);
    corrR = 0.5*g*(hR.^2 - hRs.^2);

    Fh(~face_active)   = 0;
    Fqn(~face_active)  = 0;
    Fqt(~face_active)  = 0;
    corrL(~face_active)= 0;
    corrR(~face_active)= 0;

    Fx_h(:,1:end-1)    = Fh;
    Fx_hu_L(:,1:end-1) = Fqn + corrL;
    Fx_hu_R(:,1:end-1) = Fqn + corrR;
    Fx_hv_L(:,1:end-1) = Fqt;
    Fx_hv_R(:,1:end-1) = Fqt;

    Hf_x(:,1:end-1) = max(hLs, hRs);
end

% ----------------------- Y faces: south/current -> north/up -------------
if ny > 1
    face_active = active(2:end,:) & active(1:end-1,:);

    hS = h(2:end,:);      hN = h(1:end-1,:);
    zS = z(2:end,:);      zN = z(1:end-1,:);
    etaS = eta(2:end,:);  etaN = eta(1:end-1,:);

    zf = max(zS, zN);
    hSs = max(etaS - zf, 0);
    hNs = max(etaN - zf, 0);

    uS = zeros(size(hS), 'like', h); vS = uS;
    uN = zeros(size(hN), 'like', h); vN = uN;
    wetS = hS > h_dry;
    wetN = hN > h_dry;

    % hu is east/right. hv is north/up.
    huS = hu(2:end,:);   hvS = hv(2:end,:);
    huN = hu(1:end-1,:); hvN = hv(1:end-1,:);
    uS(wetS) = huS(wetS) ./ hS(wetS);
    vS(wetS) = hvS(wetS) ./ hS(wetS);
    uN(wetN) = huN(wetN) ./ hN(wetN);
    vN(wetN) = hvN(wetN) ./ hN(wetN);

    qnS = hSs .* vS; qtS = hSs .* uS;
    qnN = hNs .* vN; qtN = hNs .* uN;

    [Fh,Fqn,Fqt] = hll_normal_flux(hSs, qnS, qtS, hNs, qnN, qtN, g, h_dry);

    corrS = 0.5*g*(hS.^2 - hSs.^2);
    corrN = 0.5*g*(hN.^2 - hNs.^2);

    Fh(~face_active)    = 0;
    Fqn(~face_active)   = 0;
    Fqt(~face_active)   = 0;
    corrS(~face_active) = 0;
    corrN(~face_active) = 0;

    Fy_h(2:end,:)    = Fh;
    Fy_hv_S(2:end,:) = Fqn + corrS;  % normal momentum = hv
    Fy_hv_N(2:end,:) = Fqn + corrN;
    Fy_hu_S(2:end,:) = Fqt;          % transverse momentum = hu
    Fy_hu_N(2:end,:) = Fqt;

    Hf_y(2:end,:) = max(hSs, hNs);
end

end

function [Fh,Fqn,Fqt] = hll_normal_flux(hL, qnL, qtL, hR, qnR, qtR, g, h_dry)

Fh  = zeros(size(hL), 'like', hL);
Fqn = zeros(size(hL), 'like', hL);
Fqt = zeros(size(hL), 'like', hL);

uL = zeros(size(hL), 'like', hL);
uR = zeros(size(hR), 'like', hR);
vtL = zeros(size(hL), 'like', hL);
vtR = zeros(size(hR), 'like', hR);

wetL = hL > h_dry;
wetR = hR > h_dry;

uL(wetL)  = qnL(wetL) ./ hL(wetL);
vtL(wetL) = qtL(wetL) ./ hL(wetL);
uR(wetR)  = qnR(wetR) ./ hR(wetR);
vtR(wetR) = qtR(wetR) ./ hR(wetR);

cL = sqrt(g .* max(hL,0));
cR = sqrt(g .* max(hR,0));

FL_h  = qnL;
FL_qn = qnL .* uL + 0.5*g*hL.^2;
FL_qt = qtL .* uL;

FR_h  = qnR;
FR_qn = qnR .* uR + 0.5*g*hR.^2;
FR_qt = qtR .* uR;

sL = min(uL - cL, uR - cR);
sR = max(uL + cL, uR + cR);

idxL = sL >= 0;
idxR = sR <= 0;
idxH = ~(idxL | idxR) & (sR > sL);

Fh(idxL)  = FL_h(idxL);
Fqn(idxL) = FL_qn(idxL);
Fqt(idxL) = FL_qt(idxL);

Fh(idxR)  = FR_h(idxR);
Fqn(idxR) = FR_qn(idxR);
Fqt(idxR) = FR_qt(idxR);

den = sR - sL;
Fh(idxH) = (sR(idxH).*FL_h(idxH) - sL(idxH).*FR_h(idxH) + ...
             sL(idxH).*sR(idxH).*(hR(idxH) - hL(idxH))) ./ den(idxH);
Fqn(idxH) = (sR(idxH).*FL_qn(idxH) - sL(idxH).*FR_qn(idxH) + ...
              sL(idxH).*sR(idxH).*(qnR(idxH) - qnL(idxH))) ./ den(idxH);
Fqt(idxH) = (sR(idxH).*FL_qt(idxH) - sL(idxH).*FR_qt(idxH) + ...
              sL(idxH).*sR(idxH).*(qtR(idxH) - qtL(idxH))) ./ den(idxH);

% Both sides dry or numerically invalid.
dry = (~wetL & ~wetR) | ~isfinite(Fh) | ~isfinite(Fqn) | ~isfinite(Fqt);
Fh(dry)  = 0;
Fqn(dry) = 0;
Fqt(dry) = 0;

end

function [Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
          Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N] = ...
    apply_draining_limiter(Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
                           Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N, ...
                           h, C_a, dx, dt, h_dry, active)

[ny,nx] = size(h);
outV = zeros(size(h), 'like', h);

if nx > 1
    f = Fx_h(:,1:end-1);
    outV(:,1:end-1) = outV(:,1:end-1) + dt*dx*max(f,0);
    outV(:,2:end)   = outV(:,2:end)   + dt*dx*max(-f,0);
end

if ny > 1
    f = Fy_h(2:end,:);
    outV(2:end,:)   = outV(2:end,:)   + dt*dx*max(f,0);
    outV(1:end-1,:) = outV(1:end-1,:) + dt*dx*max(-f,0);
end

available = max(h - h_dry, 0) .* C_a;
available(~active) = 0;

theta = ones(size(h), 'like', h);
idx = active & outV > available & outV > 0;
theta(idx) = available(idx) ./ max(outV(idx), 1e-30);
theta(~isfinite(theta)) = 0;
theta = max(min(theta,1),0);

if nx > 1
    f = Fx_h(:,1:end-1);
    scale = theta(:,1:end-1);
    thetaR = theta(:,2:end);
    neg = f < 0;
    scale(neg) = thetaR(neg);

    Fx_h(:,1:end-1)    = Fx_h(:,1:end-1)    .* scale;
    Fx_hu_L(:,1:end-1) = Fx_hu_L(:,1:end-1) .* scale;
    Fx_hu_R(:,1:end-1) = Fx_hu_R(:,1:end-1) .* scale;
    Fx_hv_L(:,1:end-1) = Fx_hv_L(:,1:end-1) .* scale;
    Fx_hv_R(:,1:end-1) = Fx_hv_R(:,1:end-1) .* scale;
end

if ny > 1
    f = Fy_h(2:end,:);
    scale = theta(2:end,:);
    thetaN = theta(1:end-1,:);
    neg = f < 0;
    scale(neg) = thetaN(neg);

    Fy_h(2:end,:)    = Fy_h(2:end,:)    .* scale;
    Fy_hu_S(2:end,:) = Fy_hu_S(2:end,:) .* scale;
    Fy_hu_N(2:end,:) = Fy_hu_N(2:end,:) .* scale;
    Fy_hv_S(2:end,:) = Fy_hv_S(2:end,:) .* scale;
    Fy_hv_N(2:end,:) = Fy_hv_N(2:end,:) .* scale;
end

end

function [hu,hv] = apply_manning_friction_semiimplicit(hu, hv, h, n2, dt, g, h_dry, active)

iswet = active & h > h_dry;

speed_q = sqrt(hu.^2 + hv.^2);  % [m^2/s]

den = ones(size(h), 'like', h);
den(iswet) = 1 + dt .* g .* n2(iswet) .* speed_q(iswet) ./ max(h(iswet),h_dry).^(7/3);
den(~isfinite(den) | den <= 0) = 1;

hu(iswet) = hu(iswet) ./ den(iswet);
hv(iswet) = hv(iswet) ./ den(iswet);

hu(~iswet) = 0;
hv(~iswet) = 0;

end

function [d_t,outlet_flow,Hf,hu_new,hv_new] = apply_outlet_face_sink( ...
    face_mask, d_t, outlet_flow, Hf, hu_new, hv_new, ...
    C_a, n2, roughness_cell, outlet_type, slope_outlet, ...
    time_step, dx, g, h_min)
%APPLY_OUTLET_FACE_SINK
% -------------------------------------------------------------------------
% Applies an outlet sink only to cells whose selected exterior face is open.
%
% The face direction was already decided during preprocessing by
% define_outlet_faces_steepest_descent.
%
% This function still removes volume from the owner cell, but only cells
% with an active outlet face are allowed to drain.
% -------------------------------------------------------------------------

if ~any(face_mask(:))
    return;
end

outlet_sub = find(face_mask);

h_out = max(d_t(outlet_sub), 0) / 1000;  % [m]
h_out(~isfinite(h_out)) = 0;

width_out = dx * ones(size(h_out), 'like', h_out);

if isscalar(C_a)
    area_out = C_a * ones(size(h_out), 'like', h_out);
else
    area_out = C_a(outlet_sub);
end
area_out(~isfinite(area_out) | area_out <= 0) = dx^2;

if outlet_type == 1

    % Normal-flow outlet
    if isscalar(slope_outlet)
        sqrtS_out = sqrt(abs(slope_outlet)) * ones(size(h_out), 'like', h_out);
    elseif isequal(size(slope_outlet), size(d_t))
        sqrtS_out = sqrt(abs(slope_outlet(outlet_sub)));
    else
        sqrtS_out = sqrt(abs(slope_outlet(1))) * ones(size(h_out), 'like', h_out);
    end

else

    % Critical-flow outlet, using the same Manning-compatible expression
    % used in the current full-momentum outlet block.
    h_safe = max(h_out, 1e-12);
    sqrtS_out = sqrt(g .* n2(outlet_sub) .* h_safe.^(-1/3));
    sqrtS_out(~isfinite(sqrtS_out)) = 0;

end

Hf3 = Hf(:,:,3);
Hf3(outlet_sub) = h_out;
Hf(:,:,3) = Hf3;

if isscalar(roughness_cell)
    n_out = roughness_cell * ones(size(h_out), 'like', h_out);
else
    n_out = roughness_cell(outlet_sub);
end

n_fallback = sqrt(max(n2(outlet_sub), 1e-12));
bad_n = ~isfinite(n_out) | n_out <= 0;
n_out(bad_n) = n_fallback(bad_n);

q_out = (1 ./ n_out) .* ...
        width_out .* ...
        h_out.^(5/3) .* ...
        sqrtS_out ./ ...
        area_out * 1000 * 3600;   % [mm/h removed from outlet cell]

% Do not apply h_dry here.
% The outlet is allowed to drain available water down to h_min.
qmax_out = max((d_t(outlet_sub) - 1000*h_min), 0) / (time_step/60);
q_out = min(q_out, qmax_out);
q_out(~isfinite(q_out)) = 0;

% Multiple calls are possible, but Option A should select only one face per
% cell. The addition keeps this safe if later we allow diagonal/corner split.
outlet_flow(outlet_sub) = outlet_flow(outlet_sub) + q_out;

h_before_out = max(d_t(outlet_sub)/1000, 0);

d_t(outlet_sub) = d_t(outlet_sub) - q_out * (time_step/60);
d_t(outlet_sub) = max(d_t(outlet_sub), 1000*h_min);

h_after_out = max(d_t(outlet_sub)/1000, 0);

% Remove proportional momentum from drained cells.
mom_scale = ones(size(h_before_out), 'like', h_before_out);
idx_scale = h_before_out > 1e-12;
mom_scale(idx_scale) = h_after_out(idx_scale) ./ h_before_out(idx_scale);

hu_new(outlet_sub) = hu_new(outlet_sub) .* mom_scale;
hv_new(outlet_sub) = hv_new(outlet_sub) .* mom_scale;

end

function [Fx_h,Fx_hu_L,Fx_hu_R,Fx_hv_L,Fx_hv_R, ...
          Fy_h,Fy_hu_S,Fy_hu_N,Fy_hv_S,Fy_hv_N,Hf_x,Hf_y] = ...
    compute_raw_hll_fluxes_no_hr(h, hu, hv, active, g, h_dry)
%COMPUTE_RAW_HLL_FLUXES_NO_HR
% -------------------------------------------------------------------------
% Raw HLL shallow-water fluxes without hydrostatic reconstruction.
%
% This is useful for simple plane-slope tests where we want to diagnose
% whether hydrostatic reconstruction is over-diffusing or suppressing
% acceleration.
%
% Bed slope is NOT included here. It must be added explicitly later through:
%   dhu/dt = -g h dz/dx
%   dhv/dt = -g h dz/dy
% -------------------------------------------------------------------------

[ny,nx] = size(h);

Fx_h    = zeros(ny,nx,'like',h);
Fx_hu_L = zeros(ny,nx,'like',h);
Fx_hu_R = zeros(ny,nx,'like',h);
Fx_hv_L = zeros(ny,nx,'like',h);
Fx_hv_R = zeros(ny,nx,'like',h);

Fy_h    = zeros(ny,nx,'like',h);
Fy_hu_S = zeros(ny,nx,'like',h);
Fy_hu_N = zeros(ny,nx,'like',h);
Fy_hv_S = zeros(ny,nx,'like',h);
Fy_hv_N = zeros(ny,nx,'like',h);

Hf_x = zeros(ny,nx,'like',h);
Hf_y = zeros(ny,nx,'like',h);

%% X faces: left cell -> right cell
if nx > 1

    face_active = active(:,1:end-1) & active(:,2:end);

    hL  = h(:,1:end-1);
    hR  = h(:,2:end);

    huL = hu(:,1:end-1);
    hvL = hv(:,1:end-1);

    huR = hu(:,2:end);
    hvR = hv(:,2:end);

    qnL = huL;   % normal discharge in x
    qtL = hvL;   % transverse discharge

    qnR = huR;
    qtR = hvR;

    [Fh,Fqn,Fqt] = hll_normal_flux(hL, qnL, qtL, hR, qnR, qtR, g, h_dry);

    Fh(~face_active)  = 0;
    Fqn(~face_active) = 0;
    Fqt(~face_active) = 0;

    Fx_h(:,1:end-1)    = Fh;
    Fx_hu_L(:,1:end-1) = Fqn;
    Fx_hu_R(:,1:end-1) = Fqn;
    Fx_hv_L(:,1:end-1) = Fqt;
    Fx_hv_R(:,1:end-1) = Fqt;

    Hf_x(:,1:end-1) = max(hL,hR);

end

%% Y faces: south/current cell -> north/up cell
if ny > 1

    face_active = active(2:end,:) & active(1:end-1,:);

    hS  = h(2:end,:);
    hN  = h(1:end-1,:);

    huS = hu(2:end,:);
    hvS = hv(2:end,:);

    huN = hu(1:end-1,:);
    hvN = hv(1:end-1,:);

    qnS = hvS;   % normal discharge in y/north direction
    qtS = huS;   % transverse discharge

    qnN = hvN;
    qtN = huN;

    [Fh,Fqn,Fqt] = hll_normal_flux(hS, qnS, qtS, hN, qnN, qtN, g, h_dry);

    Fh(~face_active)  = 0;
    Fqn(~face_active) = 0;
    Fqt(~face_active) = 0;

    Fy_h(2:end,:)    = Fh;
    Fy_hv_S(2:end,:) = Fqn;
    Fy_hv_N(2:end,:) = Fqn;
    Fy_hu_S(2:end,:) = Fqt;
    Fy_hu_N(2:end,:) = Fqt;

    Hf_y(2:end,:) = max(hS,hN);

end

end

function [dzdx,dzdy] = compute_centered_bed_gradients(z, active, dx)
%COMPUTE_CENTERED_BED_GRADIENTS
% -------------------------------------------------------------------------
% Computes robust centered/one-sided bed gradients.
%
% x is positive to the right/east.
% y is positive upward/north, which corresponds to decreasing row index.
% -------------------------------------------------------------------------

[ny,nx] = size(z);

dzdx = zeros(ny,nx,'like',z);
dzdy = zeros(ny,nx,'like',z);

%% Neighbor elevations and active masks
activeR = false(ny,nx);
activeL = false(ny,nx);
activeU = false(ny,nx);
activeD = false(ny,nx);

zR = z; zL = z; zU = z; zD = z;
zR(:) = NaN; zL(:) = NaN; zU(:) = NaN; zD(:) = NaN;

if nx > 1
    activeR(:,1:end-1) = active(:,2:end);
    activeL(:,2:end)   = active(:,1:end-1);

    zR(:,1:end-1) = z(:,2:end);
    zL(:,2:end)   = z(:,1:end-1);
end

if ny > 1
    activeU(2:end,:)   = active(1:end-1,:);
    activeD(1:end-1,:) = active(2:end,:);

    zU(2:end,:)   = z(1:end-1,:);
    zD(1:end-1,:) = z(2:end,:);
end

%% dz/dx, x positive right
bothX = active & activeR & activeL;
onlyR = active & activeR & ~activeL;
onlyL = active & activeL & ~activeR;

dzdx(bothX) = (zR(bothX) - zL(bothX)) ./ (2*dx);
dzdx(onlyR) = (zR(onlyR) - z(onlyR)) ./ dx;
dzdx(onlyL) = (z(onlyL)  - zL(onlyL)) ./ dx;

%% dz/dy, y positive upward/north
bothY = active & activeU & activeD;
onlyU = active & activeU & ~activeD;
onlyD = active & activeD & ~activeU;

dzdy(bothY) = (zU(bothY) - zD(bothY)) ./ (2*dx);
dzdy(onlyU) = (zU(onlyU) - z(onlyU)) ./ dx;
dzdy(onlyD) = (z(onlyD)  - zD(onlyD)) ./ dx;

dzdx(~active | ~isfinite(dzdx)) = 0;
dzdy(~active | ~isfinite(dzdy)) = 0;

end

function [hu,hv] = apply_manning_friction_exact_implicit(hu_star, hv_star, h, n2, dt, g, h_dry, active)
%APPLY_MANNING_FRICTION_EXACT_IMPLICIT
% -------------------------------------------------------------------------
% Exact pointwise implicit Manning friction update for shallow-water momentum.
%
% Momentum vector:
%   q = [hu, hv]
%
% Manning friction source:
%   dq/dt = - g n^2 |q| q / h^(7/3)
%
% Fully implicit scalar form:
%   q_new = q_star / (1 + dt*C*|q_new|)
%
% where:
%   C = g n^2 / h^(7/3)
%
% Let:
%   Q = |q_star|
%   r = |q_new|
%   a = dt*C
%
% Then:
%   a*r^2 + r - Q = 0
%
% Stable analytical solution:
%   r = 2Q / (1 + sqrt(1 + 4aQ))
%
% and:
%   q_new = q_star * r/Q
% -------------------------------------------------------------------------

hu = zeros(size(hu_star), 'like', hu_star);
hv = zeros(size(hv_star), 'like', hv_star);

iswet = active & h > h_dry;

Qstar = sqrt(hu_star.^2 + hv_star.^2);
Qstar(~isfinite(Qstar)) = 0;

Cfric = zeros(size(h), 'like', h);
Cfric(iswet) = g .* n2(iswet) ./ max(h(iswet), h_dry).^(7/3);
Cfric(~isfinite(Cfric) | Cfric < 0) = 0;

a = dt .* Cfric;

r = zeros(size(h), 'like', h);

idx = iswet & Qstar > 0;

% Numerically stable quadratic solution:
% r = (-1 + sqrt(1 + 4*a*Qstar)) / (2*a)
% rewritten as:
% r = 2*Qstar / (1 + sqrt(1 + 4*a*Qstar))
r(idx) = 2 .* Qstar(idx) ./ ...
         (1 + sqrt(1 + 4 .* a(idx) .* Qstar(idx)));

scale = zeros(size(h), 'like', h);
scale(idx) = r(idx) ./ max(Qstar(idx), 1e-30);

hu(idx) = hu_star(idx) .* scale(idx);
hv(idx) = hv_star(idx) .* scale(idx);

hu(~iswet) = 0;
hv(~iswet) = 0;

hu(~isfinite(hu)) = 0;
hv(~isfinite(hv)) = 0;

end

function [wall_left, wall_right, wall_up, wall_down] = ...
    build_reflective_wall_face_masks(active, Outlet_Properties)
%BUILD_REFLECTIVE_WALL_FACE_MASKS
% -------------------------------------------------------------------------
% Builds face masks for reflective wall boundary conditions.
%
% A face is reflective if:
%   1) the owner cell is active, and
%   2) the neighboring cell across that face is inactive/outside, and
%   3) the face is NOT marked as an outlet face in Outlet_Properties.
%
% This is face-based, so a corner cell can have one outlet face and one
% reflective wall face.
% -------------------------------------------------------------------------

[ny,nx] = size(active);

false_mask = active;
false_mask(:) = false;

% Neighbor activity masks from the owner-cell perspective.
activeR = false_mask;
activeL = false_mask;
activeU = false_mask;
activeD = false_mask;

if nx > 1
    activeR(:,1:end-1) = active(:,2:end);
    activeL(:,2:end)   = active(:,1:end-1);
end

if ny > 1
    activeU(2:end,:)   = active(1:end-1,:);
    activeD(1:end-1,:) = active(2:end,:);
end

% Exterior faces of active cells.
ext_right = active & ~activeR;
ext_left  = active & ~activeL;
ext_up    = active & ~activeU;
ext_down  = active & ~activeD;

% Outlet face masks, if available.
out_right = false_mask;
out_left  = false_mask;
out_up    = false_mask;
out_down  = false_mask;

if ~isempty(Outlet_Properties) && isstruct(Outlet_Properties)

    if isfield(Outlet_Properties,'face_right') && ~isempty(Outlet_Properties.face_right)
        out_right = logical(Outlet_Properties.face_right);
    end

    if isfield(Outlet_Properties,'face_left') && ~isempty(Outlet_Properties.face_left)
        out_left = logical(Outlet_Properties.face_left);
    end

    if isfield(Outlet_Properties,'face_up') && ~isempty(Outlet_Properties.face_up)
        out_up = logical(Outlet_Properties.face_up);
    end

    if isfield(Outlet_Properties,'face_down') && ~isempty(Outlet_Properties.face_down)
        out_down = logical(Outlet_Properties.face_down);
    end

    % GPU compatibility if active is already on GPU.
    if isa(active,'gpuArray')
        out_right = gpuArray(out_right);
        out_left  = gpuArray(out_left);
        out_up    = gpuArray(out_up);
        out_down  = gpuArray(out_down);
    end

    % Protect against accidental nonmatching outlet masks.
    if ~isequal(size(out_right), size(active)); out_right = false_mask; end
    if ~isequal(size(out_left),  size(active)); out_left  = false_mask; end
    if ~isequal(size(out_up),    size(active)); out_up    = false_mask; end
    if ~isequal(size(out_down),  size(active)); out_down  = false_mask; end
end

% Reflective walls are exterior faces that are not outlet faces.
wall_right = ext_right & ~out_right;
wall_left  = ext_left  & ~out_left;
wall_up    = ext_up    & ~out_up;
wall_down  = ext_down  & ~out_down;

end

function [hu_new,hv_new] = apply_reflective_wall_pressure_fluxes( ...
    hu_new, hv_new, h, C_a, dt, dx, g, h_dry, active, ...
    wall_left, wall_right, wall_up, wall_down)
%APPLY_REFLECTIVE_WALL_PRESSURE_FLUXES
% -------------------------------------------------------------------------
% Adds the missing hydrostatic pressure flux at reflective wall faces.
%
% For a wall face, the normal mass flux is zero, but the normal momentum
% flux is the hydrostatic pressure term:
%
%   p_wall = 1/2 g h^2
%
% Signs follow HydroPol2D's momentum convention:
%   hu > 0 to the right/east
%   hv > 0 upward/north
% -------------------------------------------------------------------------

iswet = active & h > h_dry;

p_wall = 0.5 .* g .* h.^2;      % [m^3/s^2] unit-width momentum flux
p_wall(~iswet) = 0;
p_wall(~isfinite(p_wall)) = 0;

coef = dt .* dx ./ C_a;
coef(~isfinite(coef) | coef < 0) = 0;

% Left/west wall: pressure pushes water to the right/east.
hu_new(wall_left) = hu_new(wall_left) + coef(wall_left) .* p_wall(wall_left);

% Right/east wall: pressure pushes water to the left/west.
hu_new(wall_right) = hu_new(wall_right) - coef(wall_right) .* p_wall(wall_right);

% Up/north wall: pressure pushes water downward/southward, i.e., negative hv.
hv_new(wall_up) = hv_new(wall_up) - coef(wall_up) .* p_wall(wall_up);

% Down/south wall: pressure pushes water upward/northward, i.e., positive hv.
hv_new(wall_down) = hv_new(wall_down) + coef(wall_down) .* p_wall(wall_down);

hu_new(~isfinite(hu_new)) = 0;
hv_new(~isfinite(hv_new)) = 0;

end
