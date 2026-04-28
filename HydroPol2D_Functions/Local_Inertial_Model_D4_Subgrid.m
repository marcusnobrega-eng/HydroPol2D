function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a,eta_t,V_t] = ...
    Local_Inertial_Model_D4_Subgrid( ...
    flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2, ...
    flag_reservoir,z,d_tot,d_p,roughness_cell,roughness_squared,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet, ...
    row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,nc,nf,River_Width,River_Depth, ...
    Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev,Subgrid_Properties,flag_inflow,SubgridTables)
%--------------------------------------------------------------------------
% Local_Inertial_Model_D4_Subgrid
%
% PURPOSE
% -------------------------------------------------------------------------
% Standalone local-inertial solver for the shared-face subgrid
% formulation.
%
% This solver is intended for isolated testing of the new subgrid method.
% It includes:
%   - shared-face subgrid hydraulics
%   - continuity solved in the VOLUME domain
%   - reservoir boundary condition logic
%   - outlet flow logic
%   - same style of outputs used by HydroPol2D
%
% IMPORTANT STATE INTERPRETATION
% -------------------------------------------------------------------------
% In this standalone subgrid solver:
%
%   d_tot, d_p, d_t   -> representative depth above cell invert [mm]
%
% and the absolute water surface elevation is:
%
%   eta = invert_el + d_rep
%
% where:
%   invert_el = minimum fine elevation inside each coarse cell [m]
%
% Therefore, this function does NOT interpret d_tot as depth above the DEM.
%
% CONTINUITY
% -------------------------------------------------------------------------
% The continuity equation is solved in stored volume:
%
%   V^(n+1) = V^n + Vol_Flux
%
% and then the updated representative depth is obtained by inverting the
% preprocessed cell storage curve.
%
% REQUIRED EXTERNAL HELPER
% -------------------------------------------------------------------------
% This function calls the helper:
%
%   subgrid_channel_function(...)
%
%
% NOTES
% -------------------------------------------------------------------------
% 1) This function is subgrid-only. No legacy/non-subgrid branch is present.
% 2) River_Width, River_Depth, Qc_prev, Qf_prev, Qci_prev, Qfi_prev, nf are
%    kept in the signature for future integration compatibility, but are not
%    used in this paper-style shared-face formulation.
% 3) Internally, SI units are used for hydraulic computations:
%       depth  -> m
%       eta    -> m
%       V      -> m3
%       q      -> m2/s
%       Q      -> m3/s
%--------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% 1) DOMAIN AND BASIC SETTINGS
% -------------------------------------------------------------------------
if isgpuarray(cell_area)
    ny = size(z,1);
    nx = size(z,2);
else
    ny = size(z,1);
    nx = size(z,2);
end

dt = time_step * 60;  % [s]
g  = 9.81;

% Minimum operative depth, consistent with the previous D4 logic
if flag_inflow ~= 1
    h_min = 1e-6;     % [m]
else
    h_min = 0;
end

% Outlet cells should not behave as river-channel cells later
River_Width(logical(outlet_index)) = 0;
River_Depth(logical(outlet_index)) = 0;

%% ------------------------------------------------------------------------
% 2) CURRENT CELL STATE: REPRESENTATIVE DEPTH, INVERT, eta
% -------------------------------------------------------------------------
invert_el = Subgrid_Properties.invert_el;   % [m]

% Current representative depth above invert [m]
drep_n = max(d_tot/1000, 0);

% Absolute water surface elevation [m]
eta_n = invert_el + drep_n;

%% ------------------------------------------------------------------------
% 3) CURRENT STORED VOLUME AND WETTED AREA FROM LOOKUP TABLES
% -------------------------------------------------------------------------
V_n = hp2d_lookup_uniform( ...
    SubgridTables.volume_cell, drep_n, ...
    SubgridTables.dz, SubgridTables.maxDepth);

C_a = hp2d_lookup_uniform( ...
    SubgridTables.area_cell, drep_n, ...
    SubgridTables.dz, SubgridTables.maxDepth);

current_volume = nansum(nansum(V_n));

%% ------------------------------------------------------------------------
% 4) PREVIOUS FACE q [m2/s]
% -------------------------------------------------------------------------
% The incoming outflow variable from HydroPol2D convention is in mm/h.
% Convert only the x and y components back to q = Q/W [m2/s].
q_prev = outflow(:,:,1:2) / 1000 / 3600 * Resolution^2 / Resolution;

%% ------------------------------------------------------------------------
% 5) SUBGRID SHARED-FACE LOCAL-INERTIAL UPDATE
% -------------------------------------------------------------------------
% This helper will be built separately next.
%
% Expected outputs:
%   q_face   [ny x nx x 2]   updated face discharge per unit width [m2/s]
%   Hf_x     [ny x nx]       effective x-face depth [m]
%   Hf_y     [ny x nx]       effective y-face depth [m]
%   Wf_x     [ny x nx]       wetted x-face width [m]
%   Wf_y     [ny x nx]       wetted y-face width [m]
%
[q_face,Hf_x,Hf_y,Wf_x,Wf_y] = subgrid_topography_model( ...
    flag_numerical_scheme,eta_n,q_prev,nc,Resolution,dt,h_min,idx_nan,Subgrid_Properties,SubgridTables);

% Pack Hf in the same style as the previous D4 function
Hf = zeros(ny,nx,3,'like',eta_n);
Hf(:,:,1) = Hf_x;
Hf(:,:,2) = Hf_y;

%% ------------------------------------------------------------------------
% 6) OPTIONAL VELOCITY LIMITERS
% -------------------------------------------------------------------------
if flag_critical == 1
    critical_velocity = Hf(:,:,1:2) .* sqrt(g * Hf(:,:,1:2));
    q_face = min(max(q_face, -critical_velocity), critical_velocity);
end

max_velocity = 10;  % [m/s]
threshold_velocity = Hf(:,:,1:2) * max_velocity;
q_face = min(max(q_face, -threshold_velocity), threshold_velocity);

%% ------------------------------------------------------------------------
% 7) CONVERT q -> Q
% -------------------------------------------------------------------------
Q_face = zeros(size(q_face), 'like', q_face);   % [m3/s]
Q_face(:,:,1) = q_face(:,:,1) * Resolution;
Q_face(:,:,2) = q_face(:,:,2) * Resolution;

%% ------------------------------------------------------------------------
% 8) STORE OUTFLOW IN MODEL COMPATIBLE UNITS (mm/h)
% -------------------------------------------------------------------------
outflow = zeros(ny,nx,3,'like',q_face);
outflow(:,:,1:2) = Q_face ./ cell_area * 1000 * 3600;

matrix_store = outflow(:,:,1:2);

%% ------------------------------------------------------------------------
% 9) CONTINUITY UPDATE IN VOLUME DOMAIN
% -------------------------------------------------------------------------
% Sign convention:
%   Q_face(i,j,1) positive from cell (i,j) to cell (i,j+1)
%   Q_face(i,j,2) positive from cell (i,j) to cell (i+1,j)
%
Qwest  = [zeros(ny,1,'like',Q_face(:,:,1)), Q_face(:,1:(nx-1),1)];
Qeast  = Q_face(:,:,1);

Qnorth = [zeros(1,nx,'like',Q_face(:,:,2)); Q_face(1:(ny-1),:,2)];
Qsouth = Q_face(:,:,2);

Vol_Flux = dt * (Qwest - Qeast + Qnorth - Qsouth);   % [m3]

%% ------------------------------------------------------------------------
% 10) RESERVOIR BOUNDARY CONDITION
% -------------------------------------------------------------------------
% Keep the same structural logic as in the previous D4 function.
% Here it acts on representative depth and adds/removes storage volume.
if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = drep_n(reservoir_y(ii),reservoir_x(ii)); % [m]
            dt_h  = time_step/60;                            % [h]

            available_volume = 1000 * max(dtsup - h1(ii), 0) / dt_h; % [mm/h]
            dh = min( ...
                k1(ii) * (max(dtsup - h1(ii),0))^k2(ii) / cell_area * 1000 * 3600, ...
                available_volume) * dt_h;                                  % [mm]

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = ...
                Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000 * cell_area;

            drep_n(yds1(ii),xds1(ii)) = drep_n(yds1(ii),xds1(ii)) + dh/1000;
        else
            dh = 0;
        end

        if ~isnan(yds2(ii))
            available_volume = 1000 * max(dtsup - h2(ii), 0) / dt_h; % [mm/h]
            dh = min( ...
                k3(ii) * (max(dtsup - h2(ii),0))^k4(ii) / cell_area * 1000 * 3600, ...
                available_volume) * dt_h;                                  % [mm]

            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = ...
                Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000 * cell_area;

            drep_n(yds2(ii),xds2(ii)) = drep_n(yds2(ii),xds2(ii)) + dh/1000;
        end
    end
end

%% ------------------------------------------------------------------------
% 11) UPDATE CELL VOLUME
% -------------------------------------------------------------------------
V_t = V_n + Vol_Flux;
V_t = max(V_t, 0);

%% ------------------------------------------------------------------------
% 12) INVERT VOLUME -> UPDATED REPRESENTATIVE DEPTH
% -------------------------------------------------------------------------
drep_t = hp2d_inverse_volume_uniform( ...
    SubgridTables.volume_cell, V_t, ...
    SubgridTables.dz, SubgridTables.maxDepth, Resolution);

%% ------------------------------------------------------------------------
% 13) UPDATED WATER SURFACE ELEVATION
% -------------------------------------------------------------------------
eta_t = invert_el + drep_t;

%% ------------------------------------------------------------------------
% 14) UPDATED WETTED AREA
% -------------------------------------------------------------------------
C_a = hp2d_lookup_uniform( ...
    SubgridTables.area_cell, drep_t, ...
    SubgridTables.dz, SubgridTables.maxDepth);

%% ------------------------------------------------------------------------
% 15) OUTPUT FLUXES IN SAME FORMAT AS THE PREVIOUS D4
% -------------------------------------------------------------------------
qout_left  = -[zeros(ny,1,'like',matrix_store(:,:,1)), matrix_store(:,1:end-1,1)];
qout_right =  matrix_store(:,:,1);
qout_up    =  matrix_store(:,:,2);
qout_down  = -[matrix_store(2:end,:,2); zeros(1,nx,'like',matrix_store(:,:,2))];

%% ------------------------------------------------------------------------
% 16) FINAL REPRESENTATIVE DEPTH IN MODEL UNITS [mm]
% -------------------------------------------------------------------------
d_t = drep_t * 1000;

%% ------------------------------------------------------------------------
% 17) OUTLET FLOW
% -------------------------------------------------------------------------
% Keep the same structure as in the previous D4 function, but use the
% representative depth and current wetted area from the subgrid formulation.
outlet_flow = zeros(size(d_t), 'like', d_t);
outflow(:,:,3) = 0;

if ~isempty(row_outlet)
    outlet_sub = sub2ind(size(d_t), row_outlet(:), col_outlet(:));

    % Representative outlet depth [m]
    h_out = max(d_t(outlet_sub), 0) / 1000;

    % Outlet width derived from current wetted area
    if isscalar(C_a)
        width_out = (C_a / Resolution) * ones(size(h_out), 'like', h_out);
        area_out  = C_a * ones(size(h_out), 'like', h_out);
    else
        width_out = C_a(outlet_sub) ./ Resolution;
        area_out  = C_a(outlet_sub);
    end

    % Outlet slope term
    if outlet_type == 1
        if isscalar(slope_outlet)
            sqrtS_out = sqrt(abs(slope_outlet)) * ones(size(h_out), 'like', h_out);
        else
            sqrtS_out = sqrt(abs(slope_outlet(outlet_sub)));
        end
    else
        h_safe = max(h_out, 1e-12);
        sqrtS_out = sqrt(g .* roughness_squared(outlet_sub) .* h_safe.^(-1/3));
    end

    % Save outlet hydraulic depth in Hf(:,:,3)
    Hf3 = Hf(:,:,3);
    Hf3(outlet_sub) = h_out;
    Hf(:,:,3) = Hf3;

    % Outlet discharge in model-compatible units [mm/h]
    q_out = (1 ./ roughness_cell(outlet_sub)) .* ...
            width_out .* ...
            h_out.^(5/3) .* ...
            sqrtS_out ./ ...
            area_out * 1000 * 3600;

    % Limit by available representative depth
    qmax_out = max((d_t(outlet_sub) - 1000*h_min), 0) / (time_step/60);
    q_out = min(q_out, qmax_out);

    % Store outlet flow
    outlet_flow(outlet_sub) = q_out;
    outflow(:,:,3) = outlet_flow;

    % Remove outlet discharge from updated storage
    dV_out = q_out * (time_step/60) / 1000 .* area_out;   % [m3]
    V_t(outlet_sub) = max(V_t(outlet_sub) - dV_out, 0);

    % Recompute drep, eta, and wetted area after outlet removal
    drep_t = hp2d_inverse_volume_uniform( ...
        SubgridTables.volume_cell, V_t, ...
        SubgridTables.dz, SubgridTables.maxDepth, Resolution);

    eta_t = invert_el + drep_t;

    C_a = hp2d_lookup_uniform( ...
        SubgridTables.area_cell, drep_t, ...
        SubgridTables.dz, SubgridTables.maxDepth);

    d_t = drep_t * 1000;
end

%% ------------------------------------------------------------------------
% 18) MASS BALANCE DIAGNOSTIC
% -------------------------------------------------------------------------
final_volume = nansum(nansum(V_t)); %#ok<NASGU>
% error_vol = final_volume - current_volume; %#ok<NASGU>

%% ------------------------------------------------------------------------
% 19) TOTAL FLOW THAT LEAVES THE CELL
% -------------------------------------------------------------------------
mask = outflow;
I_tot_end_cell = abs(sum(mask,3)) * dt / 1000 / 3600 * Resolution^2;  % [m3]

%% ------------------------------------------------------------------------
% 20) PLACEHOLDER OUTPUTS KEPT FOR COMPATIBILITY
% -------------------------------------------------------------------------
Qc  = 0;
Qf  = 0;
Qci = 0;
Qfi = 0;

end

function Vq = hp2d_lookup_uniform(Vtab, q, dz, maxDepth)
%--------------------------------------------------------------------------
% PURPOSE
%   Fast and clear linear interpolation on a uniform depth axis.
%
% INPUTS
%   Vtab      [ny x nx x nz]  lookup table versus representative depth
%   q         [ny x nx]       query depth above invert [m]
%   dz        scalar          table spacing [m]
%   maxDepth  scalar          maximum tabulated depth [m]
%
% OUTPUT
%   Vq        [ny x nx]       interpolated value
%--------------------------------------------------------------------------

    [ny, nx, nz] = size(Vtab);

    % 1) Clamp query depth to table range
    q = max(q, 0);
    q = min(q, maxDepth);

    % 2) Compute lower bracketing index
    %    If q = 0       -> idxL = 1
    %    If q = dz      -> idxL = 2
    %    If q = maxDepth -> handled below
    idxL = floor(q ./ dz) + 1;

    % Keep idxL valid so idxU = idxL + 1 stays inside table
    idxL = max(idxL, 1);
    idxL = min(idxL, nz - 1);

    idxU = idxL + 1;

    % 3) Compute interpolation weight between lower and upper levels
    qL = (idxL - 1) .* dz;
    w  = (q - qL) ./ dz;

    % If q is exactly at or above maxDepth, force top segment with w = 1
    atTop = (q >= maxDepth);
    idxL(atTop) = nz - 1;
    idxU(atTop) = nz;
    w(atTop)    = 1;

    % 4) Build spatial indices
    [I, J] = ndgrid(1:ny, 1:nx);

    % 5) Convert (i,j,k) -> linear indices for lower and upper levels
    indL = sub2ind([ny, nx, nz], I, J, idxL);
    indU = sub2ind([ny, nx, nz], I, J, idxU);

    % 6) Gather lower and upper values
    vL = Vtab(indL);
    vU = Vtab(indU);

    % 7) Linear interpolation
    Vq = vL + w .* (vU - vL);
end

%==========================================================================
% INVERSE VOLUME LOOKUP ON A UNIFORM DEPTH AXIS
%==========================================================================
function drep = hp2d_inverse_volume_uniform(Vtab, Vq, dz, maxDepth, Resolution)
%--------------------------------------------------------------------------
% PURPOSE
%   Invert the storage curve:
%
%       Vtab(drep) -> drep
%
%   using simple linear interpolation between adjacent tabulated depths.
%
% INPUTS
%   Vtab        [ny x nx x nz]  cell volume table [m^3]
%   Vq          [ny x nx]       queried volume [m^3]
%   dz          scalar          tabulated depth spacing [m]
%   maxDepth    scalar          maximum tabulated depth [m]
%   Resolution  scalar          coarse cell size [m]
%
% OUTPUT
%   drep        [ny x nx]       representative depth above invert [m]
%
% NOTES
%   1) Assumes Vtab is monotonic increasing in the 3rd dimension.
%   2) If Vq <= 0, returns drep = 0.
%   3) If Vq exceeds the top tabulated volume, continues above the table
%      using a prism of plan area Resolution^2:
%
%          drep = maxDepth + (Vq - Vtop) / Resolution^2
%--------------------------------------------------------------------------

    [ny, nx, nz] = size(Vtab);

    %----------------------------------------------------------------------
    % Depth levels associated with the table
    %----------------------------------------------------------------------
    depth_axis = (0:nz-1) * dz;

    %----------------------------------------------------------------------
    % Clamp negative queried volumes to zero
    %----------------------------------------------------------------------
    Vq = max(Vq, 0);

    %----------------------------------------------------------------------
    % Top volume in each cell
    %----------------------------------------------------------------------
    Vtop = Vtab(:,:,end);

    % For interpolation inside the table, clamp to the top tabulated volume.
    % Above-top cases are handled afterward with prism extrapolation.
    Vq_clamped = min(Vq, Vtop);

    %----------------------------------------------------------------------
    % Find lower bracketing index idxL such that:
    %
    %   Vtab(:,:,idxL) <= Vq_clamped <= Vtab(:,:,idxL+1)
    %
    % Since MATLAB does not directly do this per-pixel on a 3D array,
    % we use a logical count along the 3rd dimension.
    %----------------------------------------------------------------------
    idxL = sum(Vtab <= Vq_clamped, 3);

    % Keep idxL in valid range so that idxU = idxL + 1 stays inside the table
    idxL = max(idxL, 1);
    idxL = min(idxL, nz - 1);

    idxU = idxL + 1;

    %----------------------------------------------------------------------
    % Build spatial indices
    %----------------------------------------------------------------------
    [I, J] = ndgrid(1:ny, 1:nx);

    %----------------------------------------------------------------------
    % Convert (i,j,k) into linear indices for lower and upper bracketing
    % table values
    %----------------------------------------------------------------------
    indL = sub2ind([ny, nx, nz], I, J, idxL);
    indU = sub2ind([ny, nx, nz], I, J, idxU);

    %----------------------------------------------------------------------
    % Gather lower and upper tabulated volumes
    %----------------------------------------------------------------------
    V1 = Vtab(indL);
    V2 = Vtab(indU);

    %----------------------------------------------------------------------
    % Corresponding lower and upper depths
    %----------------------------------------------------------------------
    d1 = depth_axis(idxL);
    d2 = depth_axis(idxU);

    %----------------------------------------------------------------------
    % Linear interpolation in volume space
    %
    %   a = (Vq - V1) / (V2 - V1)
    %   d = d1 + a * (d2 - d1)
    %----------------------------------------------------------------------
    denom = V2 - V1;
    bad = abs(denom) <= eps(classUnderlyingLike(Vtab));

    % Prevent divide-by-zero in degenerate cases
    denom_safe = denom;
    denom_safe(bad) = 1;

    a = (Vq_clamped - V1) ./ denom_safe;
    a = max(min(a, 1), 0);
    a(bad) = 0;

    drep = d1 + a .* (d2 - d1);

    %----------------------------------------------------------------------
    % Above-top extrapolation using a prism of plan area Resolution^2
    %----------------------------------------------------------------------
    above_top = Vq > Vtop;
    if any(above_top(:))
        drep(above_top) = maxDepth + ...
            (Vq(above_top) - Vtop(above_top)) ./ (Resolution^2);
    end

    %----------------------------------------------------------------------
    % Final cleanup
    %----------------------------------------------------------------------
    drep(~isfinite(Vtab(:,:,1))) = 0;
    drep = max(drep, 0);
end


%==========================================================================
% HELPER: eps with correct underlying class for gpuArray / numeric arrays
%==========================================================================
function e = classUnderlyingLike(x)
% Returns a numeric class name suitable for eps(...)
% Works for regular arrays and gpuArray.

    if isa(x, 'gpuArray')
        e = classUnderlying(x);
    else
        e = class(x);
    end
end