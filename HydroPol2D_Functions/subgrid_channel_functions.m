function [outflow, C_a, Hf_x, Hf_y, Wf_x, Wf_y] = subgrid_channel_functions( ...
    depth_cell, River_Width, zbed, Resolution, nc, Qc_prev, Qci_prev, ...
    g, dt, idx_rivers, Subgrid_Properties, SubgridTables)
%SUBGRID_CHANNEL_FUNCTIONS
%
% PURPOSE
%   Compute local-inertial fluxes using SHARED-FACE subgrid lookup tables,
%   compute cell storage area for continuity, and return face hydraulic
%   quantities to the caller so they are not reconstructed twice.
%
% INPUTS
%   depth_cell   [ny x nx]     water depth above zbed [m]
%   River_Width  [ny x nx]     kept for compatibility (not used here)
%   zbed         [ny x nx]     bed elevation used to compute WSE [m]
%   Resolution   scalar        coarse grid size [m]
%   nc           scalar/array  Manning n
%   Qc_prev      [ny x nx x 2] previous discharge on faces [m3/s]
%   Qci_prev     placeholder, unused
%   g, dt        scalars
%   idx_rivers   [ny x nx]     logical mask for active donor cells
%   Subgrid_Properties.invert_el [ny x nx] cell invert elevation [m]
%   SubgridTables shared-face uniform-dz lookup tables
%   Sx_face      [ny x (nx-1)] hydrostatic x-face slope from caller
%   Sy_face      [(ny-1) x nx] hydrostatic y-face slope from caller
%
% OUTPUTS
%   outflow      [ny x nx x 2] q = Q/Wface [m2/s]
%   C_a          [ny x nx]     storage area [m2]
%   Hf_x,Hf_y    [ny x nx]     hydraulic depth on padded faces [m]
%   Wf_x,Wf_y    [ny x nx]     wetted width on padded faces [m]
%
% NOTES
%   - x-faces are between (i,j) and (i,j+1), size [ny x (nx-1)]
%   - y-faces are between (i,j) and (i+1,j), size [(ny-1) x nx]
%   - geometry is reconstructed ONCE here and passed back to caller

% -------------------------------------------------------------------------
% 0) Sizes and parameters
% -------------------------------------------------------------------------
ny = size(depth_cell,1);
nx = size(depth_cell,2);

n_chan = nc;

% -------------------------------------------------------------------------
% 1) Absolute WSE and storage depth
% -------------------------------------------------------------------------
z_inv_cell = Subgrid_Properties.invert_el;   % [ny x nx]
eta = zbed + depth_cell;                     % [ny x nx] absolute WSE

% -------------------------------------------------------------------------
% 1b) Hydrostatic face slopes (computed INSIDE this function)
% -------------------------------------------------------------------------

% ===== X faces: [ny x (nx-1)] =====
zf_x = max(zbed(:,1:end-1), zbed(:,2:end));
etaL_x = max(eta(:,1:end-1), zf_x);
etaR_x = max(eta(:,2:end),   zf_x);
Sx_face = (etaR_x - etaL_x) ./ Resolution;

% ===== Y faces: [(ny-1) x nx] =====
zf_y = max(zbed(1:end-1,:), zbed(2:end,:));
etaS_y = max(eta(1:end-1,:), zf_y);
etaN_y = max(eta(2:end,:),   zf_y);
Sy_face = (etaS_y - etaN_y) ./ Resolution;

d_store = max(eta - z_inv_cell, 0);          % [ny x nx]
C_a = hp2d_lookup_uniform_shared(SubgridTables.area, d_store, SubgridTables.dz, SubgridTables.maxDepth);
C_a(~isfinite(C_a)) = 0;
C_a = max(C_a, 0);

% -------------------------------------------------------------------------
% 2) Shared-face depths and hydraulic properties (computed ONCE)
% -------------------------------------------------------------------------
% =========================
% X faces: [ny x (nx-1)]
% =========================
etaMax_x = max(eta(:,1:end-1), eta(:,2:end));
d_x      = max(etaMax_x - SubgridTables.invert_x, 0);

Aface_x  = hp2d_lookup_uniform_shared(SubgridTables.area_x,  d_x, SubgridTables.dz, SubgridTables.maxDepth);
Wface_x  = hp2d_lookup_uniform_shared(SubgridTables.width_x, d_x, SubgridTables.dz, SubgridTables.maxDepth);
Rhface_x = hp2d_lookup_uniform_shared(SubgridTables.Rh_x,    d_x, SubgridTables.dz, SubgridTables.maxDepth);

% =========================
% Y faces: [(ny-1) x nx]
% =========================
etaMax_y = max(eta(1:end-1,:), eta(2:end,:));
d_y      = max(etaMax_y - SubgridTables.invert_y, 0);

Aface_y  = hp2d_lookup_uniform_shared(SubgridTables.area_y,  d_y, SubgridTables.dz, SubgridTables.maxDepth);
Wface_y  = hp2d_lookup_uniform_shared(SubgridTables.width_y, d_y, SubgridTables.dz, SubgridTables.maxDepth);
Rhface_y = hp2d_lookup_uniform_shared(SubgridTables.Rh_y,    d_y, SubgridTables.dz, SubgridTables.maxDepth);

% Hygiene
Aface_x(~isfinite(Aface_x))   = 0; Aface_x  = max(Aface_x, 0);
Wface_x(~isfinite(Wface_x))   = 0; Wface_x  = max(Wface_x, 0);
Rhface_x(~isfinite(Rhface_x)) = 0; Rhface_x = max(Rhface_x, 0);

Aface_y(~isfinite(Aface_y))   = 0; Aface_y  = max(Aface_y, 0);
Wface_y(~isfinite(Wface_y))   = 0; Wface_y  = max(Wface_y, 0);
Rhface_y(~isfinite(Rhface_y)) = 0; Rhface_y = max(Rhface_y, 0);

% Padded outputs back to caller
Hf_x = zeros(ny,nx,'like',depth_cell);
Hf_y = zeros(ny,nx,'like',depth_cell);
Wf_x = zeros(ny,nx,'like',depth_cell);
Wf_y = zeros(ny,nx,'like',depth_cell);

Hf_x(:,1:end-1) = max(Aface_x ./ Resolution, 0);
Hf_y(2:end,:)   = max(Aface_y ./ Resolution, 0);

Wf_x(:,1:end-1) = Wface_x;
Wf_y(2:end,:)   = Wface_y;

% -------------------------------------------------------------------------
% 3) Previous discharges and donor-cell roughness
% -------------------------------------------------------------------------
Qx_old = Qc_prev(:,1:end-1,1);   % [ny x (nx-1)] m3/s
Qy_old = Qc_prev(1:end-1,:,2);   % [(ny-1) x nx] m3/s

n_x = n_chan(:,1:end-1);
n_y = n_chan(1:end-1,:);

% -------------------------------------------------------------------------
% 4) Active masks
% -------------------------------------------------------------------------
if isempty(idx_rivers)
    active_x = true(size(Aface_x));
    active_y = true(size(Aface_y));
else
    active_x = logical(idx_rivers(:,1:end-1));
    active_y = logical(idx_rivers(1:end-1,:));
end

active_x = active_x & (Aface_x > 0) & (Wface_x > 0) & (d_x > 0);
active_y = active_y & (Aface_y > 0) & (Wface_y > 0) & (d_y > 0);

% -------------------------------------------------------------------------
% 5) Local-inertial update using caller-provided hydrostatic slopes
% -------------------------------------------------------------------------
AxRh_x = max(Aface_x .* max(Rhface_x,1e-12).^(4/3), 1e-12);
AxRh_y = max(Aface_y .* max(Rhface_y,1e-12).^(4/3), 1e-12);

den_x = 1 + g*dt .* (n_x.^2) .* abs(Qx_old) ./ AxRh_x;
den_y = 1 + g*dt .* (n_y.^2) .* abs(Qy_old) ./ AxRh_y;

Qx_new = zeros(size(Qx_old), 'like', Qx_old);
Qy_new = zeros(size(Qy_old), 'like', Qy_old);

Qx_new(active_x) = (Qx_old(active_x) - g*dt .* Aface_x(active_x) .* Sx_face(active_x)) ./ den_x(active_x);
Qy_new(active_y) = (Qy_old(active_y) - g*dt .* Aface_y(active_y) .* Sy_face(active_y)) ./ den_y(active_y);

% -------------------------------------------------------------------------
% 6) Convert to q = Q / Wface and pack
% -------------------------------------------------------------------------
qx = zeros(size(Qx_new), 'like', Qx_new);
qy = zeros(size(Qy_new), 'like', Qy_new);

qx(active_x) = Qx_new(active_x) ./ max(Wface_x(active_x), 1e-12);
qy(active_y) = Qy_new(active_y) ./ max(Wface_y(active_y), 1e-12);

outflow = zeros(ny,nx,2,'like',depth_cell);
outflow(:,1:end-1,1) = qx;
outflow(1:end-1,:,2) = qy;

outflow(~isfinite(outflow)) = 0;

end