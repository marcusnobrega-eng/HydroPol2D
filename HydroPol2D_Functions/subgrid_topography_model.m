function [q_face,Hf_x,Hf_y,Wf_x,Wf_y] = subgrid_topography_model( ...
    flag_numerical_scheme,eta_n,q_prev,nc,Resolution,dt,h_min,idx_nan,Subgrid_Properties,SubgridTables)
%==========================================================================
% SUBGRID_TOPOGRAPHY_MODEL
%==========================================================================
% PURPOSE
%   Compute subgrid-aware interface hydraulic properties and update the
%   local-inertial face discharges using the shared Inertial_Solver.
%
% INPUTS
%   flag_numerical_scheme : 1 Bates, 2 s-upwind, 3 s-centered
%   eta_n                 : [ny x nx] absolute water surface elevation [m]
%   q_prev                : [ny x nx x 2] previous discharge per unit width [m^2/s]
%                           q_prev(:,:,1) = x-face
%                           q_prev(:,:,2) = y-face
%   nc                    : scalar or [ny x nx] Manning n
%   Resolution            : coarse cell size / interface length [m]
%   dt                    : time step [s]
%   h_min                 : minimum operative flow depth [m]
%   idx_nan               : [ny x nx] inactive-domain mask
%   Subgrid_Properties    : currently kept for compatibility
%   SubgridTables         : lookup tables built in preprocessing
%
% OUTPUTS
%   q_face                : [ny x nx x 2] updated discharge per unit width [m^2/s]
%   Hf_x                  : [ny x nx] effective x-face flow depth [m]
%   Hf_y                  : [ny x nx] effective y-face flow depth [m]
%   Wf_x                  : [ny x nx] effective x-face width [m]
%   Wf_y                  : [ny x nx] effective y-face width [m]
%
% PAPER CONSISTENCY
%   The effective interface depth is computed as:
%
%       Hf = wetted_face_area / Resolution
%
%   which follows the paper's definition:
%
%       h_flow = wetted face area at eta_max / DeltaX
%
%   The interface width is therefore taken as:
%
%       Wf = Resolution
%
% NOTES
%   - This function only builds Hf, S, n^2, q_prev and calls Inertial_Solver.
%   - Continuity / volume update remains outside this function.
%   - We keep the 3D packed format so that all numerical schemes inside
%     Inertial_Solver continue to work as designed.
%==========================================================================

    [ny,nx] = size(eta_n);
    small = cast(1e-12,'like',eta_n);

    %----------------------------------------------------------------------
    % Inactive mask
    %----------------------------------------------------------------------
    if isempty(idx_nan)
        inactive = false(ny,nx);
    else
        inactive = logical(idx_nan);
    end

    %----------------------------------------------------------------------
    % Manning n field
    %----------------------------------------------------------------------
    if isscalar(nc)
        n_field = nc .* ones(ny,nx,'like',eta_n);
    else
        n_field = nc;
    end

    %----------------------------------------------------------------------
    % Preallocate packed arrays expected by Inertial_Solver
    %----------------------------------------------------------------------
    Hf   = zeros(ny,nx,2,'like',eta_n);
    S    = zeros(ny,nx,2,'like',eta_n);
    n_sq = zeros(ny,nx,2,'like',eta_n);
    q_p  = zeros(ny,nx,2,'like',eta_n);

    Wf_x = zeros(ny,nx,'like',eta_n);
    Wf_y = zeros(ny,nx,'like',eta_n);

    %======================================================================
    % X-FACES
    %======================================================================
    if nx > 1
        etaL = eta_n(:,1:nx-1);
        etaR = eta_n(:,2:nx);

        % Paper uses eta_max across the interface
        eta_face_x = max(etaL, etaR);

        % Shared-face wetted area from lookup table
        Ax = interp_subgrid_uniform_clear( ...
            eta_face_x, ...
            SubgridTables.area_x, ...
            SubgridTables.invert_x, ...
            SubgridTables.dz, ...
            SubgridTables.maxDepth);

        % Paper-style effective flow depth
        Hx = Ax ./ Resolution;

        % Effective interface width = full coarse interface length
        Wx = Resolution .* ones(size(Ax),'like',Ax);

        % Free-surface slope across interface
        Sx = (etaR - etaL) ./ Resolution;

        % Face Manning n = average of the two neighboring coarse cells
        nx_face = 0.5 .* (n_field(:,1:nx-1) + n_field(:,2:nx));

        % Inactive / dry x-faces
        inactive_x = inactive(:,1:nx-1) | inactive(:,2:nx);
        dry_x = (Hx <= h_min) | (Ax <= small) | inactive_x;

        % Zero-out dry/inactive locations
        Hx(dry_x)      = 0;
        Sx(dry_x)      = 0;
        Wx(dry_x)      = 0;
        nx_face(dry_x) = 0;

        % Pack arrays for Inertial_Solver
        Hf(:,1:nx-1,1)   = Hx;
        S(:,1:nx-1,1)    = Sx;
        n_sq(:,1:nx-1,1) = nx_face.^2;
        q_p(:,1:nx-1,1)  = q_prev(:,1:nx-1,1);

        Wf_x(:,1:nx-1)   = Wx;
    end

    %======================================================================
    % Y-FACES
    %======================================================================
    if ny > 1
        etaU = eta_n(1:ny-1,:);
        etaD = eta_n(2:ny,:);

        % Paper uses eta_max across the interface
        eta_face_y = max(etaU, etaD);

        % Shared-face wetted area from lookup table
        Ay = interp_subgrid_uniform_clear( ...
            eta_face_y, ...
            SubgridTables.area_y, ...
            SubgridTables.invert_y, ...
            SubgridTables.dz, ...
            SubgridTables.maxDepth);

        % Paper-style effective flow depth
        Hy = Ay ./ Resolution;

        % Effective interface width = full coarse interface length
        Wy = Resolution .* ones(size(Ay),'like',Ay);

        % Free-surface slope across interface
        Sy = (etaD - etaU) ./ Resolution;

        % Face Manning n = average of the two neighboring coarse cells
        ny_face = 0.5 .* (n_field(1:ny-1,:) + n_field(2:ny,:));

        % Inactive / dry y-faces
        inactive_y = inactive(1:ny-1,:) | inactive(2:ny,:);
        dry_y = (Hy <= h_min) | (Ay <= small) | inactive_y;

        % Zero-out dry/inactive locations
        Hy(dry_y)      = 0;
        Sy(dry_y)      = 0;
        Wy(dry_y)      = 0;
        ny_face(dry_y) = 0;

        % Pack arrays for Inertial_Solver
        Hf(1:ny-1,:,2)   = Hy;
        S(1:ny-1,:,2)    = Sy;
        n_sq(1:ny-1,:,2) = ny_face.^2;
        q_p(1:ny-1,:,2)  = q_prev(1:ny-1,:,2);

        Wf_y(1:ny-1,:)   = Wy;
    end

    %----------------------------------------------------------------------
    % Call the shared inertial solver
    %----------------------------------------------------------------------
    q_face = Inertial_Solver( ...
        flag_numerical_scheme, q_p, dt, Hf, S, n_sq, Resolution, idx_nan);

    %----------------------------------------------------------------------
    % Final cleanup
    %----------------------------------------------------------------------
    q_face(~isfinite(q_face)) = 0;
    Hf(~isfinite(Hf)) = 0;

    Hf_x = Hf(:,:,1);
    Hf_y = Hf(:,:,2);

    Hf_x(Hf_x <= small) = 0;
    Hf_y(Hf_y <= small) = 0;

    Wf_x(~isfinite(Wf_x)) = 0;
    Wf_y(~isfinite(Wf_y)) = 0;
end


%==========================================================================
% CLEAR, UNIFORM-AXIS INTERPOLATOR FOR SUBGRID TABLES
%==========================================================================
function out = interp_subgrid_uniform_clear(eta, value_table, invert_el, dz, maxDepth)
%--------------------------------------------------------------------------
% PURPOSE
%   Clear matrix-wise interpolation for subgrid face/cell tables defined on
%   a uniform depth axis above invert.
%
% INPUTS
%   eta         [m x n]      query water-surface elevation [m]
%   value_table [m x n x nz] tabulated quantity versus depth-above-invert
%   invert_el   [m x n]      invert elevation [m]
%   dz          scalar       uniform depth spacing [m]
%   maxDepth    scalar       maximum tabulated depth [m]
%
% OUTPUT
%   out         [m x n]      interpolated value
%
% INTERPRETATION
%   The table is assumed to live on:
%
%       depth = 0, dz, 2*dz, ..., maxDepth
%
%   so the query depth is:
%
%       d = eta - invert_el
%--------------------------------------------------------------------------

    [m,n,nz] = size(value_table);

    % Query depth above invert
    d = eta - invert_el;

    % Dry points
    dry = (d <= 0);

    % Clamp to table range
    d = max(d, 0);
    d = min(d, maxDepth);

    % Lower and upper indices
    idxL = floor(d ./ dz) + 1;
    idxL = max(idxL, 1);
    idxL = min(idxL, nz-1);
    idxU = idxL + 1;

    % Interpolation weight
    dL = (idxL - 1) .* dz;
    w  = (d - dL) ./ dz;

    % Top handling
    atTop = (d >= maxDepth);
    idxL(atTop) = nz - 1;
    idxU(atTop) = nz;
    w(atTop)    = 1;

    % Build spatial indices
    [I,J] = ndgrid(1:m, 1:n);

    % Convert to linear indices in the 3D table
    indL = sub2ind([m,n,nz], I, J, idxL);
    indU = sub2ind([m,n,nz], I, J, idxU);

    % Gather lower and upper values
    vL = value_table(indL);
    vU = value_table(indU);

    % Linear interpolation
    out = vL + w .* (vU - vL);

    % Dry points -> zero
    out(dry) = 0;
end