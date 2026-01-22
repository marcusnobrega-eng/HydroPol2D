function [h_t, u_x, u_y, q_exf, q_river, error] = Boussinesq_2D_explicit(dt,dx,dy,h,z0,Sy,R,K,river_mask,K_river,h_river,z_river,Courant,h_soil,catchment_mask,dirichlet_mask,h_dirichlet,perimeter)
%% ═══════════════════════════════════════════════════════════════════════
%  Function: Boussinesq_2D_explicit
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 03/06/2025
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Solve the 2D Boussinesq equation using an explicit finite difference method
%      with adaptive time-stepping based on the Courant stability condition.
%      This function computes the updated groundwater head, velocities, exfiltration
%      rates, and river interactions. It enforces Dirichlet and no-flow boundary
%      conditions while performing mass balance checks.
%
%  ➤ Inputs:
%      • dt            - Initial model time-step [sec]
%      • dx            - Grid spacing in x-direction [m]
%      • dy            - Grid spacing in y-direction [m]
%      • h             - Initial hydraulic head [m] (Nx-by-Ny matrix)
%      • z0            - Ground elevation [m] (Nx-by-Ny matrix)
%      • Sy            - Specific yield coefficient [1/m]
%      • R             - Recharge rate [m/s] (Nx-by-Ny matrix)
%      • K             - Hydraulic conductivity [m/s] (Nx-by-Ny matrix)
%      • river_mask    - Logical mask indicating river locations (Nx-by-Ny)
%      • K_river       - Hydraulic conductivity of the river [m/s]
%      • h_river       - Hydraulic head at the rivers (matrix matching DEM size)
%      • z_river       - Topographic river elevation [m or mm]
%      • Courant       - Stability factor (0 < Courant ≤ 1)
%      • h_soil        - Soil depth [m] (Nx-by-Ny matrix)
%      • catchment_mask- Binary mask for the computational domain (Nx-by-Ny)
%      • dirichlet_mask- Binary mask for Dirichlet boundary conditions (Nx-by-Ny)
%      • h_dirichlet   - Prescribed hydraulic head at Dirichlet boundaries [m] (Nx-by-Ny)
%      • perimeter     - Binary mask for no-flow boundaries (Nx-by-Ny)
%
%  ➤ Outputs:
%      • h_t           - Final hydraulic head [m] (Nx-by-Ny matrix)
%      • u_x           - Groundwater velocity in x-direction [m/s] (Nx-1-by-Ny)
%      • u_y           - Groundwater velocity in y-direction [m/s] (Nx-by-Ny-1)
%      • q_exf         - Exfiltration rate at cells [m/s] (Nx-by-Ny)
%      • q_river       - River exfiltration rate [m/s] (positive: leaves, negative: enters aquifer)
%      • error         - Mass balance error (e.g., in m³)
%
%  ➤ Local Functions:
%      • compute_gradients         - Computes second derivatives of h using finite differences.
%      • compute_stable_dt         - Determines a stable time-step based on groundwater velocities.
%      • apply_free_flow_bc        - Applies free-flow boundary conditions.
%      • applyNoFlowBC             - Enforces no-flow (Neumann) boundary conditions along the perimeter.
%      • compute_slopes_from_dem   - Computes slopes in x and y from a DEM.
%      • compute_boussinesq_velocities - Computes groundwater velocities using Darcy's law.
%
%  ➤ Notes:
%      • Adaptive time-stepping is used to ensure numerical stability (Courant condition).
%      • Both Dirichlet and no-flow boundary conditions are enforced.
%      • Mass balance is monitored to flag any potential numerical issues.
%      • Some river flux calculations are commented out, indicating options for extended
%        functionality.
%      • Local functions support gradient, velocity, and boundary condition calculations.
% ═══════════════════════════════════════════════════════════════════════

perimeter = logical(perimeter);
% h_dirichlet = z0; % DELETEEE
% K(:,end) = 1000*K(:,end);
% 
% dirichlet_mask = perimeter; % DELETE
% dirichlet_mask(:,end-1) = 0;
if isempty(dirichlet_mask) || isempty(h_dirichlet) 
    dirichlet_mask = false(size(z0,1),size(z0,2));
    h_dirichlet = false(size(z0,1),size(z0,2));
end

%% Surface Water Level
h_surf = z0 + h_soil;  % Compute the surface elevation from topography

%% INITIAL CONDITION
q_exf = zeros(size(h)); % Initialize exfiltration rate matrix

%% River Flux Length 
L_river = 1/2*(dx+dy);

%% TIME STEPPING LOOP
t = 0; % Initialize time
T = dt; % Set total simulation time to initial dt

%% Loop
h_0 = h;
h_t = h_0; % Copy current head for updating
Nx = size(h,2); Ny = size(h,1);

% Velocity Field
% Compute groundwater velocities
[u_x, u_y] = compute_boussinesq_velocities(h_t, K, dx, dy); % m/s
 

while t < T
    n_steps = ceil(T/dt);
    kk = 1;

    % Compute new stable dt at each step
    dt_GW = nanmin(nanmin(compute_stable_dt(u_x, u_y, dx, dy, Courant))); % Adaptive time step
    
    % Ensure dt does not exceed remaining simulation time
    dt = min(dt_GW,dt);

    % Compute groundwater discharge to rivers (maybe this has to be
    % deactivated)
    q_river = zeros(Ny, Nx);  % Initialize discharge array

    % Compute second derivatives using central differencing
    depth_threshold = 0.00; % m

    % [d2h_dx2, d2h_dy2] = compute_gradients(h_t, z0, K, dx, dy, depth_threshold, perimeter);
    % === Predictor step ===
    [Fx1, Fy1] = compute_fluxes_conservative(h_t, z0, K, dx, dy, dt, Sy);
    div1 = compute_flux_divergence(Fx1, Fy1, dx, dy);
    h_predict = h_t + dt ./ Sy .* (-div1);

    % Apply Dirichlet and physical constraints to predictor
    h_predict(dirichlet_mask) = h_dirichlet(dirichlet_mask);
    h_predict = max(h_predict, z0);  % no negative storage

    % === Corrector step ===
    [Fx2, Fy2] = compute_fluxes_conservative(h_predict, z0, K, dx, dy, dt, Sy);

    % === Average interface fluxes ===
    Fx_avg = 0.5 * (Fx1 + Fx2);
    Fy_avg = 0.5 * (Fy1 + Fy2);

    % === Final update using averaged flux divergence ===
    div_avg = compute_flux_divergence(Fx_avg, Fy_avg, dx, dy);
    h_t = h_t + dt ./ Sy .* (-div_avg + R);


    % Apply Dirichlet boundary conditions
    h_t(dirichlet_mask) = h_dirichlet(dirichlet_mask); 

    % Lost Mass
    lost_mass = nansum(nansum(Sy .* dx * dy * (-1) .* min(h_t-z0,0)));

    % Ensure groundwater head does not fall below ground level
    h_t = max(h_t, z0); % m

    % Compute exfiltration where groundwater head exceeds surface elevation 
    q_exf = max(0, (h_t - h_surf)) .* Sy / dt;  % [m/s]

    exf_vol = dt*nansum(nansum(q_exf*dx*dy)); % exiltration volume [m3]

    % Groundwater Head Considering Exfiltration
    h_t = min(h_t, h_surf); % m

    % Compute groundwater velocities
    [u_x, u_y] = compute_boussinesq_velocities(h_t, K, dx, dy); % m/s

    % Apply catchment mask (ensuring calculations remain within domain)
    h_t(~catchment_mask) = nan;
    u_x(~catchment_mask) = nan;
    u_y(~catchment_mask) = nan;
    q_exf(~catchment_mask) = nan;

    % Visualization (Remove for efficiency in large simulations)
    % surf(h_t);
    % title(sprintf('2D Groundwater Flow (Time = %.2f s)', t));
    % view(0,90);
    % colormap('turbo');
    % colorbar;
    % pause(0.00001);

    % Mass balance check
    cell_area = dx*dy;
    % === Correct per-timestep mass balance tracking ===
    
    % Step 1: Recharge input this step [m³]
    recharge_mass = nansum(nansum(R .* dt)) * cell_area;
    
    % Step 2: Exfiltration loss this step [m³]
    exfil_mass = nansum(nansum(q_exf .* dt)) * cell_area;
    
    % Step 3: Storage change this step [m³]
    storage_prev = nansum(nansum(Sy .* max(h_0 - z0, 0))) * cell_area;   % before update
    storage_curr = nansum(nansum(Sy .* max(h_t - z0, 0))) * cell_area; % after update
    storage_change = storage_curr - storage_prev;
       
    % Step 4: Mass balance error [m³]
    error = recharge_mass - exfil_mass - storage_change;
   

    % Update variables for next time step
    h = h_t;
    kk = kk + 1;
    if kk == (n_steps -1) % Last time-step
        dt = T - kk*dt;
    end
    t = t + dt;
end


end

%% LOCAL FUNCTIONS

% Compute harmonic mean of two values
function h_mean = mean_function(K1, K2,flag_mean,threshold)
eps = 1e-12;
if K1 + K2 == 0
    h_mean = 0; % Avoid division by zero
else
    if flag_mean == 1
        h_mean = 2 * (K1 .* K2) ./ (K1 + K2) + eps;

    else
        h_mean = 1/2*(K1 + K2);
    end
    if ~isempty(threshold)
        diff = abs(K1 - K2);
        h_mean(diff < threshold) = 0;
    end
end
end


function [d2h_dx2, d2h_dy2] = compute_gradients(h, z0, K, dx, dy, depth_threshold, perimeter)
    % COMPUTE_GRADIENTS
    % Computes second derivatives (divergence of groundwater fluxes) for explicit solver.
    % Assumes:
    %   x = columns (horizontal)
    %   y = rows (vertical)

    % Water depth above base
    H = max(h - z0, 0);  % [m]
    H(isnan(z0)) = nan;
    [nRows, nCols] = size(h);  % y = rows, x = columns

    % Gradient limiter (inactive by default)
    alpha = 1;
    L = 0.5 * (dx + dy);
    G_max = max(alpha * (H / L), 1);

    % === X-DIRECTION (horizontal / across columns) ===
    Hx = 0.5 * (H(:,1:end-1) + H(:,2:end));         
    Kx = 0.5 * (K(:,1:end-1) + K(:,2:end));         
    dHdx = (h(:,2:end) - h(:,1:end-1)) / dx;        
    Fx = -Kx .* Hx .* dHdx;                        

    % Enforce zero outflow: check if h(:,j+1) is NaN
    right_is_nan = isnan(h(:,2:end));  % size = (nRows, nCols-1)
    Fx(right_is_nan) = 0;

    % Divergence in x
    divFx = zeros(nRows, nCols);
    divFx(:,2:nCols-1) = (Fx(:,2:end) - Fx(:,1:end-1)) / dx;
    divFx(:,1) = Fx(:,1) / dx;
    divFx(:,nCols) = -Fx(:,end) / dx;

    % === Y-DIRECTION (vertical / across rows) ===
    Hy = 0.5 * (H(1:end-1,:) + H(2:end,:));         
    Ky = 0.5 * (K(1:end-1,:) + K(2:end,:));
    dHdy = (h(2:end,:) - h(1:end-1,:)) / dy;        
    Fy = -Ky .* Hy .* dHdy;

    % Enforce zero outflow: check if h(i+1,:) is NaN
    down_is_nan = isnan(h(2:end,:));  % size = (nRows-1, nCols)
    Fy(down_is_nan) = 0;

    % Divergence in y
    divFy = zeros(nRows, nCols);
    divFy(2:nRows-1,:) = (Fy(2:end,:) - Fy(1:end-1,:)) / dy;
    divFy(1,:) = Fy(1,:) / dy;
    divFy(nRows,:) = -Fy(end,:) / dy;

    % Output second derivatives
    d2h_dx2 = divFx;
    d2h_dy2 = divFy;

    % Optional cleanup (could help visualizations)
    d2h_dx2(isnan(z0)) = nan;
    d2h_dy2(isnan(z0)) = nan;
end

function [Fx, Fy] = compute_fluxes_conservative(h, z0, K, dx, dy, dt, Sy)
    % Returns conservative, flux-limited interface fluxes
    H = max(h - z0, 0);
    [nRows, nCols] = size(h);

    % === X-direction fluxes ===
    Hx = 0.5 * (H(:,1:end-1) + H(:,2:end));
    Kx = 0.5 * (K(:,1:end-1) + K(:,2:end));
    dHdx = (h(:,2:end) - h(:,1:end-1)) / dx;
    Fx = -Kx .* Hx .* dHdx;

    donor_Hx = H(:,1:end-1);
    max_flux_x = Sy(:,1:end-1) .* donor_Hx / dt;
    Fx = sign(Fx) .* min(abs(Fx), max_flux_x);
    Fx(isnan(H(:,2:end))) = 0;

    % === Y-direction fluxes ===
    Hy = 0.5 * (H(1:end-1,:) + H(2:end,:));
    Ky = 0.5 * (K(1:end-1,:) + K(2:end,:));
    dHdy = (h(2:end,:) - h(1:end-1,:)) / dy;
    Fy = -Ky .* Hy .* dHdy;

    donor_Hy = H(1:end-1,:);
    max_flux_y = Sy(1:end-1,:) .* donor_Hy / dt;
    Fy = sign(Fy) .* min(abs(Fy), max_flux_y);
    Fy(isnan(H(2:end,:))) = 0;
end

function div = compute_flux_divergence(Fx, Fy, dx, dy)
    % Converts interface fluxes into divergence (per cell)
    [nRows, nCols_minus1] = size(Fx);
    nCols = nCols_minus1 + 1;
    div_x = zeros(nRows, nCols);
    div_x(:,2:end-1) = (Fx(:,2:end) - Fx(:,1:end-1)) / dx;
    div_x(:,1) = Fx(:,1) / dx;
    div_x(:,end) = -Fx(:,end) / dx;

    [nRows_minus1, nCols] = size(Fy);
    nRows = nRows_minus1 + 1;
    div_y = zeros(nRows, nCols);
    div_y(2:end-1,:) = (Fy(2:end,:) - Fy(1:end-1,:)) / dy;
    div_y(1,:) = Fy(1,:) / dy;
    div_y(end,:) = -Fy(end,:) / dy;

    div = div_x + div_y;
end


% function [d2h_dx2, d2h_dy2, mass_error] = compute_gradients_conservative(h, z0, K, dx, dy, dt, Sy)
%     % Mass-conservative groundwater gradient computation
%     % Ensures symmetric fluxes and performs internal mass balance check
% 
%     [nRows, nCols] = size(h);
%     H = max(h - z0, 0);  
%     H(isnan(z0)) = nan;
% 
%     % Initialize interface fluxes
%     Fx = zeros(nRows, nCols-1);  % x-direction (cols)
%     Fy = zeros(nRows-1, nCols);  % y-direction (rows)
% 
%     %% === X-direction fluxes (left to right) ===
%     for j = 1:nCols-1
%         hL = h(:, j); hR = h(:, j+1);
%         HL = H(:, j); HR = H(:, j+1);
%         Kx = 0.5 * (K(:, j) + K(:, j+1));
%         dHdx = (hR - hL) / dx;
% 
%         raw_flux = -Kx .* 0.5 .* (HL + HR) .* dHdx;
%         max_flux = Sy(:, j) .* HL / dt;
% 
%         capped_flux = sign(raw_flux) .* min(abs(raw_flux), max_flux);
%         capped_flux(isnan(HR)) = 0;  % No flux into NaN
% 
%         Fx(:, j) = capped_flux;
%     end
% 
%     %% === Y-direction fluxes (top to bottom) ===
%     for i = 1:nRows-1
%         hT = h(i, :); hB = h(i+1, :);
%         HT = H(i, :); HB = H(i+1, :);
%         Ky = 0.5 * (K(i, :) + K(i+1, :));
%         dHdy = (hB - hT) / dy;
% 
%         raw_flux = -Ky .* 0.5 .* (HT + HB) .* dHdy;
%         max_flux = Sy(i, :) .* HT / dt;
% 
%         capped_flux = sign(raw_flux) .* min(abs(raw_flux), max_flux);
%         capped_flux(isnan(HB)) = 0;
% 
%         Fy(i, :) = capped_flux;
%     end
% 
%     %% === Compute divergence from fluxes ===
%     d2h_dx2 = zeros(nRows, nCols);
%     d2h_dy2 = zeros(nRows, nCols);
% 
%     d2h_dx2(:, 2:end-1) = (Fx(:, 2:end) - Fx(:, 1:end-1)) / dx;
%     d2h_dx2(:, 1) = Fx(:, 1) / dx;
%     d2h_dx2(:, end) = -Fx(:, end) / dx;
% 
%     d2h_dy2(2:end-1, :) = (Fy(2:end, :) - Fy(1:end-1, :)) / dy;
%     d2h_dy2(1, :) = Fy(1, :) / dy;
%     d2h_dy2(end, :) = -Fy(end, :) / dy;
% 
%     %% === Mass Balance Check ===
%     total_flux_out = nansum(Fx(:)) * dy + nansum(Fy(:)) * dx;
%     total_divergence = nansum(nansum((d2h_dx2 + d2h_dy2) * dx * dy));
%     mass_error = total_flux_out - total_divergence;
% 
%     if abs(mass_error) > 1e-10
%         warning("⚠️ Mass imbalance in compute_gradients_conservative: %.3e m³", mass_error);
%     end
% 
%     % Optional cleanup
%     d2h_dx2(isnan(z0)) = nan;
%     d2h_dy2(isnan(z0)) = nan;
% end



function [d2h_dx2, d2h_dy2] = compute_gradients_conservative(h, z0, K, dx, dy, dt, Sy)
    % Computes divergence of groundwater flux with per-interface flux limiting
    % Ensures no cell loses more mass than it contains (local conservation)

    [nRows, nCols] = size(h);
    H = max(h - z0, 0);  % Depth above bedrock

    % === X-direction (across columns) ===
    Hx = 0.5 * (H(:,1:end-1) + H(:,2:end));
    Kx = 0.5 * (K(:,1:end-1) + K(:,2:end));
    dHdx = (h(:,2:end) - h(:,1:end-1)) / dx;
    Fx = -Kx .* Hx .* dHdx;  % Darcy flux from col j to j+1 (same row)

    % Limit flux by available mass in donor cells (column j)
    donor_Hx = H(:,1:end-1);  % depth in donor cell
    max_flux = Sy(:,1:end-1) .* donor_Hx / dt;  % [m/s] maximum outflow allowed
    Fx = sign(Fx) .* min(abs(Fx), max_flux);  % cap the flux but preserve direction

    % Enforce no flux into NaNs
    Fx(isnan(H(:,2:end))) = 0;

    % Divergence in x-direction
    d2h_dx2 = zeros(nRows, nCols);
    d2h_dx2(:,2:end-1) = (Fx(:,2:end) - Fx(:,1:end-1)) / dx;
    d2h_dx2(:,1) = Fx(:,1) / dx;
    d2h_dx2(:,end) = -Fx(:,end) / dx;

    % === Y-direction (across rows) ===
    Hy = 0.5 * (H(1:end-1,:) + H(2:end,:));
    Ky = 0.5 * (K(1:end-1,:) + K(2:end,:));
    dHdy = (h(2:end,:) - h(1:end-1,:)) / dy;
    Fy = -Ky .* Hy .* dHdy;  % Darcy flux from row i to i+1 (same column)

    % Limit flux by available mass in donor cells (row i)
    donor_Hy = H(1:end-1,:);
    max_flux_y = Sy(1:end-1,:) .* donor_Hy / dt;
    Fy = sign(Fy) .* min(abs(Fy), max_flux_y);

    Fy(isnan(H(2:end,:))) = 0;

    % Divergence in y-direction
    d2h_dy2 = zeros(nRows, nCols);
    d2h_dy2(2:end-1,:) = (Fy(2:end,:) - Fy(1:end-1,:)) / dy;
    d2h_dy2(1,:) = Fy(1,:) / dy;
    d2h_dy2(end,:) = -Fy(end,:) / dy;

    % Mask NaNs
    d2h_dx2(isnan(z0)) = nan;
    d2h_dy2(isnan(z0)) = nan;
end



function dt_max = compute_stable_dt(u_x, u_y, dx, dy, Courant)
    dt_max = min(Courant*dx./abs(u_x), Courant*dy./abs(u_y));
    dt_max = min(min(dt_max));
    dt_max(isinf(dt_max)) = nan;
end
% 
function h = apply_free_flow_bc(h, mask)
[Ny, Nx] = size(h);

for i = 2:Ny-1
    for j = 2:Nx-1
        if ~mask(i, j) % If it's outside the catchment
            % free_flow BC: Set to nearest interior neighbor
            if mask(i-1, j), h(i, j) = h(i-1, j); end
            if mask(i+1, j), h(i, j) = h(i+1, j); end
            if mask(i, j-1), h(i, j) = h(i, j-1); end
            if mask(i, j+1), h(i, j) = h(i, j+1); end
        end
    end
end
end


function h = applyNoFlowBC(h, perimeterMask)
    % APPLYNOFLOWBC Enforces no-flow boundary conditions using a perimeter mask.
    %
    % This function modifies the groundwater table height matrix h to satisfy
    % zero-gradient (Neumann) boundary conditions **only at the perimeter** of the domain.
    %
    % INPUTS:
    %   h            - Nx-by-Ny matrix representing the groundwater table height (m).
    %   perimeterMask - Nx-by-Ny logical matrix where 1 indicates the perimeter cells.
    %
    % OUTPUT:
    %   h - Updated groundwater table with no-flow boundary conditions applied.
    %
    % Example usage:
    %   h = applyNoFlowBC(h, perimeterMask);
    
    % Find indices of perimeter points
    [rows, cols] = find(perimeterMask);
    
    % Loop through perimeter cells and apply no-flow condition
    for k = 1:length(rows)
        i = rows(k);
        j = cols(k);
        
        % Identify the nearest inner neighbor
        if i > 1 && ~perimeterMask(i-1, j) % Check above
            h(i, j) = h(i-1, j);
        elseif i < size(h, 1) && ~perimeterMask(i+1, j) % Check below
            h(i, j) = h(i+1, j);
        elseif j > 1 && ~perimeterMask(i, j-1) % Check left
            h(i, j) = h(i, j-1);
        elseif j < size(h, 2) && ~perimeterMask(i, j+1) % Check right
            h(i, j) = h(i, j+1);
        end
    end
end


function [Sx, Sy] = compute_slopes_from_dem(DEM, dx, dy)
% COMPUTE_SLOPES_FROM_DEM - Computes the slopes in the x and y directions from a DEM.
%
% Inputs:
% - DEM: Digital Elevation Model matrix [m] (nx, ny)
% - dx: Grid spacing in the x-direction [m]
% - dy: Grid spacing in the y-direction [m]
%
% Outputs:
% - Sx: Slope in the x-direction (nx, ny)
% - Sy: Slope in the y-direction (nx, ny)

% Get size of DEM
[nx, ny] = size(DEM);

% Initialize slope matrices
Sx = zeros(nx, ny);
Sy = zeros(nx, ny);

% Compute slope in x-direction (central difference)
% Sx(2:nx-1, :) = (DEM(3:nx, :) - DEM(1:nx-2, :)) / (2 * dx);
Sx(2:nx-1, :) = (DEM(2:nx-1, :) - DEM(1:nx-2, :)) / (dx);

% Compute slope in y-direction (central difference)
% Sy(:, 2:ny-1) = (DEM(:, 3:ny) - DEM(:, 1:ny-2)) / (2 * dy);
Sy(:, 2:ny-1) = (DEM(:, 2:ny-1) - DEM(:, 1:ny-2)) / (dy);

% Boundary conditions (one-sided difference at the boundaries)
Sx(1, :) = (DEM(2, :) - DEM(1, :)) / dx;  % Forward difference for left boundary
Sx(nx, :) = (DEM(nx, :) - DEM(nx-1, :)) / dx;  % Backward difference for right boundary
Sy(:, 1) = (DEM(:, 2) - DEM(:, 1)) / dy;  % Forward difference for bottom boundary
Sy(:, ny) = (DEM(:, ny) - DEM(:, ny-1)) / dy;  % Backward difference for top boundary

end

function [vx, vy] = compute_boussinesq_velocities(h, K, dx, dy)
% COMPUTE_BOUSSINESQ_VELOCITIES - Computes groundwater velocities using Darcy's law
% where x = columns (horizontal) and y = rows (vertical).
%
% Inputs:
%   h  - Hydraulic head [m], size (ny, nx) => [rows, cols]
%   K  - Hydraulic conductivity [m/s], size (ny, nx)
%   dx - Grid spacing in x-direction [m] (columns)
%   dy - Grid spacing in y-direction [m] (rows)
%
% Outputs:
%   vx - Velocity in x-direction [m/s], size (ny, nx-1) — between columns
%   vy - Velocity in y-direction [m/s], size (ny-1, nx) — between rows

[ny, nx] = size(h);  % rows (Y), cols (X)

% Arithmetic mean of K at interfaces
Kx = 0.5 * (K(:, 1:nx-1) + K(:, 2:nx));   % Between columns
Ky = 0.5 * (K(1:ny-1, :) + K(2:ny, :));   % Between rows

% Head gradients (forward difference on staggered grid)
dhdx = (h(:, 2:nx) - h(:, 1:nx-1)) / dx;  % Gradient in x (columns)
dhdy = (h(2:ny, :) - h(1:ny-1, :)) / dy;  % Gradient in y (rows)

% Darcy velocity components
vx = -Kx .* dhdx;  % x-component of velocity (cols)
vy = -Ky .* dhdy;  % y-component of velocity (rows)

% Pad with NaNs for full-grid compatibility (optional)
vx(:, end+1) = nan;  % pad last column
vy(end+1, :) = nan;  % pad last row

end
