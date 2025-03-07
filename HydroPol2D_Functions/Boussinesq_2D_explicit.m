function [h_t, u_x, u_y, q_exf, q_river] = Boussinesq_2D_explicit(dt,dx,dy,h,z0,Sy,R,K,river_mask,K_river,h_river,z_river,Courant,h_soil,catchment_mask,dirichlet_mask,h_dirichlet,perimeter)
% Solves the 2D Boussinesq equation using an explicit finite difference method
% with adaptive time-stepping based on the Courant stability condition.
%
% Inputs:
% - dt: Initial model time-step [sec]
% - dx: Grid spacing in x-direction [m]
% - dy: Grid spacing in y-direction [m]
% - h: Initial hydraulic head [m] (Nx, Ny)
% - z0: Ground elevation [m] (Nx, Ny)
% - Sy: Specific yield coefficient [1/m]
% - R: Recharge rate [m/s] (Nx, Ny)
% - K: Hydraulic conductivity [m/s] (Nx, Ny)
% - river_mask: logical mask indicating where rivers exist
% - K_river: hydraulic conductivity of the river [m/s]
% - h_river: hydraulic head at the rivers [matrix in DEM size]
% - z_river: topographic river elevtion [mm]
% - Courant: Stability factor (0 < Courant <= 1)
% - h_soil: Soil depth [m] (Nx, Ny)
% - catchment_mask: Binary mask indicating computational domain (Nx, Ny)
% - dirichlet_mask: Binary mask for Dirichlet boundary conditions (Nx, Ny)
% - h_dirichlet: Prescribed hydraulic head at Dirichlet boundaries [m] (Nx, Ny)
% - perimeter: Binary mask for no-flow boundaries (Nx, Ny)
%
% Outputs:
% - h_t: Final hydraulic head [m] (Nx, Ny)
% - u_x: Groundwater velocity in x-direction [m/s] (Nx-1, Ny)
% - u_y: Groundwater velocity in y-direction [m/s] (Nx, Ny-1)
% - q_exf: Exfiltration rate at cells [m/s] (Nx, Ny)
% - q_river: River exfiltrtion [m/s]. Positive leaves the aquifer. Negative,
% enters the aquifer
%

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
Nx = size(h,1); Ny = size(h,2);

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
    q_river = zeros(Nx, Ny);  % Initialize discharge array
    % river_cells = (river_mask == 1);  % Identify river locations
    % q_river(river_cells) = K_river(river_cells) .* ((h(river_cells) - h_river(river_cells)) ./ L_river); % m/s
    % 
    % % River water depth constraint
    % river_depth = h_river - z_river; % [m]
    % q_river(q_river<0 & river_depth <= 0) = 0; % No influx in the aquifer from zero-depth rivers
    % 
    % % River Maximum Outflow Constraint
    % idx = q_river < 0;
    % q_river(idx) = (-1)*max(abs(q_river(idx)), river_depth(idx)/dt);
    % % q_river(idx) = max(q_river(idx), - river_depth(idx)/dt);
    % 
    % if max(max(abs(q_river))) > 0 % At least one river cell has fluxes
    %     % Updated Groundwater Head
    %     h(river_cells) = h(river_cells) - dt./Sy(river_cells).*q_river(river_cells);
    % end
    
    
    % Apply no-flow boundary conditions
    % h_t = applyNoFlowBC(h_t, perimeter);

    % Compute second derivatives using central differencing
    depth_threshold = 0.01; % m
    [d2h_dx2, d2h_dy2] = compute_gradients(h_t, z0, K, dx, dy, depth_threshold);

    % Update hydraulic head using explicit finite difference
    h_t(2:end-1,:) = h_t(2:end-1,:) + dt./Sy(2:end-1,:) .* d2h_dx2; % x contribution
    h_t(:,2:end-1) = h_t(:,2:end-1) + dt./Sy(:,2:end-1) .* d2h_dy2; % y contribution

    % Apply Dirichlet boundary conditions
    h_t(dirichlet_mask) = h_dirichlet(dirichlet_mask); 

    % Recharge
    h_t = h_t + dt./Sy .* R; % Recharge contribution

    % Lost Mass
    lost_mass = nansum(nansum(dx*dy*(-1)*min(h_t-z0,0)));

    % Ensure groundwater head does not fall below ground level
    h_t = max(h_t, z0); % m

    % Compute exfiltration where groundwater head exceeds surface elevation 
    q_exf = max(0, Sy.*(h_t - h_surf) / dt); % m/s 

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
    error = nansum(nansum(cell_area*Sy.*(h_t - h_0))) - nansum(nansum(cell_area*dt.*R)) + nansum(nansum(cell_area*(dt*((q_river + q_exf)))));
    if error > 10
        ttt = 1;
    end

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


% % Compute the second derivatives (gradient) using matrix operations
% function [d2h_dx2, d2h_dy2] = compute_gradients(h, z0, K, dx, dy, depth_threshold)
% % Use finite difference with convolution for second derivatives
% 
% % Compute harmonic mean conductivities at interfaces
% Kx_p = mean_function(K(2:end-1, :), K(1:end-2, :),1,[]); % Right side conductivity
% Kx_m = mean_function(K(1:end-2, :), K(2:end-1, :),1,[]); % Left side conductivity
% hx_p = mean_function(h(2:end-1, :) - z0(2:end-1, :), h(1:end-2, :) - z0(1:end-2, :),2,depth_threshold); % Right side conductivity
% hx_m = mean_function(h(1:end-2, :) - z0(1:end-2, :), h(2:end-1, :) - z0(2:end-1, :),2,depth_threshold); % Left side conductivity
% 
% % Central difference for y-direction
% Ky_p = mean_function(K(:, 2:end-1), K(:, 1:end-2),1,[]); % Upward conductivity
% Ky_m = mean_function(K(:, 1:end-2), K(:, 2:end-1),1,[]); % Downward conductivity
% hy_p = mean_function(h(:, 2:end-1) - z0(:, 2:end-1), h(:, 1:end-2) - z0(:, 1:end-2),2,depth_threshold); % Upward conductivity
% hy_m = mean_function(h(:, 1:end-2) - z0(:, 1:end-2), h(:, 2:end-1) - z0(:, 2:end-1),2,depth_threshold); % Downward conductivity
% 
% % Compute second derivatives using correct finite differences
% Grad_1 = (h(3:end,:) - h(2:end-1,:))/dx; Grad_1(isnan(Grad_1)) = 0;
% Grad_2 = (h(2:end-1,:) - h(1:end-2,:))/dx; Grad_2(isnan(Grad_2)) = 0;
% Grad_3 = (h(:, 3:end) - h(:,2:end-1))/dy; Grad_3(isnan(Grad_3)) = 0;
% Grad_4 = (h(:,2:end-1) - h(:,1:end-2))/dy; Grad_4(isnan(Grad_4)) = 0;
% 
% d2h_dx2 = 1/dx*(Kx_p .*hx_p.*Grad_1 - Kx_m .*hx_m .* Grad_2);
% d2h_dx2(isnan(d2h_dx2)) = 0;
% 
% d2h_dy2 = 1/dy*(Ky_p .* hy_p.*Grad_3 - Ky_m .*hy_m .* Grad_4);
% d2h_dy2(isnan(d2h_dy2)) = 0;
% end

function [d2h_dx2, d2h_dy2] = compute_gradients(h, z0, K, dx, dy, depth_threshold)
    % Ensure non-negative depths

    H = max(h - z0, depth_threshold); 

    % Compute regularized harmonic mean for conductivities
    reg_eps = 1e-8; % Small value to avoid division by zero
    Kx_p = 2 * (K(2:end-1, :) .* K(1:end-2, :)) ./ (K(2:end-1, :) + K(1:end-2, :) + reg_eps);
    Kx_m = 2 * (K(1:end-2, :) .* K(2:end-1, :)) ./ (K(1:end-2, :) + K(2:end-1, :) + reg_eps);
    
    Ky_p = 2 * (K(:, 2:end-1) .* K(:, 1:end-2)) ./ (K(:, 2:end-1) + K(:, 1:end-2) + reg_eps);
    Ky_m = 2 * (K(:, 1:end-2) .* K(:, 2:end-1)) ./ (K(:, 1:end-2) + K(:, 2:end-1) + reg_eps);

    % Compute gradients with mixed central-upwind scheme
    hx_p = mean_function(H(2:end-1, :), H(1:end-2, :),2,depth_threshold); % Right side conductivity
    hx_m = mean_function(h(1:end-2, :), H(2:end-1, :),2,depth_threshold); % Left side conductivity
    
    hy_p = mean_function(H(:, 2:end-1), H(:, 1:end-2) ,2,depth_threshold); % Upward conductivity
    hy_m = mean_function(H(:, 1:end-2), H(:, 2:end-1),2,depth_threshold); % Downward conductivity

    % Apply hybrid central-upwind differencing
    Grad_x_fwd = (H(3:end,:) - H(2:end-1,:)) / dx;  
    Grad_x_bwd = (H(2:end-1,:) - H(1:end-2,:)) / dx;  
    Grad_y_fwd = (H(:, 3:end) - H(:, 2:end-1)) / dy;  
    Grad_y_bwd = (H(:, 2:end-1) - H(:, 1:end-2)) / dy;  

    % Use central difference where gradients are small, upwind where steep
    Grad_x = 0.5 * (Grad_x_fwd + Grad_x_bwd) .* (abs(Grad_x_fwd - Grad_x_bwd) < depth_threshold) ...
           + Grad_x_fwd .* (abs(Grad_x_fwd - Grad_x_bwd) >= depth_threshold & Grad_x_fwd > 0) ...
           + Grad_x_bwd .* (abs(Grad_x_fwd - Grad_x_bwd) >= depth_threshold & Grad_x_fwd <= 0);

    Grad_y = 0.5 * (Grad_y_fwd + Grad_y_bwd) .* (abs(Grad_y_fwd - Grad_y_bwd) < depth_threshold) ...
           + Grad_y_fwd .* (abs(Grad_y_fwd - Grad_y_bwd) >= depth_threshold & Grad_y_fwd > 0) ...
           + Grad_y_bwd .* (abs(Grad_y_fwd - Grad_y_bwd) >= depth_threshold & Grad_y_fwd <= 0);

    % Compute second derivatives
    d2h_dx2 = (Kx_p .* hx_p .* Grad_x - Kx_m .* hx_m .* Grad_x) / dx;
    d2h_dy2 = (Ky_p .* hy_p .* Grad_y - Ky_m .* hy_m .* Grad_y) / dy;

    % Ensure no NaN values
    d2h_dx2(isnan(d2h_dx2)) = 0;
    d2h_dy2(isnan(d2h_dy2)) = 0;
end

% Compute the maximum stable time step based on the Courant condition 
% (not quite correct) 
% function dt_max = compute_stable_dt(K, Sy, dx, dy, Courant)
%     dt_max = (Sy ./ K) .* ((dx^2 * dy^2) / (2 * (dx^2 + dy^2)));
%     dt_max = Courant * dt_max; % Apply a safety factor
%     dt_max(isinf(dt_max)) = nan;
%     dt_max = min(min(dt_max));
% end

function dt_max = compute_stable_dt(u_x, u_y, dx, dy, Courant)
    dt_max = min(Courant*dx./abs(u_x), Courant*dy./abs(u_y));
    dt_max = min(min(dt_max));
    dt_max(isinf(dt_max)) = nan;
end

function h = apply_free_flow_bc(h, mask)
[Nx, Ny] = size(h);

for i = 2:Nx-1
    for j = 2:Ny-1
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
    % This function modifies the groundwater table height matrix `h` to satisfy
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
% COMPUTE_BOUSSINESQ_VELOCITIES - Computes groundwater velocities from the Boussinesq equation
% using Darcy's law.
%
% Inputs:
% - h: Hydraulic head matrix [m] (nx, ny)
% - K: Hydraulic conductivity matrix [m/s] (nx, ny)
% - dx: Grid spacing in x-direction [m]
% - dy: Grid spacing in y-direction [m]
%
% Outputs:
% - vx: Groundwater velocity in x-direction [m/s] (nx-1, ny)
% - vy: Groundwater velocity in y-direction [m/s] (nx, ny-1)

[nx, ny] = size(h);

% Compute harmonic mean of hydraulic conductivity at cell interfaces
Kx = 2 ./ (1./K(1:nx-1, :) + 1./K(2:nx, :)); % Interface values in x-direction
Ky = 2 ./ (1./K(:, 1:ny-1) + 1./K(:, 2:ny)); % Interface values in y-direction

% Compute hydraulic head gradients using central differences
dhdx = (h(2:nx, :) - h(1:nx-1, :)) / dx;  % Gradient in x-direction
dhdy = (h(:, 2:ny) - h(:, 1:ny-1)) / dy;  % Gradient in y-direction

% Compute velocities using Darcy's law
vx = -Kx .* dhdx; % x-component of groundwater velocity
vy = -Ky .* dhdy; % y-component of groundwater velocity

% Include boundary faces
vx(end+1,:) = nan;
vy(:,end+1) = nan;

end