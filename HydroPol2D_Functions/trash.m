% Initialize parameters
Nx = 100;   % Number of grid points in the x-direction
Ny = 100;   % Number of grid points in the y-direction
dx = 1;     % Grid spacing in the x-direction
dy = 1;     % Grid spacing in the y-direction
dt = 0.1;   % Time step
Nt = 1000;  % Number of time steps
g = 9.81;   % Gravity (m/s^2)

% Create grid
x = (1:Nx) * dx;    % x-grid
y = (1:Ny) * dy;    % y-grid

% Initial conditions (water depth and velocity)
h = ones(Nx, Ny);   % Initial water depth (flat surface)
u = zeros(Nx, Ny);  % Initial velocity in x-direction
v = zeros(Nx, Ny);  % Initial velocity in y-direction

% Manning's roughness coefficient (example, should be spatially varying)
n = 0.03 * ones(Nx, Ny);

% Friction slopes (calculated from the terrain)
Sfx = zeros(Nx, Ny);  % Friction slope in x-direction (example)
Sfy = zeros(Nx, Ny);  % Friction slope in y-direction (example)

% Construct matrices A_x and B_y
A_x = construct_matrix_A(Nx, Ny, dx, dt);
B_y = construct_matrix_B(Nx, Ny, dy, dt, g, n, h, u, v);

% Time-stepping loop
for t = 1:Nt
    % Step 1: Solve in the x-direction (momentum equation for u)
    rhs_x = calculate_rhs_x(h, u, v, dx, dt, g);  % Right-hand side for x-direction
    u = A_x \ rhs_x;  % Solve for u in the x-direction (implicit)

    % Step 2: Solve in the y-direction (momentum equation for v)
    rhs_y = calculate_rhs_y(h, u, v, dy, dt, g);  % Right-hand side for y-direction
    v = B_y \ rhs_y;  % Solve for v in the y-direction (implicit)

    % Step 3: Update the water depth using continuity equation
    h = update_h(h, u, v, dt, dx, dy);  % Update water depth based on the continuity equation
    
    % Plot or store results
    if mod(t, 100) == 0
        figure(1);
        surf(x, y, h);
        title(['Water Depth at t = ', num2str(t*dt)]);
        colorbar;
        pause(0.1);
    end
end

%% Constructing Matrix A_x (Implicit in x-direction)
function A_x = construct_matrix_A(Nx, Ny, dx, dt)
    % Create the sparse matrix for the implicit system in the x-direction (velocity)
    
    % Number of unknowns (Nx * Ny)
    N = Nx * Ny;
    
    % Main diagonal and off-diagonals for the system
    main_diag = 1 + 2 * dt / dx^2 * ones(N, 1);  % Main diagonal
    left_diag = -dt / dx^2 * ones(N - 1, 1);     % Left diagonal (off by 1)
    right_diag = -dt / dx^2 * ones(N - 1, 1);    % Right diagonal (off by 1)
    
    % Create the sparse matrix A_x using sparse matrix construction
    A_x = spdiags([left_diag, main_diag, right_diag], [-1, 0, 1], N, N);
end

%% Constructing Matrix B_y (Implicit in y-direction)
function B_y = construct_matrix_B(Nx, Ny, dy, dt)
    % Create the sparse matrix for the implicit system in the y-direction (velocity)
    
    % Number of unknowns (Nx * Ny)
    N = Nx * Ny;
    
    % Main diagonal and off-diagonals for the system
    main_diag = 1 + 2 * dt / dy^2 * ones(N, 1);  % Main diagonal
    down_diag = -dt / dy^2 * ones(N - Ny, 1);   % Downward diagonal (off by Ny)
    up_diag = -dt / dy^2 * ones(N - Ny, 1);     % Upward diagonal (off by Ny)
    
    % Create the sparse matrix B_y using sparse matrix construction
    B_y = spdiags([down_diag, main_diag, up_diag], [-Ny, 0, Ny], N, N);
end

%% Right-hand side for x-direction (momentum equation for u)
function rhs_x = calculate_rhs_x(h, u, v, dx, dt, g)
    % Compute the right-hand side for the x-direction momentum equation
    rhs_x = zeros(size(h));  % Initialize right-hand side (same size as h)
    
    % Calculate the spatial derivatives of q_x = h * u
    q_x = h .* u;
    
    % Central difference approximation for the spatial derivative
    rhs_x(2:end-1, :) = -dt / dx^2 * (q_x(3:end, :) - q_x(1:end-2, :));
end

%% Right-hand side for y-direction (momentum equation for v)
function rhs_y = calculate_rhs_y(h, u, v, dy, dt, g)
    % Compute the right-hand side for the y-direction momentum equation
    rhs_y = zeros(size(h));  % Initialize right-hand side (same size as h)
    
    % Calculate the spatial derivatives of q_y = h * v
    q_y = h .* v;
    
    % Central difference approximation for the spatial derivative
    rhs_y(:, 2:end-1) = -dt / dy^2 * (q_y(:, 3:end) - q_y(:, 1:end-2));
end

%% Update water depth using continuity equation
function h_new = update_h(h, u, v, dt, dx, dy)
    % Update the water depth using the continuity equation
    % Central difference for spatial derivatives of q_x and q_y
    q_x = h .* u;
    q_y = h .* v;
    
    h_new = h - dt * ( (q_x(2:end, :) - q_x(1:end-1, :)) / dx + (q_y(:, 2:end) - q_y(:, 1:end-1)) / dy );
end
