function q = Inertial_Solver(flag_numerical_scheme,q_p,dt,Hf,S,n_sq,dx,idx_nan)
% Calculates solution of Local Inertial Model
%
% Input:
% flag_numerical_scheme: chooses which numerical scheme is used
%   1: original Bates
%   2: s-upwind
%   3: s-centered
%
% q_p: discharge per face [m2/s]
% dt: time-step [s]
% Hf: effective flow depth at interface [m]
% S: water surface elevation slope [m/m]
% n_sq: Manning roughness coefficient squared
% dx: space discretization [m]
% idx_nan: cells outside domain
%
% Output:
% q: discharge per face [m2/s]

g = 9.81;

%% Original Bates formulation
if flag_numerical_scheme == 1

    q = (q_p - g*Hf*dt.*S) ./ ...
        (1 + g*dt.*n_sq.*abs(q_p)./(Hf.^(7/3)));

end

%% S-upwind scheme
if flag_numerical_scheme == 2

    total_flux = q_p;

    q_upwind = 0*q_p;

    mask_negative(:,:,1) = q_p(:,:,1) < 0;
    mask_negative(:,:,2) = q_p(:,:,2) < 0;

    ny = size(q_p,1);
    nx = size(q_p,2);

    mid_x = 2:(nx-1);
    mid_y = 2:(ny-1);

    % Positive flow direction
    q_upwind(:,mid_x,1) = q_p(:,mid_x-1,1);
    q_upwind(mid_y,:,2) = q_p(mid_y-1,:,2);

    % Negative flow direction
    q_upwind_negative = 0*q_p;

    q_upwind_negative(:,mid_x,1) = q_p(:,mid_x+1,1);
    q_upwind_negative(mid_y,:,2) = q_p(mid_y+1,:,2);

    % Select upwind value based on sign of q
    q_upwind(mask_negative) = q_upwind_negative(mask_negative);

    % Diffusivity factor
    theta = 1 - dt/dx*min(abs(q_p)./Hf, sqrt(g*Hf));

    mask = isnan(theta);
    theta(mask) = 1;

    % Flow
    q = (theta.*q_p + (1-theta).*q_upwind - g*Hf*dt.*S) ./ ...
        (1 + g*dt.*n_sq.*abs(total_flux)./(Hf.^(7/3)));

end

%% S-centered scheme
if flag_numerical_scheme == 3

    % Diffusivity factor
    theta = 1 - dt/dx*min(abs(q_p)./Hf, sqrt(g*Hf));

    mask = isnan(theta);
    theta(mask) = 1;

    q_average = 0*q_p;

    % x-direction centered average
    q_average(:,2:(end-1),1) = 0.5 * ...
        (q_p(:,1:(end-2),1) + q_p(:,3:end,1));

    % y-direction centered average
    q_average(2:(end-1),:,2) = 0.5 * ...
        (q_p(1:(end-2),:,2) + q_p(3:end,:,2));

    total_flux = q_p;

    % Flow
    q = (theta.*q_p + (1-theta).*q_average - g*Hf*dt.*S) ./ ...
        (1 + g*dt.*n_sq.*abs(total_flux)./(Hf.^(7/3)));

end

%% Domain cleanup
q(~isfinite(q)) = 0;

if ~isempty(idx_nan)
    q(logical(idx_nan)) = 0;
end

end