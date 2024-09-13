function [q] = Inertial_Solver(flag_numerical_scheme,q_p,dt,Hf,S,n,dx,idx_nan)
% Calculates solution of Local Inertial Model
%
% Input:
% flag_numerical_scheme: chooses which numerical scheme is used
% 1: original bates
% 2: s-upwind
% 3: s-centered
% q_p: discharge per face in [m2/s]
% dt: time-step in seconds
% Hf: effective flow depth in the interface [m]
% S: water surface elevation slope [m/m] at
% n: manning's roughness coefficient [sm-1/3]
% dx: space discretizaion [m]
% idx_nan: cells outside of the domain
%
% Output:
% q: discharge per face in [m2/s] at the end of the time-step

g = 9.81; % Gravity acceleration [m2/s]

%% Original Bates Formulation
if flag_numerical_scheme == 1
    q = (q_p - g*Hf*dt.*S)./ ...
        (1 + g*dt.*n.^2.*abs(q_p)./(Hf.^(7/3))); % m2 per sec (Hf can be simplified)
end

%% Upwind Scheme
if flag_numerical_scheme == 2
    total_flux = q_p;
    q_upwind = 0*q_p;
    mask_negative(:,:,1) = q_p(:,:,1)  < 0;
    mask_negative(:,:,2) = q_p(:,:,2)  < 0; % This one might be positve
    nx = size(n,2);
    ny = size(n,1);
    mid_x = 2:(nx-1);
    mid_y = 2:(ny-1);

    %%% Positive
    % x-x
    % q_upwind(:,mid_x,1) = q_p(:,mid_x-1,1); % Left
    q_upwind(:,mid_x,1) = q_p(:,mid_x-1,1); % Right

    % y-y
    q_upwind(mid_y,:,2) = q_p(mid_y-1,:,2); % Up
    % q_upwind(mid_y,:,4) = q_p(mid_y-1,:,4); % Down

    %%% Negative
    q_upwind_negative = 0*q_p;
    % q_upwind_negative(:,mid_x,1) = q_p(:,mid_x+1,1); % Left
    q_upwind_negative(:,mid_x,1) = q_p(:,mid_x+1,1); % Right

    % y-y
    q_upwind_negative(mid_y,:,2) = q_p(mid_y+1,:,2); % Up
    % q_upwind_negative(mid_y,:,4) = q_p(mid_y+1,:,4); % Down

    % Final q_upwind
    q_upwind(mask_negative) = q_upwind_negative(mask_negative);

    %%% Diffusivity factor
    theta = 1 - dt/dx*min(abs(q_p)./Hf,sqrt(g*Hf));

    %%% Effect of borders
    % theta(1,:,1:2) = 1; theta(size(q_p,1),:,1:2) = 1;
    % theta(:,1,1:2) = 1; theta(:,size(q_p,2),1:2) = 1;
    % %
    % theta(:,2,1) = 1;
    % theta(:,end-1,1) = 1;
    % theta(2,:,2) = 1;
    % theta(end-1,:,2) = 1;
    % mask = isnan(theta);
    % theta = max(theta,0.7);
    % theta(mask) = 1;

    %%% Flow
    q = (theta.*q_p + (1-theta).*(q_upwind) - g*Hf*dt.*S)./ ...
        (1 + g*dt.*n.^2.*abs(total_flux)./(Hf.^(7/3))); % m2 per sec
end

%% ------ s-centered scheme ----- %
if flag_numerical_scheme == 3

    %%% theta
    % theta = 0.8; % Controls artificial diffusivity
    theta = 1 - dt/dx*min(abs(q_p)./Hf,sqrt(g*Hf));
    mask = isnan(theta);
    % theta = max(theta,0.7);
    theta(mask) = 1;

    q_average = 0*q_p;

    % x-x direction
    % q_average(:,2:(end-1),1) = 1/2*(q_p(:,1:(end-2),1) + q_p(:,3:(end),1));        % Left
    q_average(:,2:(end-1),1) = 1/2*(q_p(:,1:(end-2),1) + q_p(:,3:(end),1));          % Right


    % y-y direction
    q_average(2:(end-1),:,2) = 1/2*(q_p(1:(end-2),:,2) + q_p(3:(end),:,2));          % Up
    % q_average(2:(end-1),:,4) = 1/2*(q_p(1:(end-2),:,4) + q_p(3:(end),:,4));        % Down


    % total flux
    total_flux = q_p;
    % total_flux = sqrt((1/2*(q_p(:,:,2) + q_p(:,:,1))).^2 + (1/2*(q_p(:,:,3) - q_p(:,:,4))).^2);

    %%% Effect of borders
    % theta(1,:,1:2) = 1; theta(size(q_p,1),:,1:2) = 1;
    % theta(:,1,1:2) = 1; theta(:,size(q_p,2),1:2) = 1;
    % %
    % theta(:,2,1) = 1;
    % theta(:,end-1,1) = 1;
    % theta(2,:,2) = 1;
    % theta(end-1,:,2) = 1;
    % mask = isnan(theta);
    % theta = max(theta,0.7);
    % theta(mask) = 1;

    % Flow
    q = (theta.*q_p + (1-theta).*(q_average) - g*Hf*dt.*S)./ ...
        (1 + g*dt.*n.^2.*abs(total_flux)./(Hf.^(7/3))); % m2 per sec (Hf can be simplified)
end

end



