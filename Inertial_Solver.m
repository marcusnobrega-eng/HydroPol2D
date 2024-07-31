function [q] = Inertial_Solver(flag_numerical_scheme,q_p,dt,Hf,S,n,dx)
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
%
% Output:
% q: discharge per face in [m2/s] at the end of the time-step

g = 9.81; % Gravity acceleration [m2/s]

%% Original Bates Formulation
if flag_numerical_scheme == 1
    q = (q_p - g*Hf*dt.*S)./ ...
        (1 + g*dt.*n.^2.*abs(q_p)./(Hf.^(7/3))); % m2 per sec (Hf can be simplified)
end

%% ------ s-centered scheme ----- %
if flag_numerical_scheme == 3

    %%% theta
    % theta = 0.8; % Controls artificial diffusivity
    theta = 1 - dt/dx*min(abs(q_p)./Hf,sqrt(g*Hf));
    q_average = 0*q_p;

    % x-x direction
    q_average(:,2:(end-1),1) = 1/2*(q_p(:,1:(end-2),1) + q_p(:,3:(end),1));        % Left
    q_average(:,2:(end-1),2) = 1/2*(q_p(:,1:(end-2),2) + q_p(:,3:(end),2));        % Right


    % y-y direction
    q_average(2:(end-1),:,3) = 1/2*(q_p(1:(end-2),:,3) + q_p(3:(end),:,3));        % Up
    q_average(2:(end-1),:,4) = 1/2*(q_p(1:(end-2),:,4) + q_p(3:(end),:,4));        % Down


    % total flux
    total_flux = q_p;

    % Flow
    q = (theta.*q_p + (1-theta).*(q_average) - g*Hf*dt.*S)./ ...
        (1 + g*dt.*n.^2.*abs(total_flux)./(Hf.^(7/3))); % m2 per sec (Hf can be simplified)
end

%% Upwind Scheme
if flag_numerical_scheme == 2
    %%% ---- Upwind Scheme ---- %%%
    q_flux = q_p;
    mask = q_flux < 0; % Negative fluxes

    %%% x-x direction
    q_flux(:,1:(end-1),1) = q_p(:,2:(end),1);        % Left
    q_flux(:,1:(end-1),2) = q_p(:,2:(end),2);        % Right


    %%% y-y direction
    q_flux(1:(end-1),:,3) = q_p(2:(end),:,3);        % Up
    q_flux(1:(end-1),:,4) = q_p(2:(end),:,4);        % Down

    q_flux(mask) = q_p(mask);

    total_flux = q_p;

    %%% Diffusivity factor
    theta = 1 - dt/dx*min(abs(q_p)./Hf,sqrt(g*Hf));

    %%% Flow
    q = (theta.*q_p + (1-theta).*(q_flux) - g*Hf*dt.*S)./ ...
        (1 + g*dt.*n.^2.*abs(total_flux)./(Hf.^(7/3))); % m2 per sec
end
end



