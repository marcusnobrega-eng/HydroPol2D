% Author - Ambarish Prashant Chandurkar
% Adapted by Marcus Nobrega, Ph.D

function [F_forward,f] = GA_Newton_Raphson(F,dt,K,psi,theta,h,ir,precision,idx_imp)
% Input:
% F: infiltrated depth at time t
% dt: time-step in hours
% K: saturated hydraulic conductivity [mm/h]
% psi: soil suction head [mm]
% theta: soil moisture deficit
% h: ponding depth [mm]
% ir: Inflow_Rate in [mm/h]
% precision: number of decimal places
% idx_imp: logical mask representing impervious cells
%
% Outputs:
% F_forward: infiltrated depth at time t + dt
% f: infiltration rate at t = t + dt


% ---- Defining the derivative ---- %
% syms x F K dt psi theta h; % Infiltrated Depth
% f=-x + F + K.*dt + (h + psi).*theta.*(log(x + (h + psi).*theta) - log(F + (h + psi).*theta));
% g=diff(f); %The Derivative of the Function

% ---- Mannually Inputed Derivatives
f_x = @(x)(-x + F + K.*dt + (h + psi).*theta.*(log(x + (h + psi).*theta) - log(F + (h + psi).*theta)));
% g = (theta.*(h + psi))./(x + theta.*(h + psi)) - 1;

% Number of decimal plates
n = precision;
epsilon = 5*10^-(n+1);

% Initial Guess
x0 = max(F,0); % Infiltrated depth at time t
for i=1:100 % Number of maximum iterations
    % --- Calculating Symbolically --- %
    % f0 = vpa(subs(f,x,x0)); % Calculating the value of function at x0
    % f0_der = vpa(subs(g,x,x0)); %Calculating the value of function derivative at x0
    % --- Calculating Analytically --- %
    % f0 = vpa(subs(f,x,x0)); % Calculating the value of function at x0
    f0 = f_x(x0);
    f0_der = (theta.*(h + psi))./(x0 + theta.*(h + psi)) - 1;

    y = x0 -f0./f0_der; % The Formula

    err=abs(y-x0); % Error
    if sum(err < epsilon) == 0 % All values below the error threshold
        break
    end    

    x0 = y; % Updating
end

F_forward = max(y,0); % mm / h

% Infiltration Rate
f = 1/dt*(F_forward - F); % mm / h

% Cells which infiltrated depth exceeds the inflow rate
idx = f > ir;

f(idx) = ir(idx); % These cells had all inflow volume infiltrated

% Constraint at cells that have impervious surface
f(idx_imp) = 0;

% Changing infiltrated depth accordingly 
F_forward(idx) = F(idx) + f(idx)*dt; % mm
end