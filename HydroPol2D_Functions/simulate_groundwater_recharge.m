function [recharge_rate, updated_soil_moisture, cumulative_recharge] = simulate_groundwater_recharge( ...
    infiltration_rate, initial_soil_moisture, alpha, dt, min_soil_moisture, max_soil_moisture, idx_imp, current_recharge)

%% ═══════════════════════════════════════════════════════════════════════
%  Function: simulate_groundwater_recharge
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 03/06/2025
%  🔄 Updated: 02/21/2026 (HydroPol2D – vG/Mualem Darcy recharge)
% ─────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Simulate groundwater recharge from a 1-bucket unsaturated storage
%      using a reduced Richards-consistent Darcy flux approximation.
%
%      This version replaces the previous linear reservoir recharge model
%      with a physics-based closure using van Genuchten–Mualem functions:
%        • Pressure head: h(θ)
%        • Hydraulic conductivity: K(θ)
%      coupled to the water table via an effective drainage length scale.
%
%  ➤ Key Equations:
%      1) Bucket water content:
%           θ = θ_i + S / z_wt
%
%      2) van Genuchten pressure head:
%           h(θ) = -(1/α_vg) * (Se^(-1/m) - 1)^(1/n)
%           m = 1 - 1/n
%
%      3) Mualem-vG conductivity:
%           K(θ) = Ksat * Se^l * [1 - (1 - Se^(1/m))^m]^2
%
%      4) Recharge flux (downward positive):
%           R = K(θ) * ( 1 + h(θ)/L_gw )
%         where L_gw is an effective distance from the representative
%         bucket state to the water table (default: L_gw = 0.5*z_wt).
%
%      5) Storage update:
%           S_new = S_old + dt*(f - R)
%
%      Recharge is then recomputed from the update to enforce exact mass
%      balance:
%           R = f - (S_new - S_old)/dt
%
%  ➤ Inputs (units consistent with HydroPol2D coupling):
%      • infiltration_rate      - Infiltration flux into bucket [m/s]
%      • initial_soil_moisture  - Bucket storage (water depth) [m]
%      • alpha                  - (unused as linear coeff now)
%                                kept for backward compatibility; vG params
%                                must exist in workspace/structs (see below).
%      • dt                     - Time step [s]
%      • min_soil_moisture      - Minimum storage bound [m]
%      • max_soil_moisture      - Maximum storage bound (UZ capacity) [m]
%      • idx_imp                - Logical mask of impervious cells
%      • current_recharge       - Cumulative recharge memory [mm]
%
%  ➤ Outputs:
%      • recharge_rate          - Recharge flux to groundwater [m/s]
%      • updated_soil_moisture  - Updated bucket storage [m]
%      • cumulative_recharge    - Updated cumulative recharge [mm]
%
%  ➤ Required parameters (expected to be accessible here):
%      This function assumes the following fields exist in the caller
%      workspace as Soil_Properties (or equivalent globals you will wire):
%
%        Soil_Properties.theta_i
%        Soil_Properties.theta_sat
%        Soil_Properties.theta_r
%        Soil_Properties.ksat          [mm/h] (surface soil Ksat)
%        Soil_Properties.alpha_vg      [1/m]
%        Soil_Properties.n_vg          [-]
%        Soil_Properties.l_vg          [-]  (default 0.5)
%
%      If you prefer, you can later pass a Soil struct directly and remove
%      the dependency; for now we keep the interface unchanged.
% ═══════════════════════════════════════════════════════════════════════

% -------------------------------------------------------------------------
% 0) Trivial guards
% -------------------------------------------------------------------------
dt = max(dt, 1e-12);
infiltration_rate(isnan(infiltration_rate)) = 0;

% -------------------------------------------------------------------------
% 1) Impervious: no recharge dynamics (keep storage unchanged)
% -------------------------------------------------------------------------
updated_soil_moisture = initial_soil_moisture;
recharge_rate         = 0 * initial_soil_moisture;

% Only compute over pervious cells
idx = ~idx_imp;
if ~any(idx(:))
    cumulative_recharge = current_recharge; % [mm]
    return
end

% -------------------------------------------------------------------------
% 2) Compute effective unsaturated thickness (z_wt) from capacity
%    max_soil_moisture = z_wt * (theta_s - theta_i)
% -------------------------------------------------------------------------
% Pull soil parameters (you said you'll wire these in your structs)
% NOTE: this reads Soil_Properties from the caller workspace.
Soil_Properties = evalin('caller','Soil_Properties');

theta_i = Soil_Properties.theta_i;
theta_s = Soil_Properties.theta_sat;
theta_r = Soil_Properties.theta_r;

por = max(theta_s - theta_i, 1e-12);
zwt = max_soil_moisture ./ por;        % [m]
zwt = max(zwt, 0);

% Effective drainage length from representative bucket state to WT
Lgw = 0.5 * zwt;                       % [m]
Lgw = max(Lgw, 1e-6);

% -------------------------------------------------------------------------
% 3) Bucket storage -> theta (representative)
%    theta = theta_i + S/zwt (bounded)
% -------------------------------------------------------------------------
S = initial_soil_moisture;             % [m water]
zwt_safe = max(zwt, 1e-6);

theta = theta_i + (S ./ zwt_safe);
theta = max(theta, theta_i);
theta = min(theta, theta_s);

% Effective saturation
Se = (theta - theta_r) ./ max(theta_s - theta_r, 1e-12);
Se = min(max(Se, 1e-6), 1);

% -------------------------------------------------------------------------
% 4) van Genuchten pressure head h(theta) [m]
% -------------------------------------------------------------------------
n = Soil_Properties.n_vg;
m = 1 - 1./n;
a = Soil_Properties.alpha_vg;          % [1/m]

h_soil = -(1./a) .* ((Se.^(-1./m) - 1).^(1./n));
h_soil(Se >= 0.999999) = 0;            % saturated limit

% -------------------------------------------------------------------------
% 5) Mualem–vG K(theta) [m/s]
% -------------------------------------------------------------------------
if ~isfield(Soil_Properties,'l_vg') || isempty(Soil_Properties.l_vg)
    ell = 0.5;
else
    ell = Soil_Properties.l_vg;
end

Ks = Soil_Properties.ksat / 1000 / 3600;  % [mm/h] -> [m/s]

term = (1 - Se.^(1./m));
Kr = Se.^ell .* (1 - term.^m).^2;
Kr = min(max(Kr, 0), 1);

Ksoil = Ks .* Kr;                       % [m/s]

% -------------------------------------------------------------------------
% 6) Darcy-like recharge (downward positive), allows capillary rise
%    R = K(θ) * ( 1 + h(θ)/Lgw )
% -------------------------------------------------------------------------
grad = 1 + (h_soil ./ Lgw);

% Bound gradient to avoid absurd fluxes (stable but still physical)
% - Lower bound allows limited capillary rise (upward flux)
grad = max(grad, -1.0);
% - Upper bound prevents unrealistically strong drainage gradients
grad = min(grad,  2.0);

R = Ksoil .* grad;                      % [m/s]

% -------------------------------------------------------------------------
% 7) Update storage with mass conservation and bounds
% -------------------------------------------------------------------------
S_new = S + dt .* (infiltration_rate - R);   % [m]

% Enforce physical bounds
S_new = max(S_new, min_soil_moisture);
S_new = min(S_new, max_soil_moisture);

% Recompute recharge from exact mass balance:
%   dS/dt = f - R  ->  R = f - (S_new - S)/dt
R_mb = infiltration_rate - (S_new - S) ./ dt;

% Mask impervious and invalid
R_mb(~idx) = 0;

% Outputs
updated_soil_moisture = S_new;
recharge_rate         = R_mb;

% -------------------------------------------------------------------------
% 8) Cumulative recharge [mm]
% -------------------------------------------------------------------------
cumulative_recharge = current_recharge + recharge_rate * dt * 1000; % [mm]

end