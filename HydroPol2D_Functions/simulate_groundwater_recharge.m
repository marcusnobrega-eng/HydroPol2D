function [recharge_rate, updated_soil_moisture, cumulative_recharge] = simulate_groundwater_recharge( ...
    infiltration_rate, initial_soil_moisture, alpha, dt, min_soil_moisture, max_soil_moisture, idx_imp, current_recharge)

%% ═══════════════════════════════════════════════════════════════════════
%  Function: simulate_groundwater_recharge
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 03/06/2025
%  🔄 Updated: 03/06/2026 (HydroPol2D – conservative vG/Mualem Darcy recharge)
% ─────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Simulate groundwater recharge from a 1-bucket unsaturated storage
%      using a reduced Richards-consistent Darcy flux approximation with
%      strict storage conservation under a moving water table.
%
%      Key modeling choices in this conservative update:
%        • Impervious cells may still drain/recharge from stored vadose
%          water; only infiltration should be zero there upstream.
%        • The maximum vadose storage is allowed to vary with water-table
%          depth through max_soil_moisture.
%        • If the water table rises and vadose capacity shrinks, storage is
%          clipped to the new admissible capacity. This is treated as a
%          bookkeeping transfer from unsaturated to saturated storage, not
%          as a surface-water source.
%        • Recharge returned by this routine is recomputed from the final
%          storage update to guarantee exact bucket mass balance:
%
%              dS/dt = f - R
%              R = f - (S_new - S_old)/dt
%
%  ➤ Inputs:
%      • infiltration_rate      - Infiltration flux into bucket [m/s]
%      • initial_soil_moisture  - Bucket storage above θ_i [m]
%      • alpha                  - Unused; kept for backward compatibility
%      • dt                     - Time step [s]
%      • min_soil_moisture      - Lower storage bound [m]
%      • max_soil_moisture      - Upper storage bound / UZ capacity [m]
%      • idx_imp                - Impervious mask (not used to zero R)
%      • current_recharge       - Cumulative recharge memory [mm]
%
%  ➤ Outputs:
%      • recharge_rate          - Recharge flux to groundwater [m/s]
%      • updated_soil_moisture  - Updated bucket storage [m]
%      • cumulative_recharge    - Updated cumulative recharge [mm]
%
%  ➤ Required parameters in caller workspace:
%        Soil_Properties.theta_i
%        Soil_Properties.theta_sat
%        Soil_Properties.theta_r
%        Soil_Properties.ksat          [mm/h]
%        Soil_Properties.alpha_vg      [1/m]
%        Soil_Properties.n_vg          [-]
%        Soil_Properties.l_vg          [-]   optional, default 0.5
% ═══════════════════════════════════════════════════════════════════════

% Keep backward compatibility
alpha = alpha; %#ok<NASGU>

%% -----------------------------------------------------------------------
% 0) Guards and initialization
% ------------------------------------------------------------------------
dt = max(dt, 1e-12);

% Pull soil parameters from caller workspace
Soil_Properties = evalin('caller','Soil_Properties');

theta_i = Soil_Properties.theta_i;
theta_s = Soil_Properties.theta_sat;
theta_r = Soil_Properties.theta_r;

n = Soil_Properties.n_vg;
m = 1 - 1 ./ n;
a = Soil_Properties.alpha_vg;                          % [1/m]

if ~isfield(Soil_Properties,'l_vg') || isempty(Soil_Properties.l_vg)
    ell = 0.5;
else
    ell = Soil_Properties.l_vg;
end

Ks = Soil_Properties.ksat / 1000 / 3600;              % [m/s]

% Initialize outputs
recharge_rate         = zeros(size(initial_soil_moisture), 'like', initial_soil_moisture);
updated_soil_moisture = initial_soil_moisture;
cumulative_recharge   = current_recharge;

% Sanitize inputs
infiltration_rate(~isfinite(infiltration_rate))         = 0;
initial_soil_moisture(~isfinite(initial_soil_moisture)) = 0;
min_soil_moisture(~isfinite(min_soil_moisture))         = 0;
max_soil_moisture(~isfinite(max_soil_moisture))         = 0;

% Effective admissible lower bound cannot exceed upper bound
min_eff = min(min_soil_moisture, max_soil_moisture);

% Active cells: do not suppress drainage/recharge on impervious cells
idx = isfinite(initial_soil_moisture) & isfinite(max_soil_moisture) & isfinite(min_eff);
if ~any(idx(:))
    return
end

%% -----------------------------------------------------------------------
% 1) Enforce physically admissible initial storage for current WT
%    If WT rose since previous step, max_soil_moisture may have shrunk.
%    Any excess above the new capacity is removed from UZ bookkeeping here.
% ------------------------------------------------------------------------
S = initial_soil_moisture;                             % [m]
S(~idx) = 0;

S = max(S, min_eff);
S = min(S, max_soil_moisture);

%% -----------------------------------------------------------------------
% 2) Compute effective unsaturated thickness from current capacity
%    max_soil_moisture = z_wt * (theta_s - theta_i)
% ------------------------------------------------------------------------
por_eff = max(theta_s - theta_i, 1e-12);
zwt = max_soil_moisture ./ por_eff;                    % [m]
zwt = max(zwt, 0);

% Effective drainage length
Lgw = max(0.5 .* zwt, 1e-6);                           % [m]

%% -----------------------------------------------------------------------
% 3) Bucket storage -> representative water content
%    Storage S is defined above theta_i over the current unsaturated depth
% ------------------------------------------------------------------------
zwt_safe = max(zwt, 1e-6);

theta = theta_i + S ./ zwt_safe;
theta = max(theta, theta_i);
theta = min(theta, theta_s);

Se = (theta - theta_r) ./ max(theta_s - theta_r, 1e-12);
Se = min(max(Se, 1e-6), 1);

%% -----------------------------------------------------------------------
% 4) van Genuchten pressure head h(theta) [m]
% ------------------------------------------------------------------------
h_soil = -(1 ./ a) .* ((Se .^ (-1 ./ m) - 1) .^ (1 ./ n));
h_soil(Se >= 0.999999) = 0;

%% -----------------------------------------------------------------------
% 5) Mualem–vG conductivity K(theta) [m/s]
% ------------------------------------------------------------------------
term  = 1 - Se .^ (1 ./ m);
Kr    = Se .^ ell .* (1 - term .^ m) .^ 2;
Kr    = min(max(Kr, 0), 1);

Ksoil = Ks .* Kr;                                      % [m/s]

%% -----------------------------------------------------------------------
% 6) Darcy-like recharge (downward positive), allows capillary rise
%    R = K(theta) * (1 + h(theta)/Lgw)
% ------------------------------------------------------------------------
grad = 1 + h_soil ./ Lgw;

% Stable bounds
grad = max(grad, -1.0);    % limited upward flux
grad = min(grad,  2.0);    % limited strong drainage

R = Ksoil .* grad;                                      % [m/s]

% Optional: if there is no unsaturated thickness, there is no UZ bucket.
% In that case, the bucket storage is zero by definition. We set the bucket
% exchange returned by this routine to zero and leave any upward/downward
% exchange with the surface/saturated domain to the groundwater solver.
idx_no_uz = (max_soil_moisture <= 0);
R(idx_no_uz) = 0;

%% -----------------------------------------------------------------------
% 7) Conservative bucket update with current WT-based storage bounds
% ------------------------------------------------------------------------
S_new = S + dt .* (infiltration_rate - R);             % [m]

% Enforce admissible storage bounds under current WT
S_new = max(S_new, min_eff);
S_new = min(S_new, max_soil_moisture);

% If WT is at surface, UZ storage must be zero
S_new(idx_no_uz) = 0;

% Exact mass-balance recharge implied by final storage update
R_mb = infiltration_rate - (S_new - S) ./ dt;          % [m/s]

% In cells with no UZ, do not pass a bucket recharge flux
R_mb(idx_no_uz) = 0;

% Mask invalid cells
R_mb(~idx)  = 0;
S_new(~idx) = initial_soil_moisture(~idx);

%% -----------------------------------------------------------------------
% 8) Outputs
% ------------------------------------------------------------------------
updated_soil_moisture = S_new;
recharge_rate         = R_mb;

% Cumulative recharge [mm]
cumulative_recharge = current_recharge + max(recharge_rate,0) .* dt .* 1000;

end