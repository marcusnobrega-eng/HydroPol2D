function Results = compare_case1A_supply_limited_deep_soil()
%% ============================================================
% Case 1A: supply-limited infiltration, deep soil
%
% New infiltration closure:
%
%   Ksoil = wK*Ks + (1-wK)*K(theta_top)
%
% Original soil parameters:
%   Ks = 0.5 mm/h
%   wK = 0.5
%   Ltop = 0.05 m
%   dh_max = 1 m
%
% Gradient multiplier:
%   G = dh_max/Ltop + 1 = 21
%
% Expected capacity lower bound:
%   C_low = wK * Ks * G
%         = 0.5 * 0.5 * 21
%         = 5.25 mm/h
%
% Rainfall:
%   P = 5 mm/h
%
% Since P < C_low:
%   f = P = 5 mm/h
%   outlet Q should be approximately zero
% ============================================================

Cfg = struct();

Cfg.caseFolder = 'Case1A_supply_limited_deep_soil';

Cfg.flag_infiltration = 1;

Cfg.rain_mm_h = 5;
Cfg.rain_duration_min = 9 * 60;
Cfg.sim_duration_min  = 9 * 60;

Cfg.Lx_m = 2000;
Cfg.Ly_m = 1000;
Cfg.dx_m = 20;

Cfg.S0 = 0.01;
Cfg.n_manning = 0.015;

Cfg.zwt0_m = 5.0;

% Original soil 1 parameters
Cfg.theta_sat = 0.385;
Cfg.theta_r   = 0.068;
Cfg.theta_i   = 0.0997;

Cfg.ksat_mm_h = 0.5;
Cfg.Ltop_m = 0.05;
Cfg.dh_max_m = 1.0;
Cfg.Ksurf_weight = 0.5;

% For Case 1A, the new expected range collapses to f = rainfall.
Cfg.expected_f_low_mm_h  = 5.0;
Cfg.expected_f_high_mm_h = 5.0;

Cfg.expected_fill_time_h = NaN;

Results = compare_plane_case_core(Cfg);

end