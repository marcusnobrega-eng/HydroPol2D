function Results = compare_case1B_capacity_limited_deep_soil()
%% ============================================================
% Case 1B: capacity-limited infiltration, deep soil
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
% Expected capacity range:
%
%   C_low  = wK * Ks * G
%          = 0.5 * 0.5 * 21
%          = 5.25 mm/h
%
%   C_high = Ks * G
%          = 0.5 * 21
%          = 10.5 mm/h
%
% Rainfall:
%   P = 20 mm/h
%
% Therefore expected f should lie approximately between:
%   5.25 and 10.5 mm/h
%
% Expected outlet Q range:
%   using f = 10.5 -> Q = 5.278 m3/s
%   using f = 5.25 -> Q = 8.194 m3/s
%
% The model can fall anywhere in this range depending on theta_top.
% ============================================================

Cfg = struct();

Cfg.caseFolder = 'Case1B_capacity_limited_deep_soil';

Cfg.flag_infiltration = 1;

Cfg.rain_mm_h = 20;
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

% Expected range under the new Ksurf interface closure
Cfg.expected_f_low_mm_h  = 5.25;
Cfg.expected_f_high_mm_h = 10.5;

Cfg.expected_fill_time_h = NaN;

Results = compare_plane_case_core(Cfg);

end