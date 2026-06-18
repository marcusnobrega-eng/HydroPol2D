function Results = compare_case2_shallow_soil_saturation_excess()
%% ============================================================
% Case 2: shallow soil saturation-excess test
%
% Purpose:
%   Check whether infiltration collapses when shallow soil storage fills.
%
% Soil:
%   Original soil 1
%
% Soil depth / initial WT depth:
%   zwt0 = 0.30 m
%
% Initial remaining storage:
%   S_rem0 = (theta_sat - theta_i) * zwt0 * 1000
%          = (0.385 - 0.0997) * 0.30 * 1000
%          = 85.59 mm
%
% New infiltration closure expected f range:
%   5.25 to 10.5 mm/h
%
% Expected fill-time range:
%   if f = 10.5 -> 85.59/10.5 = 8.15 h
%   if f = 5.25 -> 85.59/5.25 = 16.30 h
%
% From Case 1B, a typical f may be around 8.5 mm/h:
%   85.59 / 8.5 = 10.07 h
%
% Therefore in a 12 h run, the soil may fill near the end.
% ============================================================

Cfg = struct();

Cfg.caseFolder = 'Case2_shallow_soil_saturation_excess';

Cfg.flag_infiltration = 1;

Cfg.rain_mm_h = 20;
Cfg.rain_duration_min = 12 * 60;
Cfg.sim_duration_min  = 12 * 60;

Cfg.Lx_m = 2000;
Cfg.Ly_m = 1000;
Cfg.dx_m = 20;

Cfg.S0 = 0.01;
Cfg.n_manning = 0.015;

Cfg.zwt0_m = 0.30;

% Original soil 1 parameters
Cfg.theta_sat = 0.385;
Cfg.theta_r   = 0.068;
Cfg.theta_i   = 0.0997;

Cfg.ksat_mm_h = 0.5;
Cfg.Ltop_m = 0.05;
Cfg.dh_max_m = 1.0;
Cfg.Ksurf_weight = 0.5;

% Expected early-time range under new closure
Cfg.expected_f_low_mm_h  = 5.25;
Cfg.expected_f_high_mm_h = 10.5;

% Use typical value from Case 1B for approximate marker.
Cfg.expected_fill_time_h = ((Cfg.theta_sat - Cfg.theta_i) * Cfg.zwt0_m * 1000) / 8.5;

Results = compare_plane_case_core(Cfg);

end