function Results = compare_case0_no_infiltration()
%% ============================================================
% Case 0: no-infiltration hydraulic control
% ============================================================

Cfg = struct();

Cfg.caseFolder = 'Case0_no_infiltration';

Cfg.flag_infiltration = 0;

Cfg.rain_mm_h = 20;
Cfg.rain_duration_min = 9 * 60;
Cfg.sim_duration_min  = 9 * 60;

Cfg.Lx_m = 2000;
Cfg.Ly_m = 1000;
Cfg.dx_m = 20;

Cfg.S0 = 0.01;
Cfg.n_manning = 0.015;

Cfg.zwt0_m = 5.0;
Cfg.theta_sat = 0.385;
Cfg.theta_r = 0.068;
Cfg.theta_i = 0.0997;

Cfg.expected_f_mm_h = 0;

Results = compare_plane_case_core(Cfg);

end