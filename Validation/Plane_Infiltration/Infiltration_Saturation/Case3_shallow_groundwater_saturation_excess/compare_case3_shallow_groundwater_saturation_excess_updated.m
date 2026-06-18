function Results = compare_case3_shallow_groundwater_saturation_excess()
%% ============================================================
% Case 3: near-surface groundwater saturation-excess test
%
% Purpose:
%   Evaluate whether a very shallow initial water table produces rapid
%   groundwater-driven saturation-excess overland flow.
%
% Updated Case 3 design:
%   - Soil/aquifer thickness: DTB = 0.30 m
%   - Initial water-table depth: zwt0 = 0.05 m below land surface
%   - Rainfall: 20 mm/h for 48 h
%   - Simulation duration: 48 h
%
% Storage convention:
%   HydroPol2D stores I_t as above-residual vadose storage.
%
%   Total UZ capacity       = zwt0 * (theta_sat - theta_r) * 1000
%   Initial UZ storage      = zwt0 * (theta_i   - theta_r) * 1000
%   Initial remaining store = zwt0 * (theta_sat - theta_i) * 1000
%
% Expected result:
%   The initial remaining storage is small (~14.27 mm). Infiltration should
%   be active early, then collapse rapidly once the near-surface groundwater
%   system fills the available unsaturated storage. Outlet discharge should
%   approach the no-infiltration benchmark for most of the 48 h simulation.
% ============================================================

Cfg = struct();

Cfg.caseFolder = 'Case3_shallow_groundwater_saturation_excess';

Cfg.flag_infiltration = 1;

% ------------------------------------------------------------
% Rainfall and simulation duration
% ------------------------------------------------------------
Cfg.rain_mm_h = 20;
Cfg.rain_duration_min = 48 * 60;
Cfg.sim_duration_min  = 48 * 60;

% ------------------------------------------------------------
% Plane geometry
% ------------------------------------------------------------
Cfg.Lx_m = 2000;
Cfg.Ly_m = 1000;
Cfg.dx_m = 20;

Cfg.S0 = 0.01;
Cfg.n_manning = 0.015;

% ------------------------------------------------------------
% Near-surface groundwater / thin aquifer setup
% ------------------------------------------------------------
Cfg.DTB_m  = 0.30;   % total soil/aquifer thickness [m]
Cfg.zwt0_m = 0.05;   % initial water-table depth below land surface [m]

% ------------------------------------------------------------
% Original soil 1 parameters
% ------------------------------------------------------------
Cfg.theta_sat = 0.385;
Cfg.theta_r   = 0.068;
Cfg.theta_i   = 0.0997;

Cfg.ksat_mm_h = 0.5;
Cfg.Ltop_m = 0.05;
Cfg.dh_max_m = 1.0;
Cfg.Ksurf_weight = 0.5;

% ------------------------------------------------------------
% Expected early-time infiltration range under current closure
% ------------------------------------------------------------
Cfg.expected_f_low_mm_h  = 5.25;
Cfg.expected_f_high_mm_h = 10.5;

% Use typical value from Case 1B for approximate marker.
Cfg.expected_fill_time_h = ...
    ((Cfg.theta_sat - Cfg.theta_i) * Cfg.zwt0_m * 1000) / 8.5;

% Extra diagnostic quantities used by compare_plane_case_core, if present.
Cfg.storage_capacity_total_above_residual_mm = ...
    Cfg.zwt0_m * (Cfg.theta_sat - Cfg.theta_r) * 1000;

Cfg.initial_storage_above_residual_mm = ...
    Cfg.zwt0_m * max(Cfg.theta_i - Cfg.theta_r, 0) * 1000;

Cfg.remaining_storage_initial_mm = ...
    Cfg.zwt0_m * (Cfg.theta_sat - Cfg.theta_i) * 1000;

Cfg.expected_note = ['Updated Case 3: near-surface groundwater, thin aquifer, ', ...
    '48 h rainfall. Infiltration should collapse after roughly 1.4--2.7 h ', ...
    'depending on the early infiltration rate, and outlet discharge should ', ...
    'approach the no-infiltration benchmark.'];

Results = compare_plane_case_core(Cfg);

end
