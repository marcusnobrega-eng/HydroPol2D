function test_pf_core_math()
%TEST_PF_CORE_MATH Minimal runnable checks for the PF math helpers.

case_dir = fileparts(mfilename('fullpath'));
addpath(case_dir);
rng(7);

w = pf_normalize_logweights([-1000; -1001; -1002]);
assert(abs(sum(w) - 1) < 1e-12);
assert(all(isfinite(w)));
assert(w(1) > w(2) && w(2) > w(3));

idx = pf_systematic_resample([0.34; 0.33; 0.33], 0.01);
assert(isequal(idx(:)', [1 2 3]));

spec = struct( ...
    'name', "soil1_ksat_surface_mult", ...
    'target_type', "SOIL", ...
    'class_id', 1, ...
    'property', "ksat_surface_mult", ...
    'baseline', 1, ...
    'truth', 1, ...
    'lower', 0.5, ...
    'upper', 1.5, ...
    'initial_std', 1, ...
    'perturb_std', 100, ...
    'transform', "log", ...
    'enabled', true);
theta = struct();
theta.soil(1).ksat_surface_mult = 1;
theta.lulc = struct();
theta.gw = struct();
theta = pf_perturb_theta(theta, spec, 1);
assert(theta.soil(1).ksat_surface_mult >= spec.lower);
assert(theta.soil(1).ksat_surface_mult <= spec.upper);
theta = pf_perturb_theta_inflated(theta, spec, struct('soil1_ksat_surface_mult', 10));
assert(theta.soil(1).ksat_surface_mult >= spec.lower);
assert(theta.soil(1).ksat_surface_mult <= spec.upper);

ic_spec = spec;
ic_spec.name = "ic1_initial_wtd_mult";
ic_spec.target_type = "IC";
ic_spec.property = "initial_wtd_mult";
theta.ic(1).initial_wtd_mult = 0.25;
theta2 = pf_perturb_theta_inflated(theta, ic_spec, struct('ic1_initial_wtd_mult', 10));
assert(theta2.ic(1).initial_wtd_mult == 0.25);

particles = repmat(struct('theta', [], 'state', [], 'weight', 1/3), 3, 1);
particles(1).theta.a = 1; particles(1).state.depths.d_t = 10;
particles(2).theta.a = 2; particles(2).state.depths.d_t = 20;
particles(3).theta.a = 3; particles(3).state.depths.d_t = 30;
copy_idx = [3; 1; 3];
resampled = particles(copy_idx);
assert(resampled(1).theta.a == 3);
assert(resampled(1).state.depths.d_t == 30);
assert(resampled(2).theta.a == 1);
assert(resampled(2).state.depths.d_t == 10);

forecast_state.flag = 123;
analysis_state = pf_apply_state_update_identity(forecast_state);
assert(analysis_state.flag == 123);

PF = pf_config_vtilted();
Domain = struct('valid', true(10, 10), 'cell_area_m2', 100);
obs = table( ...
    [1; 1], [15; 15], ["Q_outlet"; "GW_1"], ["discharge"; "groundwater_depth"], ...
    [1.0; 0.5], [0.001; 0.025], ["m3/s"; "m"], ...
    'VariableNames', {'window', 'time_min', 'obs_id', 'type', 'value', 'sigma', 'units'});
sim = [
    1.0, 0.5
    1.1, 0.8
    0.7, 0.2
];
[w, objective, logL, info] = pf_paper_weights(PF, Domain, obs, sim);
assert(abs(sum(w) - 1) < 1e-12);
assert(w(1) == max(w));
assert(objective(1) == min(objective));
assert(all(isfinite(logL)));
assert(info.method == "paper_beta_normal_pdf");

fprintf('test_pf_core_math passed.\n');
end
