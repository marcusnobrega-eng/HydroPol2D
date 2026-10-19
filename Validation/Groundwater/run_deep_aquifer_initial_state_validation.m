function report = run_deep_aquifer_initial_state_validation()
%RUN_DEEP_AQUIFER_INITIAL_STATE_VALIDATION Deep-aquifer starting state S0 = q0 x area x tau.
%   Checks hp2d_deep_aquifer_initial_state against the fixture shared with
%   HydroPol2D-Python (tests/fixtures/deep_aquifer_initial_state.json), that a start
%   from the mean recharge is a steady state of hp2d_deep_aquifer_step, the per-cell
%   area form, and that giving both starting values is rejected.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(root, 'HydroPol2D_Functions'));
fixture = jsondecode(fileread(fullfile(here, 'deep_aquifer_initial_state.json')));
ids = fixture.subbasin_id;            % null (outside the domain) decodes to NaN
area = fixture.cell_area_m2;
tau = fixture.response_time_days;
f = fixture.recharge_fraction;
worst = 0;
cases = fixture.cases;
if ~iscell(cases), cases = num2cell(cases); end   % differing fields decode to a cell array
for k = 1:numel(cases)
    c = cases{k};
    q0 = NaN; recharge = NaN;
    if isfield(c, 'initial_baseflow_mm_day'), q0 = c.initial_baseflow_mm_day; end
    if isfield(c, 'mean_recharge_mm_year'), recharge = c.mean_recharge_mm_year; end
    [S, Q, sorted_ids] = hp2d_deep_aquifer_initial_state(ids, area, tau, f, q0, recharge);
    assert(isequal(sorted_ids, [1; 2]));
    rel = max([abs(S - c.expected_storage_m3)./c.expected_storage_m3; ...
        abs(Q - c.expected_discharge_m3_s)./c.expected_discharge_m3_s]);
    assert(rel < 1e-12, 'Fixture case %s differs from Python by %g.', c.name, rel);
    worst = max(worst, rel);
end

% Per-cell areas give the same result as the scalar area.
[S_scalar] = hp2d_deep_aquifer_initial_state(ids, area, tau, f, 0.5, NaN);
[S_cells] = hp2d_deep_aquifer_initial_state(ids, area + zeros(size(ids)), tau, f, 0.5, NaN);
assert(max(abs(S_scalar - S_cells)) < 1e-12);

% No starting value: empty reservoirs.
assert(all(hp2d_deep_aquifer_initial_state(ids, area, tau, f, NaN, NaN) == 0));

% A start from the mean recharge stays put under that recharge (steady state).
recharge_mm_year = 300;
[S0, Q0] = hp2d_deep_aquifer_initial_state(ids, area, tau, f, NaN, recharge_mm_year);
dt_s = 86400;
active = isfinite(ids);
recharge_rate = zeros(size(ids));
recharge_rate(active) = recharge_mm_year/1000/(365.25*86400);
state = struct();
[state, ~, ~] = hp2d_deep_aquifer_step(state, zeros(size(ids)), ids, double(active), ...
    f, tau, area, dt_s);              % first call builds the subbasin index
state.deep_aquifer_storage_m3 = S0;
[state, ~, surface_mm] = hp2d_deep_aquifer_step(state, recharge_rate, ids, double(active), ...
    f, tau, area, dt_s);
steady = max(abs(state.deep_aquifer_storage_m3 - S0)./S0);
outflow = abs(sum(surface_mm(active))/1000*area/dt_s - sum(Q0))/sum(Q0);
assert(steady < 1e-12 && outflow < 1e-10, 'Mean-recharge start is not steady (%g, %g).', ...
    steady, outflow);

% Both starting values at once are rejected.
rejected = false;
try
    hp2d_deep_aquifer_initial_state(ids, area, tau, f, 0.5, 300);
catch
    rejected = true;
end
assert(rejected, 'Giving both starting values must raise an error.');

report = struct('fixture_worst_relative_error', worst, 'steady_state_drift', steady, ...
    'steady_outflow_error', outflow);
fprintf('Deep-aquifer initial state: fixture worst %.2e, steady drift %.2e, outflow %.2e\n', ...
    worst, steady, outflow);
end
