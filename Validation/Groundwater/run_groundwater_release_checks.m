function Results = run_groundwater_release_checks()
%RUN_GROUNDWATER_RELEASE_CHECKS Run every groundwater, layered, Feddes and deep-aquifer check.
%   Returns a table (check, passed, detail, seconds) and prints it. Includes the
%   Python -> MATLAB parity harnesses, the deep-aquifer checks and seven full-model
%   Excel-mode V-tilted runs (about 1 min each). A full run must close the volume
%   balance within 1e-6 m3; the deep runs must also start at S0 = q0 x area x tau.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(here, fullfile(here, 'Layered'), fullfile(root, 'HydroPol2D_Functions'));
checks = {
    'deep-aquifer step', @() run_deep_aquifer_validation()
    'deep-aquifer initial state', @() run_deep_aquifer_initial_state_validation()
    'layered kernel parity', @() run_layered_kernel_parity()
    'layered step parity', @() run_layered_step_parity()
    'Feddes legacy consistency', @() run_feddes_legacy_consistency()
    'layered mesh case', @() run_layered_mesh_case()
    };
variants = {'local', 'lateral', 'feddes', 'legacy', 'deep', 'deep_legacy', 'deep_empty'};
for k = 1:numel(variants)
    checks(end + 1, :) = {['V-tilted ' variants{k}], @() vtilted_check(variants{k})}; %#ok<AGROW>
end
checks(end + 1, :) = {'V-tilted lateral K multipliers', @() multiplier_check()};
names = checks(:, 1);
passed = false(numel(names), 1);
detail = strings(numel(names), 1);
seconds = zeros(numel(names), 1);
for k = 1:numel(names)
    started = tic;
    try
        out = checks{k, 2}();
        passed(k) = true;
        detail(k) = describe(out);
    catch problem
        detail(k) = string(problem.message);
    end
    seconds(k) = toc(started);
    fprintf('[%s] %s (%.0f s): %s\n', ternary(passed(k), 'PASS', 'FAIL'), names{k}, ...
        seconds(k), detail(k));
end
Results = table(string(names), passed, detail, seconds, ...
    'VariableNames', {'check', 'passed', 'detail', 'seconds'});
fprintf('Groundwater release checks: %d of %d passed\n', nnz(passed), numel(passed));
end

function out = vtilted_check(variant)
R = run_layered_vtilted_case(variant);
assert(abs(R.volume_error_m3) < 1e-6 && all(abs(R.module_errors_m3) < 1e-6), ...
    'Volume error %.3e m3, module errors %s m3.', R.volume_error_m3, mat2str(R.module_errors_m3, 3));
out = sprintf('volume error %.2e m3', R.volume_error_m3);
if any(strcmp(variant, {'deep', 'deep_legacy'}))
    expected = 0.5/1000*R.deep_active_area_m2*100;   % q0 = 0.5 mm/day, tau = 100 d
    S0 = sum(R.deep_initial_storage_m3);
    assert(abs(S0 - expected) <= 1e-9*expected, 'Deep S0 %.6g m3, expected %.6g m3.', S0, expected);
    out = sprintf('%s; deep S0 %.6g m3 = q0 x area x tau', out, S0);
elseif strcmp(variant, 'deep_empty')
    assert(isempty(R.deep_initial_storage_m3) || all(R.deep_initial_storage_m3 == 0), ...
        'An empty start must have zero deep storage.');
end
end

function out = multiplier_check()
% Soil x2 and regolith x0.5 lateral-conductivity multipliers scale exactly those slots' Kh.
base = run_layered_vtilted_case('lateral');
scaled = run_layered_vtilted_case('lateral_scaled');
assert(abs(scaled.volume_error_m3) < 1e-6, 'Volume error %.3e m3.', scaled.volume_error_m3);
ratio = scaled.kh_m_s./base.kh_m_s;
soil = ratio(:, 1:3); rock = ratio(:, 4:end);
assert(all(abs(soil(isfinite(soil)) - 2) < 1e-12) && all(abs(rock(isfinite(rock)) - 0.5) < 1e-12), ...
    'Kh multipliers were not applied to the soil and deep slots.');
out = sprintf('soil Kh x2, deep Kh x0.5; volume error %.2e m3', scaled.volume_error_m3);
end

function text = describe(out)
if ischar(out) || isstring(out)
    text = string(out);
elseif istable(out)
    text = sprintf('%d rows', height(out));
elseif isstruct(out)
    f = fieldnames(out);
    text = strjoin(cellfun(@(n) sprintf('%s=%s', n, mat2str(out.(n), 3)), f(1:min(3, end))', ...
        'UniformOutput', false), ', ');
else
    text = "ok";
end
end

function value = ternary(condition, a, b)
if condition, value = a; else, value = b; end
end
