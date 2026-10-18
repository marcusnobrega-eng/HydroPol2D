function Summary = run_layered_step_parity(fixture_path)
%RUN_LAYERED_STEP_PARITY Replay the Python layered step (2 x 2 V-tilted) in MATLAB.
%   The fixture is written by HydroPol2D-Python tools/export_layered_step_fixture.py from
%   core.layered_vertical.advance_layered_subsurface. Each step applies
%   hp2d_layered_vertical_step, the water-table update and the exfiltration bookkeeping
%   of Layered_Subsurface_Module, and must match Python within 1e-10 relative.
arguments
    fixture_path (1,:) char = fullfile(fileparts(mfilename('fullpath')), 'layered_steps.json')
end
root = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(root, 'HydroPol2D_Functions'));
F = jsondecode(fileread(fixture_path));
names = {}; errors = [];
    function check(name, actual, expected)
        a = double(actual(:)); e = double(expected(:));
        scale = max(1, max(abs(e)));
        err = max(abs(a - e))/scale;
        names{end+1, 1} = name; errors(end+1, 1) = err; %#ok<AGROW>
        if ~(err < 1e-10)
            error('HydroPol2D:LayeredStepParity', '%s differs from Python: %.3e', name, err);
        end
    end
scenarios = F.scenarios;
if isstruct(scenarios), scenarios = num2cell(scenarios); end
for s = 1:numel(scenarios)
    C = scenarios{s};
    column = C.column;
    column.valid = logical(column.valid);
    Pj = C.parameters;
    P = struct('valid', logical(Pj.valid), 'pervious', Pj.pervious, ...
        'root_depth_m', Pj.root_depth_m, 'ksat_near_surface_mm_h', Pj.ksat_near_surface_mm_h, ...
        'l_vg', Pj.l_vg, 'ltop_m', Pj.ltop_m, 'dh_max_m', Pj.dh_max_m, ...
        'ksurf_weight', Pj.ksurf_weight, 'drainage_substeps', Pj.drainage_substeps, ...
        'infiltration', true, 'evapotranspiration', true, 'groundwater', true);
    if isempty(C.feddes)
        P.feddes = [];
    else
        P.feddes = C.feddes;
    end
    S.water_mm = C.initial.water_mm;
    S.head_m = C.initial.head_m;
    surface = zeros(size(S.head_m));
    steps = C.steps;
    if isstruct(steps), steps = num2cell(steps); end
    tag = sprintf('s%d', s);
    for i = 1:numel(steps)
        R = steps{i};
        S.surface_depth_mm = surface + C.rain_mm_h*C.dt_s/3600;
        check([tag '.surface_in'], S.surface_depth_mm, R.surface_in_mm);
        S.open_water_evaporation_mm = C.open_water_evaporation_mm;
        S.soil_et_demand_mm = C.soil_et_demand_mm;
        S = hp2d_layered_vertical_step(column, S, P, C.dt_s);
        [S.head_m, S.water_mm, exfiltration, unmet] = hp2d_water_table_update( ...
            column, S.head_m, S.water_mm, S.recharge_mm - S.capillary_mm);
        S.capillary_mm = S.capillary_mm - unmet;
        surface = S.surface_depth_mm + exfiltration.*P.valid;
        check([tag '.water'], S.water_mm, R.water_mm);
        check([tag '.head'], S.head_m, R.head_m);
        check([tag '.surface'], surface, R.surface_depth_mm);
        check([tag '.infiltration'], S.infiltration_mm, R.infiltration_mm);
        check([tag '.evaporation'], S.surface_evaporation_mm, R.surface_evaporation_mm);
        check([tag '.soil_et'], S.soil_et_mm, R.soil_et_mm);
        check([tag '.recharge'], S.recharge_mm, R.recharge_mm);
        check([tag '.capillary'], S.capillary_mm, R.capillary_mm);
        check([tag '.exfiltration'], exfiltration, R.exfiltration_mm);
    end
end
Summary = table(names, errors, 'VariableNames', {'check', 'max_relative_error'});
fprintf('Layered step parity: %d checks over %d scenarios, worst %.3e\n', height(Summary), ...
    numel(scenarios), max(errors));
end
