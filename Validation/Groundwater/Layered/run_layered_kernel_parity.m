function Summary = run_layered_kernel_parity(fixture_path)
%RUN_LAYERED_KERNEL_PARITY Compare the MATLAB layered-subsurface kernels with Python.
%   The fixture is written by HydroPol2D-Python tools/export_layered_fixtures.py
%   (Python -> MATLAB parity). Every kernel output must match within 1e-10 relative.
arguments
    fixture_path (1,:) char = fullfile(fileparts(mfilename('fullpath')), 'layered_kernels.json')
end
root = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(root, 'HydroPol2D_Functions'));
F = jsondecode(fileread(fixture_path));
names = {}; errors = [];
    function check(name, actual, expected)
        a = double(actual); e = double(expected);
        if ~isequal(size(a), size(e))
            error('HydroPol2D:LayeredParity', '%s: size %s vs %s', name, mat2str(size(a)), mat2str(size(e)));
        end
        both_nan = isnan(a) & isnan(e);
        scale = max(1, max(abs(e(~isnan(e))), [], 'all'));
        if isempty(scale), scale = 1; end
        diff = abs(a - e); diff(both_nan) = 0;
        err = max(diff(:))/scale;
        if isempty(err), err = 0; end
        names{end+1, 1} = name; errors(end+1, 1) = err; %#ok<AGROW>
        if ~(err < 1e-10)
            error('HydroPol2D:LayeredParity', '%s differs from Python: %.3e', name, err);
        end
    end
    function C = col(S)
        C = S;
        C.valid = logical(S.valid);
    end

for i = 1:numel(F.deep_layers)
    D = F.deep_layers(i);
    check(sprintf('deep_layers_%d', i), hp2d_deep_layer_thicknesses(D.max_depth_m, D.first_m, ...
        D.ratio, D.max_m, D.explicit), D.out(:).');
end

I = F.column_inputs;
dtb = min(I.depth_to_bedrock_m, I.max_depth_m);
[z_top, z_bot] = hp2d_slot_geometry(I.root_depth_m, min(I.soil_base_depth_m, dtb), dtb, ...
    I.near_surface_depth_m, I.deep_thicknesses_m);
P = I; P.column_base_depth_m = dtb;
column = hp2d_slot_column(z_top, z_bot, P);
C = col(F.column);
for f = {'z_top_m', 'z_bot_m', 'theta_s', 'theta_r', 'alpha_m_inverse', 'n_vg', 'kv_m_s', ...
        'kh_m_s', 'theta_retained', 'ground_m', 'base_depth_m'}
    check(['column.' f{1}], column.(f{1}), C.(f{1}));
end
check('column.valid', column.valid, C.valid);

for i = 1:numel(F.diagnostics)
    D = F.diagnostics(i);
    check('depth', hp2d_water_table_depth(column, D.head_m), D.depth_m);
    check('unsaturated', hp2d_unsaturated_thickness(column, D.head_m), D.unsaturated_m);
    check('saturated', hp2d_saturated_thickness(column, D.head_m), D.saturated_m);
    check('theta', hp2d_slot_theta(column, D.head_m, D.water_mm), D.theta);
    check('storage', hp2d_column_storage_mm(column, D.head_m, D.water_mm), D.storage_mm);
    check('drainable', hp2d_drainable_storage_mm(column, D.head_m), D.drainable_mm);
    check('specific_yield', hp2d_effective_specific_yield(column, D.head_m), D.specific_yield);
    check('transmissivity', hp2d_layered_transmissivity(column, D.head_m), D.transmissivity_m2_s);
end
for i = 1:numel(F.water_table_update)
    U = F.water_table_update(i);
    [h, w, exf, unmet] = hp2d_water_table_update(column, U.head_m, U.water_mm, U.net_mm);
    check('wtu.head', h, U.out_head_m); check('wtu.water', w, U.out_water_mm);
    check('wtu.exfiltration', exf, U.exfiltration_mm); check('wtu.unmet', unmet, U.unmet_removal_mm);
end
for i = 1:numel(F.drain_cascade)
    D = F.drain_cascade(i);
    check(sprintf('drain_%d', D.substeps), hp2d_drain_cascade(column, D.head_m, D.water_mm, ...
        D.dt_s, D.substeps), D.out_water_mm);
end
for i = 1:numel(F.water_table_flux)
    D = F.water_table_flux(i);
    [w, recharge, capillary] = hp2d_water_table_flux(column, D.head_m, D.water_mm, D.dt_s);
    check('flux.water', w, D.out_water_mm); check('flux.recharge', recharge, D.recharge_mm);
    check('flux.capillary', capillary, D.capillary_rise_mm);
end
S = F.initial_slot_water;
check('initial_water', hp2d_initial_slot_water(column, S.head_m, S.soil_theta_initial), S.out_water_mm);

Fd = F.feddes;
fp = struct('h1_m', Fd.h1_m, 'h2_m', Fd.h2_m, 'h3_high_m', Fd.h3_high_m, ...
    'h3_low_m', Fd.h3_low_m, 'h4_m', Fd.h4_m);
check('feddes.alpha', hp2d_feddes_alpha(Fd.head_m, fp, Fd.rate_mm_day), Fd.alpha);
check('feddes.h3', hp2d_feddes_h3(fp, Fd.rate_mm_day), Fd.h3_m);
check('vg.head', hp2d_matric_head_from_theta(Fd.theta, Fd.theta_r, Fd.theta_s, Fd.alpha_vg, Fd.n_vg), Fd.matric_head_m);
check('vg.theta', hp2d_theta_from_matric_head(Fd.head_m, Fd.theta_r, Fd.theta_s, Fd.alpha_vg, Fd.n_vg), Fd.theta_from_head);

for i = 1:numel(F.shared_table)
    S = F.shared_table(i);
    sub = col(S.column);
    R = hp2d_solve_shared_water_table(sub, S.cell_head_m, S.water_mm, S.net_volume_m3, S.host, ...
        S.area_m2, S.ground_m, S.base_depth_m);
    check('shared.head', R.head_m, S.out_cell_head_m); check('shared.water', R.water_mm, S.out_water_mm);
    check('shared.exfiltration', R.exfiltration_m3, S.out_exfiltration_m3);
    check('shared.unmet', R.unmet_removal_m3, S.out_unmet_removal_m3);
    M = F.move_to_depth(i);
    [delta, w] = hp2d_storage_change_to_depth(col(M.column), M.head_m, M.water_mm, M.depth_m);
    check('move.delta', delta, M.out_delta_mm); check('move.water', w, M.out_water_mm);
end

Lt = F.lateral;
lat_column = col(Lt.column);
R = hp2d_layered_lateral_advance(lat_column, reshape(Lt.head_m, Lt.ny, Lt.nx), Lt.water_mm, ...
    Lt.dt_s, Lt.dx_m, Lt.dy_m, reshape(logical(Lt.fixed), Lt.ny, Lt.nx), ...
    reshape(Lt.fixed_head_m, Lt.ny, Lt.nx), Lt.courant);
check('lateral.head', R.head_m(:), Lt.out_head_m);
check('lateral.water', R.water_mm, Lt.out_water_mm);
check('lateral.exfiltration', R.exfiltration_mm(:), Lt.out_exfiltration_mm);
check('lateral.boundary', R.boundary_inflow_m3, Lt.out_boundary_inflow_m3);
check('lateral.substeps', R.substeps, Lt.out_substeps);

Summary = table(names, errors, 'VariableNames', {'check', 'max_relative_error'});
fprintf('Layered kernel parity: %d checks, worst %.3e\n', height(Summary), max(errors));
end
