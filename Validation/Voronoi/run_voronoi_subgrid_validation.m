function summary = run_voronoi_subgrid_validation()
%RUN_VORONOI_SUBGRID_VALIDATION Checks for the HEC-RAS style sub-grid closure.
%
%   Mirrors the Python regression tests so the two models are checked against the
%   same statements. Each check records the measured symptom it guards.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','..','HydroPol2D_Functions'));
summary = struct();

% ---- 1. flat terrain: the two closures MUST agree -----------------------
% Level pool is exact on flat ground and every face profile is flat, so the
% sub-grid step must reproduce the flat-prism step to round-off. The Python
% implementation gives a discharge ratio of 1.0000 at three timesteps.
AREA = 2500; L = 50; DIST = 50; NROUGH = 0.03; BED = 0;
mesh = struct('n_cells',2,'n_edges',1, ...
    'edge_owner',[1],'edge_neighbor',[2], ...
    'edge_length',[L],'edge_distance',[DIST], ...
    'surface_bed',[BED;BED],'surface_area',[AREA;AREA]);
mesh.channel = struct('n_nodes',0);

hs = linspace(0,5,4001).';
tables = struct();
tables.cell_datum_m = [BED;BED];
tables.cell_zeta_m = [0 20; 0 20];
tables.cell_volume_m3 = [0 20*AREA; 0 20*AREA];
tables.cell_wet_area_m2 = [AREA AREA; AREA AREA];
tables.cell_point_count = [2;2];
tables.cell_plan_area_m2 = [AREA;AREA];
tables.cell_is_subgrid = [true;true];
tables.face_datum_m = BED;
tables.face_zeta_m = hs.';
tables.face_flow_area_m2 = (L*hs).';
tables.face_perimeter_m = repmat(L,1,numel(hs));
tables.face_conveyance = ((L*hs)/NROUGH .* max(hs,0).^(2/3)).';
tables.face_point_count = numel(hs);
tables.face_length_m = L;

ratios = zeros(3,1); dts = [0.5 1.0 2.0];
for k = 1:3
    v0 = [0.30;0.10] * AREA;
    [~, qa] = Voronoi_Local_Inertial_Step(mesh, v0, 0, [NROUGH;NROUGH], dts(k), ...
        critical_flow=true, dry_tolerance_m=1e-5);
    [~, qb] = Voronoi_Local_Inertial_Subgrid_Step(mesh, v0, 0, [NROUGH;NROUGH], ...
        tables, dts(k), critical_flow=true, dry_tolerance_m=1e-5);
    ratios(k) = qb(1) / qa(1);
    fprintf('  dt=%4.1fs  flat-prism Q %10.5f   sub-grid Q %10.5f   ratio %7.4f\n', ...
        dts(k), qa(1)*L, qb(1)*L, ratios(k));
end
summary.flat_agreement_max_error = max(abs(ratios - 1));
assert(summary.flat_agreement_max_error < 1e-9, ...
    'HydroPol2D:SubgridFlatMismatch', ...
    'On flat terrain the sub-grid and flat-prism closures must agree; worst ratio error %.3g.', ...
    summary.flat_agreement_max_error);

% ---- 2. cell stage inversion round-trips ---------------------------------
depths = [1e-6 0.01 0.1 1.0 10.0];
err = zeros(numel(depths),1);
for k = 1:numel(depths)
    vol = depths(k) * AREA * [1;1];
    st = hp2d_voronoi_subgrid_cell_stage(tables, vol);
    err(k) = max(abs(st - (BED + depths(k))));
end
summary.stage_roundtrip_max_error_m = max(err);
fprintf('  stage inversion round-trip: max error %.3e m\n', summary.stage_roundtrip_max_error_m);
assert(summary.stage_roundtrip_max_error_m < 1e-9, 'HydroPol2D:SubgridStageRoundtrip', ...
    'Cell stage inversion must round-trip; worst error %.3g m.', summary.stage_roundtrip_max_error_m);

% ---- 3. wetted perimeter is a STEP function ------------------------------
% Interpolating it linearly returned 12 m for a 20 m wide flat channel in Python.
[~, p_shallow, ~] = hp2d_voronoi_subgrid_face_state(tables, BED + 0.02);
summary.perimeter_at_shallow_depth_m = p_shallow;
fprintf('  wetted perimeter at 2 cm depth on a 50 m flat face: %.2f m (must be 50)\n', p_shallow);
assert(abs(p_shallow - L) < 1e-9, 'HydroPol2D:SubgridPerimeterInterpolated', ...
    'Wetted perimeter must be piecewise constant; got %.4f m for a %.0f m face.', p_shallow, L);

% ---- 4. selective closure: a face needs BOTH cells --------------------
tables_gated = tables;
tables_gated.cell_is_subgrid = [true;false];
v0 = [0.30;0.10] * AREA;
[~, ~, diag_gated] = Voronoi_Local_Inertial_Subgrid_Step(mesh, v0, 0, ...
    [NROUGH;NROUGH], tables_gated, 1.0, critical_flow=true, dry_tolerance_m=1e-5);
summary.gated_subgrid_faces = diag_gated.subgrid_faces;
fprintf('  one cell gated off: sub-grid faces %d (must be 0)\n', diag_gated.subgrid_faces);
assert(diag_gated.subgrid_faces == 0, 'HydroPol2D:SubgridFaceGate', ...
    'A face must use the tables only where BOTH its cells do.');

% ---- 5. the momentum step reports the PAIRED face depth ------------------
% critical_flow caps q at h*sqrt(g*h) during the step, so dividing the stored
% discharge by a depth recomputed afterwards is a pair the limiter never
% constrained: it put 31.56% of wet faces above critical, up to 182x, and reported
% 34.996 m/s where critical was 1.69 m/s.
[~, ~, diag_pair] = Voronoi_Local_Inertial_Step(mesh, v0, 0, [NROUGH;NROUGH], 1.0, ...
    critical_flow=true, dry_tolerance_m=1e-5);
summary.reports_paired_face_depth = isfield(diag_pair,'face_flow_depth_m') && ...
    numel(diag_pair.face_flow_depth_m) == mesh.n_edges;
fprintf('  momentum step reports face_flow_depth_m: %d\n', summary.reports_paired_face_depth);
assert(summary.reports_paired_face_depth, 'HydroPol2D:MissingPairedFaceDepth', ...
    'The momentum step must report the face depth paired with the stored discharge.');

% ---- 6. intensive remap preserves a constant on partial coverage ---------
% The conservative remap scaled a 605 m water surface to 302 m at a half-covered
% raster cell and to 0 m where the mesh did not reach.
mapping = struct();
mapping.mesh_to_raster = sparse((1:4).', (1:4).', [1;1;1;0.5], 4, 4);
mapping.raster_coverage = [1;1;1;0.5];
conservative = mapping.mesh_to_raster * repmat(605,4,1);
intensive = hp2d_remap_intensive(mapping, repmat(605,4,1), 0.25, NaN);
summary.conservative_min = min(conservative);
summary.intensive_max_error = max(abs(intensive - 605));
fprintf('  constant 605 m field: conservative min %.1f m, intensive max error %.3g m\n', ...
    summary.conservative_min, summary.intensive_max_error);
assert(summary.intensive_max_error < 1e-9, 'HydroPol2D:IntensiveRemap', ...
    'The intensive remap must preserve a constant field where the mesh reaches.');

% ---- 7. a dry domain must move no water --------------------------------
% The face TABLES still report geometry on dry ground; the SOLVER gates on water
% availability. That is deliberate: clipping the face profile up to a controlling
% invert instead flattens the face into a weir and inflated conveyance 57x in the
% Python implementation, so the gate belongs in the solver.
mesh_dry = mesh;
[v_dry, q_dry, dg_dry] = Voronoi_Local_Inertial_Subgrid_Step(mesh_dry, ...
    zeros(mesh_dry.n_cells,1), zeros(mesh_dry.n_edges,1), [NROUGH;NROUGH], ...
    tables, 1.0, critical_flow=true, dry_tolerance_m=1e-4);
summary.dry_max_discharge = max(abs(q_dry));
summary.dry_total_volume = sum(v_dry);
fprintf('  dry domain: max |q| %.3e, total volume %.3e (both must be 0)\n', ...
    summary.dry_max_discharge, summary.dry_total_volume);
assert(summary.dry_max_discharge == 0 && summary.dry_total_volume == 0, ...
    'HydroPol2D:SubgridDryDomain', 'A dry domain must move no water.');

% ---- 8. interior steps conserve mass -----------------------------------
vol = [0.30;0.10] * AREA; total0 = sum(vol); qe = 0;
for k = 1:20
    [vol, qe, ~] = Voronoi_Local_Inertial_Subgrid_Step(mesh, vol, qe, ...
        [NROUGH;NROUGH], tables, 0.5, critical_flow=true, dry_tolerance_m=1e-5);
end
summary.mass_relative_change = abs(sum(vol) - total0) / total0;
fprintf('  20 interior steps: mass relative change %.3e\n', summary.mass_relative_change);
assert(summary.mass_relative_change < 1e-12, 'HydroPol2D:SubgridMassBalance', ...
    'Interior sub-grid steps must conserve mass; relative change %.3g.', ...
    summary.mass_relative_change);

% ---- write the pass/fail record the registry auditor expects -------------
out_dir = fullfile(here, 'Outputs', 'Validation');
if ~isfolder(out_dir), mkdir(out_dir); end
checks = {
    'flat_terrain_agreement_with_flat_prism', summary.flat_agreement_max_error, 1e-9, 'max |ratio-1|'
    'cell_stage_inversion_roundtrip_m',       summary.stage_roundtrip_max_error_m, 1e-9, 'max |error| (m)'
    'wetted_perimeter_step_function_m',       abs(summary.perimeter_at_shallow_depth_m - L), 1e-9, 'deviation from face length (m)'
    'face_gate_requires_both_cells',          summary.gated_subgrid_faces, 0, 'sub-grid faces with one cell gated off'
    'momentum_reports_paired_face_depth',     double(~summary.reports_paired_face_depth), 0, '0 when reported'
    'intensive_remap_preserves_constant_m',   summary.intensive_max_error, 1e-9, 'max |error| (m)'
    'dry_domain_moves_no_water',              summary.dry_max_discharge, 0, 'max |q| (m2/s)'
    'interior_mass_conservation',             summary.mass_relative_change, 1e-12, 'relative mass change'
};
fid = fopen(fullfile(out_dir, 'Voronoi_Subgrid_Pass_Fail.csv'), 'w');
fprintf(fid, 'case_id,check,measured,threshold,units,passed\n');
all_pass = true;
for k = 1:size(checks,1)
    ok = checks{k,2} <= checks{k,3};
    all_pass = all_pass && ok;
    fprintf(fid, 'VAL-VORONOI-SUBGRID-001,%s,%.6g,%.6g,%s,%d\n', ...
        checks{k,1}, checks{k,2}, checks{k,3}, checks{k,4}, ok);
end
fprintf(fid, 'VAL-VORONOI-SUBGRID-001,passes_all_screening_criteria,%d,1,boolean,%d\n', ...
    all_pass, all_pass);
fclose(fid);
summary.pass_fail_file = fullfile(out_dir, 'Voronoi_Subgrid_Pass_Fail.csv');
summary.passes_all = all_pass;
fprintf('  wrote %s\n', summary.pass_fail_file);

fprintf('\nAll Voronoi sub-grid checks passed.\n');
end
