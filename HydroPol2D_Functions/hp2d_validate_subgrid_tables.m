function summary = hp2d_validate_subgrid_tables(SubgridTables)
%HP2D_VALIDATE_SUBGRID_TABLES Validate required lookup-table fields.

if isfield(SubgridTables, 'sfincs_exact') && SubgridTables.sfincs_exact
    summary = validate_sfincs_exact_tables(SubgridTables);
    return;
end

required = {'volume_cell','area_cell','invert_el', ...
    'area_x','width_x','wetfrac_x','hrep_x','phi_x','perimeter_x','Rh_x','n_x','nrep_x','K_x','invert_x', ...
    'area_y','width_y','wetfrac_y','hrep_y','phi_y','perimeter_y','Rh_y','n_y','nrep_y','K_y','invert_y', ...
    'area_north','width_north','hrep_north','phi_north','nrep_north','invert_north', ...
    'area_south','width_south','hrep_south','phi_south','nrep_south','invert_south', ...
    'area_west','width_west','hrep_west','phi_west','nrep_west','invert_west', ...
    'area_east','width_east','hrep_east','phi_east','nrep_east','invert_east', ...
    'dz','maxDepth'};

missing = {};
for k = 1:numel(required)
    if ~isfield(SubgridTables, required{k}) || isempty(SubgridTables.(required{k}))
        missing{end+1} = required{k}; %#ok<AGROW>
    end
end

if ~isempty(missing)
    error('hp2d_validate_subgrid_tables:missingFields', ...
        'Missing required subgrid table fields: %s', strjoin(missing, ', '));
end

summary = struct();
summary.volume_monotonic_violations = count_monotonic_violations(SubgridTables.volume_cell);
summary.area_negative_count = nnz(SubgridTables.area_cell(:) < -eps(class_underlying_like(SubgridTables.area_cell)));
summary.face_area_negative_count = nnz(SubgridTables.area_x(:) < -eps(class_underlying_like(SubgridTables.area_x))) + ...
    nnz(SubgridTables.area_y(:) < -eps(class_underlying_like(SubgridTables.area_y)));
summary.face_width_negative_count = nnz(SubgridTables.width_x(:) < -eps(class_underlying_like(SubgridTables.width_x))) + ...
    nnz(SubgridTables.width_y(:) < -eps(class_underlying_like(SubgridTables.width_y)));
summary.face_hrep_negative_count = nnz(SubgridTables.hrep_x(:) < -eps(class_underlying_like(SubgridTables.hrep_x))) + ...
    nnz(SubgridTables.hrep_y(:) < -eps(class_underlying_like(SubgridTables.hrep_y)));
summary.face_phi_out_of_range_count = nnz(SubgridTables.phi_x(:) < -eps(class_underlying_like(SubgridTables.phi_x)) | SubgridTables.phi_x(:) > 1 + eps(class_underlying_like(SubgridTables.phi_x))) + ...
    nnz(SubgridTables.phi_y(:) < -eps(class_underlying_like(SubgridTables.phi_y)) | SubgridTables.phi_y(:) > 1 + eps(class_underlying_like(SubgridTables.phi_y)));
summary.face_perimeter_negative_count = nnz(SubgridTables.perimeter_x(:) < -eps(class_underlying_like(SubgridTables.perimeter_x))) + ...
    nnz(SubgridTables.perimeter_y(:) < -eps(class_underlying_like(SubgridTables.perimeter_y)));
summary.face_Rh_negative_count = nnz(SubgridTables.Rh_x(:) < -eps(class_underlying_like(SubgridTables.Rh_x))) + ...
    nnz(SubgridTables.Rh_y(:) < -eps(class_underlying_like(SubgridTables.Rh_y)));
summary.face_conveyance_negative_count = nnz(SubgridTables.K_x(:) < -eps(class_underlying_like(SubgridTables.K_x))) + ...
    nnz(SubgridTables.K_y(:) < -eps(class_underlying_like(SubgridTables.K_y)));
summary.bad_roughness_count = nnz(isfinite(SubgridTables.n_x(:)) & SubgridTables.n_x(:) <= 0) + ...
    nnz(isfinite(SubgridTables.n_y(:)) & SubgridTables.n_y(:) <= 0);

if summary.volume_monotonic_violations > 0
    error('hp2d_validate_subgrid_tables:nonmonotonicVolume', ...
        'Subgrid volume table has %d monotonicity violations.', summary.volume_monotonic_violations);
end

if summary.area_negative_count > 0 || summary.face_area_negative_count > 0 || ...
        summary.face_width_negative_count > 0 || summary.face_hrep_negative_count > 0 || ...
        summary.face_phi_out_of_range_count > 0 || summary.face_perimeter_negative_count > 0 || ...
        summary.face_Rh_negative_count > 0 || summary.face_conveyance_negative_count > 0
    error('hp2d_validate_subgrid_tables:negativeGeometry', ...
        'Subgrid tables contain invalid area, hrep, phi, perimeter, roughness, or conveyance values.');
end

wet_x = SubgridTables.area_x > 0;
wet_y = SubgridTables.area_y > 0;
bad_geometry_x = wet_x & (SubgridTables.width_x <= 0 | SubgridTables.hrep_x <= 0 | ...
    SubgridTables.phi_x <= 0 | SubgridTables.Rh_x <= 0 | SubgridTables.K_x <= 0 | ...
    SubgridTables.nrep_x <= 0);
bad_geometry_y = wet_y & (SubgridTables.width_y <= 0 | SubgridTables.hrep_y <= 0 | ...
    SubgridTables.phi_y <= 0 | SubgridTables.Rh_y <= 0 | SubgridTables.K_y <= 0 | ...
    SubgridTables.nrep_y <= 0);
summary.bad_wet_face_geometry_count = nnz(bad_geometry_x(:)) + nnz(bad_geometry_y(:));
if summary.bad_wet_face_geometry_count > 0
    error('hp2d_validate_subgrid_tables:badWetFaceGeometry', ...
        'Subgrid tables contain %d wet face entries with missing width, hydraulic radius, or conveyance.', ...
        summary.bad_wet_face_geometry_count);
end
end

function summary = validate_sfincs_exact_tables(SubgridTables)
required = {'z_zmin','z_zmax','z_volmax','z_level', ...
    'u_zmin','u_zmax','u_havg','u_nrep','u_pwet','u_navg','u_ffit', ...
    'v_zmin','v_zmax','v_havg','v_nrep','v_pwet','v_navg','v_ffit'};

missing = {};
for k = 1:numel(required)
    if ~isfield(SubgridTables, required{k}) || isempty(SubgridTables.(required{k}))
        missing{end+1} = required{k}; %#ok<AGROW>
    end
end

if ~isempty(missing)
    error('hp2d_validate_subgrid_tables:missingSfincsFields', ...
        'Missing required SFINCS subgrid table fields: %s', strjoin(missing, ', '));
end

summary = struct();
summary.z_level_monotonic_violations = count_monotonic_violations(SubgridTables.z_level);
summary.bad_cell_bounds = nnz(isfinite(SubgridTables.z_zmin(:)) & ...
    isfinite(SubgridTables.z_zmax(:)) & SubgridTables.z_zmax(:) < SubgridTables.z_zmin(:));
summary.bad_u_phi = nnz(isfinite(SubgridTables.u_pwet(:)) & ...
    (SubgridTables.u_pwet(:) < -eps(class_underlying_like(SubgridTables.u_pwet)) | ...
     SubgridTables.u_pwet(:) > 1 + eps(class_underlying_like(SubgridTables.u_pwet))));
summary.bad_v_phi = nnz(isfinite(SubgridTables.v_pwet(:)) & ...
    (SubgridTables.v_pwet(:) < -eps(class_underlying_like(SubgridTables.v_pwet)) | ...
     SubgridTables.v_pwet(:) > 1 + eps(class_underlying_like(SubgridTables.v_pwet))));
summary.bad_u_havg = nnz(isfinite(SubgridTables.u_havg(:)) & SubgridTables.u_havg(:) < -eps(class_underlying_like(SubgridTables.u_havg)));
summary.bad_v_havg = nnz(isfinite(SubgridTables.v_havg(:)) & SubgridTables.v_havg(:) < -eps(class_underlying_like(SubgridTables.v_havg)));
summary.bad_u_nrep = nnz(isfinite(SubgridTables.u_nrep(:)) & SubgridTables.u_nrep(:) <= 0);
summary.bad_v_nrep = nnz(isfinite(SubgridTables.v_nrep(:)) & SubgridTables.v_nrep(:) <= 0);

if summary.z_level_monotonic_violations > 0 || summary.bad_cell_bounds > 0 || ...
        summary.bad_u_phi > 0 || summary.bad_v_phi > 0 || ...
        summary.bad_u_havg > 0 || summary.bad_v_havg > 0 || ...
        summary.bad_u_nrep > 0 || summary.bad_v_nrep > 0
    error('hp2d_validate_subgrid_tables:badSfincsTables', ...
        'SFINCS subgrid tables contain invalid z-level, H_G, phi, or nrep values.');
end
end

function nbad = count_monotonic_violations(V)
dV = diff(V, 1, 3);
tol = 100 * eps(class_underlying_like(V));
nbad = nnz(isfinite(dV(:)) & dV(:) < -tol);
end

function c = class_underlying_like(x)
if isa(x, 'gpuArray')
    c = classUnderlying(x);
else
    c = class(x);
end
end
