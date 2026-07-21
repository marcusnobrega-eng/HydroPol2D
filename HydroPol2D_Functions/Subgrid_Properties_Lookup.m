function [SubgridTables, invert_el] = Subgrid_Properties_Lookup( ...
    DEM_raster, Roughness_raster, Reference_raster, coarse_res, varargin)
%SUBGRID_PROPERTIES_LOOKUP Build SFINCS-style lookup-subgrid tables.
%
% This routine follows the subgrid corrections for the linear inertial
% equations described by van Ormondt et al. (2025). The production solver
% uses the SFINCS fields below:
%
%   Cell continuity points:
%       z_zmin, z_zmax, z_volmax, z_dep
%
%   Velocity points:
%       u_zmin, u_zmax, u_havg, u_nrep, u_pwet, u_navg_w, u_ffit
%       v_zmin, v_zmax, v_havg, v_nrep, v_pwet, v_navg_w, v_ffit
%
% Legacy HydroPol fields are also populated as aliases for diagnostics and
% older plotting routines, but lookup-subgrid routing should query the
% SFINCS fields directly.

p = inputParser;
p.addParameter('nr_levels', 10, @(x)isscalar(x) && x >= 2);
p.addParameter('nlevels', [], @(x)isempty(x) || (isscalar(x) && x >= 2));
p.addParameter('tol_rel', 5e-3, @(x)isscalar(x) && x > 0);
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('huthresh', 0.01, @(x)isscalar(x) && x >= 0);
p.addParameter('zvolmin', -20.0, @(x)isscalar(x) && isfinite(x));
p.addParameter('max_gradient', 99999.0, @(x)isscalar(x) && x > 0);
p.addParameter('q_table_option', 2, @(x)isscalar(x) && any(x == [1 2]));
p.addParameter('weight_option', 'min', @(x)ischar(x) || isstring(x));
% Accepted only for backward compatibility with earlier HydroPol configs.
p.addParameter('dz', [], @(x)isempty(x) || (isscalar(x) && x > 0));
p.addParameter('maxDepth', [], @(x)isempty(x) || (isscalar(x) && x > 0));
p.parse(varargin{:});

if ~isempty(p.Results.nlevels)
    nr_levels = round(p.Results.nlevels);
else
    nr_levels = round(p.Results.nr_levels);
end
tol_rel = p.Results.tol_rel;
verbose = logical(p.Results.verbose);
huthresh = double(p.Results.huthresh);
zvolmin = double(p.Results.zvolmin);
max_gradient = double(p.Results.max_gradient);
q_table_option = double(p.Results.q_table_option);
weight_option = lower(char(p.Results.weight_option));
if ~ismember(weight_option, {'min', 'mean'})
    error('Subgrid_Properties_Lookup:badWeightOption', ...
        'weight_option must be ''min'' or ''mean''.');
end

DEM = double(DEM_raster.Z);
cellsize = double(DEM_raster.cellsize);
ROUGH_INPUT = double(Roughness_raster.Z);
ROUGH_INPUT(~isfinite(ROUGH_INPUT) | ROUGH_INPUT <= 0) = NaN;

roughness_is_fine = isequal(size(ROUGH_INPUT), size(DEM));
roughness_is_coarse = isequal(size(ROUGH_INPUT), size(Reference_raster.Z));
if ~roughness_is_fine && ~roughness_is_coarse
    error(['Roughness_raster.Z must either match DEM_raster.Z for fine subgrid ', ...
           'roughness or Reference_raster.Z for coarse fallback roughness.']);
end
if roughness_is_coarse
    warning('Subgrid_Properties_Lookup:coarseRoughnessFallback', ...
        ['Using coarse Manning values as a fallback for subgrid roughness. ', ...
         'Fine-resolution Manning is preferred for SFINCS-style tables.']);
end

r = coarse_res / cellsize;
r_round = round(r);
rel_err = abs(r - r_round) / max(r, eps);
if rel_err > tol_rel
    warning('Subgrid_Properties_Lookup:nonIntegerRatio', ...
        ['coarse_res/cellsize is not approximately integer (r = %.6f). ', ...
         'SFINCS-style subgrid tables require an exact or near-exact ', ...
         'coarse-to-fine ratio.'], r);
end
if mod(r_round, 2) ~= 0
    error('Subgrid_Properties_Lookup:oddRefinementRatio', ...
        ['SFINCS-style velocity-point tables require an even integer ', ...
         'coarse-to-fine refinement ratio. Got r = %d.'], r_round);
end

[nrows, ncols] = size(DEM);
nrows_coarse = size(Reference_raster.Z, 1);
ncols_coarse = size(Reference_raster.Z, 2);
if nrows ~= nrows_coarse * r_round || ncols ~= ncols_coarse * r_round
    error('Subgrid_Properties_Lookup:gridMismatch', ...
        ['Exact SFINCS lookup tables require the fine DEM dimensions to equal ', ...
         'Reference_raster dimensions times the integer refinement ratio. ', ...
         'Got fine %dx%d, coarse %dx%d, ratio %d.'], ...
        nrows, ncols, nrows_coarse, ncols_coarse, r_round);
end

SubgridTables = struct();
SubgridTables.sfincs_exact = true;
SubgridTables.nr_levels = nr_levels;
SubgridTables.table_huthresh = huthresh;
% SFINCS resets runtime huthresh to zero after reading the newer NetCDF
% subgrid tables because the threshold is already included in uv_zmin.
SubgridTables.huthresh = 0.0;
SubgridTables.zvolmin = zvolmin;
SubgridTables.max_gradient = max_gradient;
SubgridTables.q_table_option = q_table_option;
SubgridTables.weight_option = weight_option;
SubgridTables.cellsize = cellsize;
SubgridTables.coarse_res = coarse_res;
SubgridTables.cell_area = coarse_res^2;
SubgridTables.g = 9.81;

% Compatibility metadata for older helpers; exact SFINCS mode does not use
% a fixed vertical-depth axis.
SubgridTables.depth_axis = 0:(nr_levels-1);
SubgridTables.dz = 1;
SubgridTables.maxDepth = nr_levels - 1;
SubgridTables.is_uniform = false;

SubgridTables.z_zmin = NaN(nrows_coarse, ncols_coarse);
SubgridTables.z_zmax = NaN(nrows_coarse, ncols_coarse);
SubgridTables.z_volmax = NaN(nrows_coarse, ncols_coarse);
SubgridTables.z_dep = NaN(nrows_coarse, ncols_coarse, nr_levels);
SubgridTables.z_level = NaN(nrows_coarse, ncols_coarse, nr_levels);
SubgridTables.z_volume = NaN(nrows_coarse, ncols_coarse, nr_levels);
SubgridTables.z_area = NaN(nrows_coarse, ncols_coarse, nr_levels);

nxface = max(ncols_coarse - 1, 1);
nyface = max(nrows_coarse - 1, 1);
SubgridTables.u_zmin = NaN(nrows_coarse, nxface);
SubgridTables.u_zmax = NaN(nrows_coarse, nxface);
SubgridTables.u_havg = NaN(nrows_coarse, nxface, nr_levels);
SubgridTables.u_nrep = NaN(nrows_coarse, nxface, nr_levels);
SubgridTables.u_pwet = NaN(nrows_coarse, nxface, nr_levels);
SubgridTables.u_navg_w = NaN(nrows_coarse, nxface);
SubgridTables.u_navg = NaN(nrows_coarse, nxface);
SubgridTables.u_ffit = NaN(nrows_coarse, nxface);

SubgridTables.v_zmin = NaN(nyface, ncols_coarse);
SubgridTables.v_zmax = NaN(nyface, ncols_coarse);
SubgridTables.v_havg = NaN(nyface, ncols_coarse, nr_levels);
SubgridTables.v_nrep = NaN(nyface, ncols_coarse, nr_levels);
SubgridTables.v_pwet = NaN(nyface, ncols_coarse, nr_levels);
SubgridTables.v_navg_w = NaN(nyface, ncols_coarse);
SubgridTables.v_navg = NaN(nyface, ncols_coarse);
SubgridTables.v_ffit = NaN(nyface, ncols_coarse);

for side_name = ["north", "south", "west", "east"]
    s = char(side_name);
    SubgridTables.([s '_zmin']) = NaN(nrows_coarse, ncols_coarse);
    SubgridTables.([s '_zmax']) = NaN(nrows_coarse, ncols_coarse);
    SubgridTables.([s '_havg']) = NaN(nrows_coarse, ncols_coarse, nr_levels);
    SubgridTables.([s '_nrep']) = NaN(nrows_coarse, ncols_coarse, nr_levels);
    SubgridTables.([s '_pwet']) = NaN(nrows_coarse, ncols_coarse, nr_levels);
    SubgridTables.([s '_navg_w']) = NaN(nrows_coarse, ncols_coarse);
    SubgridTables.([s '_navg']) = NaN(nrows_coarse, ncols_coarse);
    SubgridTables.([s '_ffit']) = NaN(nrows_coarse, ncols_coarse);
end

row_idx_start = floor((0:nrows_coarse-1) * coarse_res / cellsize) + 1;
row_idx_end = floor((1:nrows_coarse) * coarse_res / cellsize);
col_idx_start = floor((0:ncols_coarse-1) * coarse_res / cellsize) + 1;
col_idx_end = floor((1:ncols_coarse) * coarse_res / cellsize);

cellPatch = cell(nrows_coarse, ncols_coarse);
roughPatch = cell(nrows_coarse, ncols_coarse);

invert_el = NaN(nrows_coarse, ncols_coarse);

total_cells = nrows_coarse * ncols_coarse;
cell_count = 0;
for rowc = 1:nrows_coarse
    for colc = 1:ncols_coarse
        row_idx = row_idx_start(rowc):min(nrows, row_idx_end(rowc));
        col_idx = col_idx_start(colc):min(ncols, col_idx_end(colc));
        sub_DEM = DEM(row_idx, col_idx);
        if roughness_is_fine
            sub_ROUGH = ROUGH_INPUT(row_idx, col_idx);
        else
            sub_ROUGH = ROUGH_INPUT(rowc, colc) .* ones(size(sub_DEM));
        end

        cellPatch{rowc, colc} = sub_DEM;
        roughPatch{rowc, colc} = sub_ROUGH;

        [zmin, zmax, volmax, z_level, z_volume, z_area] = ...
            sfincs_cell_table(sub_DEM, cellsize, nr_levels, zvolmin, max_gradient);

        SubgridTables.z_zmin(rowc, colc) = zmin;
        SubgridTables.z_zmax(rowc, colc) = zmax;
        SubgridTables.z_volmax(rowc, colc) = volmax;
        SubgridTables.z_dep(rowc, colc, :) = reshape(z_level, 1, 1, nr_levels);
        SubgridTables.z_level(rowc, colc, :) = reshape(z_level, 1, 1, nr_levels);
        SubgridTables.z_volume(rowc, colc, :) = reshape(z_volume, 1, 1, nr_levels);
        SubgridTables.z_area(rowc, colc, :) = reshape(z_area, 1, 1, nr_levels);
        invert_el(rowc, colc) = zmin;

        cell_count = cell_count + 1;
        if verbose && (mod(cell_count, 100) == 0 || cell_count == total_cells)
            fprintf('SFINCS subgrid cell tables: %6.2f%% (%d of %d)\n', ...
                100 * cell_count / total_cells, cell_count, total_cells);
        end
    end
end

SubgridTables.invert_el = invert_el;

if ncols_coarse > 1
    total_u = nrows_coarse * (ncols_coarse - 1);
    face_count = 0;
    for rowc = 1:nrows_coarse
        for colc = 1:(ncols_coarse-1)
            [tbl] = sfincs_velocity_table( ...
                cellPatch{rowc, colc}, cellPatch{rowc, colc+1}, ...
                roughPatch{rowc, colc}, roughPatch{rowc, colc+1}, ...
                nr_levels, huthresh, q_table_option, weight_option, 'u');
            SubgridTables = assign_velocity_table(SubgridTables, 'u', rowc, colc, tbl);
            face_count = face_count + 1;
            if verbose && (mod(face_count, 100) == 0 || face_count == total_u)
                fprintf('SFINCS subgrid u tables:    %6.2f%% (%d of %d)\n', ...
                    100 * face_count / total_u, face_count, total_u);
            end
        end
    end
end

if nrows_coarse > 1
    total_v = (nrows_coarse - 1) * ncols_coarse;
    face_count = 0;
    for rowc = 1:(nrows_coarse-1)
        for colc = 1:ncols_coarse
            [tbl] = sfincs_velocity_table( ...
                cellPatch{rowc, colc}, cellPatch{rowc+1, colc}, ...
                roughPatch{rowc, colc}, roughPatch{rowc+1, colc}, ...
                nr_levels, huthresh, q_table_option, weight_option, 'v');
            SubgridTables = assign_velocity_table(SubgridTables, 'v', rowc, colc, tbl);
            face_count = face_count + 1;
            if verbose && (mod(face_count, 100) == 0 || face_count == total_v)
                fprintf('SFINCS subgrid v tables:    %6.2f%% (%d of %d)\n', ...
                    100 * face_count / total_v, face_count, total_v);
            end
        end
    end
end

for rowc = 1:nrows_coarse
    for colc = 1:ncols_coarse
        patch = cellPatch{rowc, colc};
        rough = roughPatch{rowc, colc};
        if isempty(patch) || all(~isfinite(patch(:)))
            continue;
        end
        edge_specs = {'north'; 'south'; 'west'; 'east'};
        for iside = 1:size(edge_specs, 1)
            side = edge_specs{iside, 1};
            tbl = sfincs_velocity_table_single(patch, rough, nr_levels, ...
                huthresh, q_table_option, weight_option, side);
            SubgridTables = assign_boundary_table(SubgridTables, side, rowc, colc, tbl);
        end
    end
end

SubgridTables = add_legacy_aliases(SubgridTables);
end

function [zmin, zmax, volmax, z_level, z_volume, z_area] = sfincs_cell_table( ...
    z_patch, cellsize, nr_levels, zvolmin, max_gradient)
% Port of HydroMT-SFINCS workflows.subgrid.subgrid_v_table for square pixels.
z = z_patch(:);
z = z(isfinite(z));
if isempty(z)
    zmin = NaN;
    zmax = NaN;
    volmax = NaN;
    z_level = NaN(nr_levels, 1);
    z_volume = NaN(nr_levels, 1);
    z_area = NaN(nr_levels, 1);
    return;
end

fine_area = cellsize^2;
cell_area_fine = numel(z) * fine_area;
ele_sort = sort(max(z, zvolmin));
for j = 2:numel(ele_sort)
    if ele_sort(j) <= ele_sort(j - 1)
        ele_sort(j) = ele_sort(j - 1) + 1.0e-6;
    end
end

depth = ele_sort - ele_sort(1);
volume = zeros(size(depth));
if numel(depth) > 1
    volume(2:end) = cumsum(diff(depth) .* fine_area .* (1:(numel(depth)-1)).');
end

z_volume = linspace(0, volume(end), nr_levels).';
z_level = interp1(volume, ele_sort, z_volume, 'linear', 'extrap');
dvol = volume(end) / max(nr_levels - 1, 1);
dzdh = sfincs_get_dzdh(z_level, z_volume, cell_area_fine);
niter = 0;
while max(dzdh) > max_gradient && abs(max(dzdh) - max_gradient) > eps(max_gradient) && niter < nr_levels
    idx = find(dzdh == max(dzdh));
    idx = idx(idx < nr_levels);
    z_level(idx + 1) = z_level(idx) + max_gradient * (dvol / cell_area_fine);
    dzdh = sfincs_get_dzdh(z_level, z_volume, cell_area_fine);
    niter = niter + 1;
end

zmin = min(z);
zmax = max(z_level);
volmax = z_volume(end);
z_area = zeros(nr_levels, 1);
for k = 1:nr_levels
    z_area(k) = sum(z_level(k) > z) * fine_area;
end
end

function dzdh = sfincs_get_dzdh(z, V, cell_area)
dzdh = zeros(size(z));
if numel(z) < 2
    return;
end
dV = max(diff(V), eps);
dzdh(1:end-1) = diff(z) ./ (dV ./ cell_area);
dzdh(end) = dzdh(end-1);
end

function tbl = sfincs_velocity_table(zA, zB, nA, nB, nr_levels, ...
    huthresh, q_table_option, weight_option, direction)
[zA, nA, zB, nB] = velocity_point_support(zA, nA, zB, nB, direction);
tbl = sfincs_velocity_table_from_sides(zA(:), zB(:), nA(:), nB(:), ...
    nr_levels, huthresh, q_table_option, weight_option);
end

function tbl = sfincs_velocity_table_single(z, n, nr_levels, ...
    huthresh, q_table_option, weight_option, side)
[z, n] = boundary_velocity_support(z, n, side);
tbl = sfincs_velocity_table_from_sides(z(:), [], n(:), [], ...
    nr_levels, huthresh, q_table_option, weight_option);
end

function [zA, nA, zB, nB] = velocity_point_support(zA, nA, zB, nB, direction)
switch lower(char(direction))
    case 'u'
        halfA = max(round(size(zA, 2) / 2), 1);
        halfB = max(round(size(zB, 2) / 2), 1);
        zA = zA(:, (end-halfA+1):end);
        nA = nA(:, (end-halfA+1):end);
        zB = zB(:, 1:halfB);
        nB = nB(:, 1:halfB);
    case 'v'
        halfA = max(round(size(zA, 1) / 2), 1);
        halfB = max(round(size(zB, 1) / 2), 1);
        zA = zA((end-halfA+1):end, :);
        nA = nA((end-halfA+1):end, :);
        zB = zB(1:halfB, :);
        nB = nB(1:halfB, :);
    otherwise
        error('Subgrid_Properties_Lookup:badVelocityDirection', ...
            'Velocity direction must be u or v.');
end
end

function [z, n] = boundary_velocity_support(z, n, side)
switch lower(char(side))
    case 'north'
        half = max(round(size(z, 1) / 2), 1);
        z = z(1:half, :);
        n = n(1:half, :);
    case 'south'
        half = max(round(size(z, 1) / 2), 1);
        z = z((end-half+1):end, :);
        n = n((end-half+1):end, :);
    case 'west'
        half = max(round(size(z, 2) / 2), 1);
        z = z(:, 1:half);
        n = n(:, 1:half);
    case 'east'
        half = max(round(size(z, 2) / 2), 1);
        z = z(:, (end-half+1):end);
        n = n(:, (end-half+1):end);
    otherwise
        error('Subgrid_Properties_Lookup:badBoundarySide', ...
            'Boundary side must be north, south, west, or east.');
end
end

function tbl = sfincs_velocity_table_from_sides(zA, zB, nA, nB, ...
    nr_levels, huthresh, q_table_option, weight_option)
% Port of HydroMT-SFINCS workflows.subgrid.subgrid_q_table, option 2 by default.
tbl = empty_velocity_table(nr_levels);
validA = isfinite(zA) & isfinite(nA) & nA > 0;
zA = zA(validA);
nA = nA(validA);
if isempty(zB)
    zB = zA;
    nB = nA;
else
    validB = isfinite(zB) & isfinite(nB) & nB > 0;
    zB = zB(validB);
    nB = nB(validB);
end
if isempty(zA) || isempty(zB)
    return;
end

zminA = min(zA);
zmaxA = max(zA);
zminB = min(zB);
zmaxB = max(zB);
zmin = max(zminA, zminB) + huthresh;
zmax = max(max(zmaxA, zmaxB), zmin + 0.01);
dlevel = (zmax - zmin) / (nr_levels - 1);

havg = zeros(nr_levels, 1);
nrep = zeros(nr_levels, 1);
pwet = zeros(nr_levels, 1);
z_all = [zA; zB];
n_all = [nA; nB];

for ibin = 1:nr_levels
    zbin = zmin + (ibin - 1) * dlevel;
    h = max(zbin - z_all, 0);
    hA = max(zbin - zA, 0);
    hB = max(zbin - zB, 0);

    qA = mean(hA.^(5/3) ./ nA, 'omitnan');
    hAmean = mean(hA, 'omitnan');
    qB = mean(hB.^(5/3) ./ nB, 'omitnan');
    hBmean = mean(hB, 'omitnan');
    qAll = mean(h.^(5/3) ./ n_all, 'omitnan');
    hAll = mean(h, 'omitnan');
    qMin = min(qA, qB);
    hMin = min(hAmean, hBmean);

    if q_table_option == 1
        w = (ibin - 1) / (nr_levels - 1);
        q = (1 - w) * qMin + w * qAll;
        hmean = hAll;
        pwet(ibin) = sum(zbin > z_all + huthresh) / numel(z_all);
    else
        pwetA = sum(zbin > zA) / numel(zA);
        pwetB = sum(zbin > zB) / numel(zB);
        if ibin == 1
            if pwetA < pwetB
                pwetA = 0;
            else
                pwetB = 0;
            end
        elseif ibin == nr_levels
            pwetA = 1;
            pwetB = 1;
        end

        if strcmp(weight_option, 'mean')
            w = 2 * min(pwetA, pwetB) / max(pwetA + pwetB, 1.0e-9);
            q = (1 - w) * qMin + w * qAll;
            hmean = (1 - w) * hMin + w * hAll;
        else
            if qA < qB
                q = qA;
                hmean = hAmean;
            else
                q = qB;
                hmean = hBmean;
            end
        end
        pwet(ibin) = 0.5 * (pwetA + pwetB);
    end

    havg(ibin) = hmean;
    if q > 0 && hmean > 0 && isfinite(q)
        nrep(ibin) = hmean^(5/3) / q;
    else
        nrep(ibin) = max(mean(n_all, 'omitnan'), 1e-4);
    end
end

nrep_top = nrep(end);
havg_top = havg(end);
zfit = zmax + zmax - zmin;
hfit = havg_top + zmax - zmin;
if strcmp(weight_option, 'mean')
    hfit_all = max(zfit - z_all, 0);
    qfit = mean(hfit_all.^(5/3) ./ n_all, 'omitnan');
    navg = mean(n_all, 'omitnan');
else
    hfitA = max(zfit - zA, 0);
    hfitB = max(zfit - zB, 0);
    qfitA = mean(hfitA.^(5/3) ./ nA, 'omitnan');
    qfitB = mean(hfitB.^(5/3) ./ nB, 'omitnan');
    if qfitA < qfitB
        qfit = qfitA;
        navg = mean(nA, 'omitnan');
    else
        qfit = qfitB;
        navg = mean(nB, 'omitnan');
    end
end
if qfit > 0 && isfinite(qfit)
    nfit = hfit^(5/3) / qfit;
else
    nfit = nrep_top;
end
gnavg_w = 9.81 * max(navg, 1e-4)^2;
gnavg_top = 9.81 * max(nrep_top, 1e-4)^2;
if gnavg_w / max(gnavg_top, eps) > 0.99 && gnavg_w / max(gnavg_top, eps) < 1.01
    ffit = 0;
else
    if navg > nrep_top
        if nfit > navg
            nfit = nrep_top + 0.9 * (navg - nrep_top);
        end
        if nfit < nrep_top
            nfit = nrep_top + 0.1 * (navg - nrep_top);
        end
    else
        if nfit < navg
            nfit = nrep_top + 0.9 * (navg - nrep_top);
        end
        if nfit > nrep_top
            nfit = nrep_top + 0.1 * (navg - nrep_top);
        end
    end
    gnfit = 9.81 * max(nfit, 1e-4)^2;
    zfit = max(zfit, zmax + 1.0e-6);
    gnavg_w = max(gnavg_w, gnfit + 1.0e-8);
    ffit = (((gnavg_w - gnavg_top) / (gnavg_w - gnfit)) - 1) / (zfit - zmax);
end
if ~isfinite(ffit)
    ffit = 0;
end

tbl.zmin = zmin;
tbl.zmax = zmax;
tbl.havg = havg;
tbl.nrep = 9.81 .* max(nrep, 1e-4).^2;
tbl.pwet = pwet;
tbl.navg = gnavg_w;
tbl.ffit = ffit;
end

function tbl = empty_velocity_table(nr_levels)
tbl.zmin = NaN;
tbl.zmax = NaN;
tbl.havg = NaN(nr_levels, 1);
tbl.nrep = NaN(nr_levels, 1);
tbl.pwet = NaN(nr_levels, 1);
tbl.navg = NaN;
tbl.ffit = NaN;
end

function S = assign_velocity_table(S, prefix, rowc, colc, tbl)
S.([prefix '_zmin'])(rowc, colc) = tbl.zmin;
S.([prefix '_zmax'])(rowc, colc) = tbl.zmax;
S.([prefix '_havg'])(rowc, colc, :) = reshape(tbl.havg, 1, 1, []);
S.([prefix '_nrep'])(rowc, colc, :) = reshape(tbl.nrep, 1, 1, []);
S.([prefix '_pwet'])(rowc, colc, :) = reshape(tbl.pwet, 1, 1, []);
S.([prefix '_navg_w'])(rowc, colc) = tbl.navg;
S.([prefix '_navg'])(rowc, colc) = tbl.navg;
S.([prefix '_ffit'])(rowc, colc) = tbl.ffit;
end

function S = assign_boundary_table(S, side, rowc, colc, tbl)
S.([side '_zmin'])(rowc, colc) = tbl.zmin;
S.([side '_zmax'])(rowc, colc) = tbl.zmax;
S.([side '_havg'])(rowc, colc, :) = reshape(tbl.havg, 1, 1, []);
S.([side '_nrep'])(rowc, colc, :) = reshape(tbl.nrep, 1, 1, []);
S.([side '_pwet'])(rowc, colc, :) = reshape(tbl.pwet, 1, 1, []);
S.([side '_navg_w'])(rowc, colc) = tbl.navg;
S.([side '_navg'])(rowc, colc) = tbl.navg;
S.([side '_ffit'])(rowc, colc) = tbl.ffit;
end

function S = add_legacy_aliases(S)
S.eta_cell = S.z_level;
S.volume_cell = S.z_volume;
S.area_cell = S.z_area;
S.eta_top_cell = S.z_zmax;
S.Vmax_cell = S.z_volmax;

S.invert_x = S.u_zmin;
S.eta_top_x = S.u_zmax;
S.hrep_x = S.u_havg;
S.phi_x = S.u_pwet;
S.wetfrac_x = S.u_pwet;
S.n_x = S.u_nrep;
S.nrep_x = S.u_nrep;
S.area_x = S.u_havg .* S.coarse_res;
S.width_x = S.u_pwet .* S.coarse_res;
S.K_x = sqrt(S.g ./ S.u_nrep) .* max(S.u_havg, 0).^(5/3) .* S.coarse_res;
S.perimeter_x = NaN(size(S.u_havg));
S.Rh_x = NaN(size(S.u_havg));

S.invert_y = S.v_zmin;
S.eta_top_y = S.v_zmax;
S.hrep_y = S.v_havg;
S.phi_y = S.v_pwet;
S.wetfrac_y = S.v_pwet;
S.n_y = S.v_nrep;
S.nrep_y = S.v_nrep;
S.area_y = S.v_havg .* S.coarse_res;
S.width_y = S.v_pwet .* S.coarse_res;
S.K_y = sqrt(S.g ./ S.v_nrep) .* max(S.v_havg, 0).^(5/3) .* S.coarse_res;
S.perimeter_y = NaN(size(S.v_havg));
S.Rh_y = NaN(size(S.v_havg));

for side_name = ["north", "south", "west", "east"]
    side = char(side_name);
    S.(['invert_' side]) = S.([side '_zmin']);
    S.(['eta_top_' side]) = S.([side '_zmax']);
    S.(['hrep_' side]) = S.([side '_havg']);
    S.(['phi_' side]) = S.([side '_pwet']);
    S.(['wetfrac_' side]) = S.([side '_pwet']);
    S.(['n_' side]) = S.([side '_nrep']);
    S.(['nrep_' side]) = S.([side '_nrep']);
    S.(['area_' side]) = S.([side '_havg']) .* S.coarse_res;
    S.(['width_' side]) = S.([side '_pwet']) .* S.coarse_res;
    S.(['K_' side]) = sqrt(S.g ./ S.([side '_nrep'])) .* ...
        max(S.([side '_havg']), 0).^(5/3) .* S.coarse_res;
    S.(['perimeter_' side]) = NaN(size(S.([side '_havg'])));
    S.(['Rh_' side]) = NaN(size(S.([side '_havg'])));
end
end
