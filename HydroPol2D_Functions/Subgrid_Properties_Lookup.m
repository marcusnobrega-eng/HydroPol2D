function [SubgridTables, invert_el] = Subgrid_Properties_Lookup( ...
    DEM_raster, Roughness_raster, Reference_raster, coarse_res, varargin)
%SUBGRID_PROPERTIES_LOOKUP Build SFINCS-style lookup-subgrid tables.
%
% This routine follows the subgrid corrections for the linear inertial
% equations described by van Ormondt et al. (2025). The production solver
% uses the SFINCS fields below:
%
%   Cell continuity points:
%       z_zmin, z_zmax, z_volmax, z_level
%
%   Velocity points:
%       u_zmin, u_zmax, u_havg, u_nrep, u_pwet, u_navg, u_ffit
%       v_zmin, v_zmax, v_havg, v_nrep, v_pwet, v_navg, v_ffit
%
% Legacy HydroPol fields are also populated as aliases for diagnostics and
% older plotting routines, but lookup-subgrid routing should query the
% SFINCS fields directly.

p = inputParser;
p.addParameter('nr_levels', 20, @(x)isscalar(x) && x >= 2);
p.addParameter('nlevels', [], @(x)isempty(x) || (isscalar(x) && x >= 2));
p.addParameter('tol_rel', 5e-3, @(x)isscalar(x) && x > 0);
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
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

SubgridTables = struct();
SubgridTables.sfincs_exact = true;
SubgridTables.nr_levels = nr_levels;
SubgridTables.cellsize = cellsize;
SubgridTables.coarse_res = coarse_res;
SubgridTables.cell_area = coarse_res^2;

% Compatibility metadata for older helpers; exact SFINCS mode does not use
% a fixed vertical-depth axis.
SubgridTables.depth_axis = 0:(nr_levels-1);
SubgridTables.dz = 1;
SubgridTables.maxDepth = nr_levels - 1;
SubgridTables.is_uniform = false;

SubgridTables.z_zmin = NaN(nrows_coarse, ncols_coarse);
SubgridTables.z_zmax = NaN(nrows_coarse, ncols_coarse);
SubgridTables.z_volmax = NaN(nrows_coarse, ncols_coarse);
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
SubgridTables.u_navg = NaN(nrows_coarse, nxface);
SubgridTables.u_ffit = NaN(nrows_coarse, nxface);

SubgridTables.v_zmin = NaN(nyface, ncols_coarse);
SubgridTables.v_zmax = NaN(nyface, ncols_coarse);
SubgridTables.v_havg = NaN(nyface, ncols_coarse, nr_levels);
SubgridTables.v_nrep = NaN(nyface, ncols_coarse, nr_levels);
SubgridTables.v_pwet = NaN(nyface, ncols_coarse, nr_levels);
SubgridTables.v_navg = NaN(nyface, ncols_coarse);
SubgridTables.v_ffit = NaN(nyface, ncols_coarse);

for side_name = ["north", "south", "west", "east"]
    s = char(side_name);
    SubgridTables.([s '_zmin']) = NaN(nrows_coarse, ncols_coarse);
    SubgridTables.([s '_zmax']) = NaN(nrows_coarse, ncols_coarse);
    SubgridTables.([s '_havg']) = NaN(nrows_coarse, ncols_coarse, nr_levels);
    SubgridTables.([s '_nrep']) = NaN(nrows_coarse, ncols_coarse, nr_levels);
    SubgridTables.([s '_pwet']) = NaN(nrows_coarse, ncols_coarse, nr_levels);
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
            sfincs_cell_table(sub_DEM, cellsize, nr_levels);

        SubgridTables.z_zmin(rowc, colc) = zmin;
        SubgridTables.z_zmax(rowc, colc) = zmax;
        SubgridTables.z_volmax(rowc, colc) = volmax;
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
                nr_levels, 'u');
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
                nr_levels, 'v');
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
            tbl = sfincs_velocity_table_single(patch, rough, nr_levels, side);
            SubgridTables = assign_boundary_table(SubgridTables, side, rowc, colc, tbl);
        end
    end
end

SubgridTables = add_legacy_aliases(SubgridTables);
end

function [zmin, zmax, volmax, z_level, z_volume, z_area] = sfincs_cell_table(z_patch, cellsize, nr_levels)
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
zmin = min(z);
zmax = max(z);
volmax = sum(max(zmax - z, 0)) * fine_area;
Vlevels = linspace(0, volmax, nr_levels).';
z_level = zeros(nr_levels, 1);
z_area = zeros(nr_levels, 1);

for k = 1:nr_levels
    z_level(k) = invert_volume_exact(z, Vlevels(k), zmin, zmax, volmax, fine_area);
    z_area(k) = sum(z_level(k) > z) * fine_area;
end
z_volume = Vlevels;
end

function zs = invert_volume_exact(z, V, zmin, zmax, volmax, fine_area)
if ~isfinite(V) || V <= 0
    zs = zmin;
    return;
end
if volmax <= eps || zmax <= zmin
    zs = zmax + V / (numel(z) * fine_area);
    return;
end
if V >= volmax
    zs = zmax + (V - volmax) / (numel(z) * fine_area);
    return;
end
lo = zmin;
hi = zmax;
for iter = 1:60
    mid = 0.5 * (lo + hi);
    Vm = sum(max(mid - z, 0)) * fine_area;
    if Vm < V
        lo = mid;
    else
        hi = mid;
    end
end
zs = 0.5 * (lo + hi);
end

function tbl = sfincs_velocity_table(zA, zB, nA, nB, nr_levels, direction)
[zA, nA, zB, nB] = velocity_point_support(zA, nA, zB, nB, direction);
zA = zA(:);
zB = zB(:);
nA = nA(:);
nB = nB(:);
validA = isfinite(zA) & isfinite(nA) & nA > 0;
validB = isfinite(zB) & isfinite(nB) & nB > 0;
zA = zA(validA);
zB = zB(validB);
nA = nA(validA);
nB = nB(validB);

if isempty(zA) || isempty(zB)
    tbl = empty_velocity_table(nr_levels);
    return;
end

zmin = max(min(zA), min(zB));
zmax = max(max(zA), max(zB));
z = [zA; zB];
n = [nA; nB];
tbl = sfincs_velocity_table_from_vectors(z, n, zmin, zmax, nr_levels);
end

function tbl = sfincs_velocity_table_single(z, n, nr_levels, side)
[z, n] = boundary_velocity_support(z, n, side);
z = z(:);
n = n(:);
valid = isfinite(z) & isfinite(n) & n > 0;
z = z(valid);
n = n(valid);
if isempty(z)
    tbl = empty_velocity_table(nr_levels);
    return;
end
zmin = min(z);
zmax = max(z);
tbl = sfincs_velocity_table_from_vectors(z, n, zmin, zmax, nr_levels);
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

function tbl = sfincs_velocity_table_from_vectors(z, n, zmin, zmax, nr_levels)
tbl = empty_velocity_table(nr_levels);
if isempty(z) || isempty(n) || ~isfinite(zmin) || ~isfinite(zmax)
    return;
end

if zmax > zmin
    levels = linspace(zmin, zmax, nr_levels).';
else
    levels = zmin .* ones(nr_levels, 1);
end

havg = zeros(nr_levels, 1);
nrep = Inf(nr_levels, 1);
pwet = zeros(nr_levels, 1);

for k = 1:nr_levels
    [havg(k), nrep(k), pwet(k)] = sfincs_velocity_values(levels(k), z, n, zmin);
end

navg = mean(n, 'omitnan');
nM = nrep(end);
HG_M = havg(end);
dzfit = max(zmax - zmin, eps);
zfit = zmax + dzfit;
denom_fit = mean((max(zfit - max(z, zmin), 0).^(5/3)) ./ n, 'omitnan');
if denom_fit > 0 && isfinite(denom_fit)
    nfit = (HG_M + zfit - zmax)^(5/3) / denom_fit;
else
    nfit = nM;
end
if isfinite(navg) && isfinite(nM) && isfinite(nfit) && abs(navg - nfit) > eps
    beta = ((navg - nM) / (navg - nfit) - 1) / dzfit;
else
    beta = 0;
end
if ~isfinite(beta)
    beta = 0;
end

tbl.zmin = zmin;
tbl.zmax = zmax;
tbl.havg = havg;
tbl.nrep = nrep;
tbl.pwet = pwet;
tbl.navg = navg;
tbl.ffit = beta;
end

function [havg, nrep, pwet] = sfincs_velocity_values(zu, z, n, zmin)
h = max(zu - z, 0);
pwet = mean(h > 0);
havg = mean(h);
denom = mean((max(zu - max(z, zmin), 0).^(5/3)) ./ n, 'omitnan');
if havg > 0 && denom > 0 && isfinite(denom)
    nrep = havg^(5/3) / denom;
else
    nrep = Inf;
end
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
S.([prefix '_navg'])(rowc, colc) = tbl.navg;
S.([prefix '_ffit'])(rowc, colc) = tbl.ffit;
end

function S = assign_boundary_table(S, side, rowc, colc, tbl)
S.([side '_zmin'])(rowc, colc) = tbl.zmin;
S.([side '_zmax'])(rowc, colc) = tbl.zmax;
S.([side '_havg'])(rowc, colc, :) = reshape(tbl.havg, 1, 1, []);
S.([side '_nrep'])(rowc, colc, :) = reshape(tbl.nrep, 1, 1, []);
S.([side '_pwet'])(rowc, colc, :) = reshape(tbl.pwet, 1, 1, []);
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
S.K_x = (1 ./ S.u_nrep) .* max(S.u_havg, 0).^(5/3) .* S.coarse_res;
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
S.K_y = (1 ./ S.v_nrep) .* max(S.v_havg, 0).^(5/3) .* S.coarse_res;
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
    S.(['K_' side]) = (1 ./ S.([side '_nrep'])) .* ...
        max(S.([side '_havg']), 0).^(5/3) .* S.coarse_res;
    S.(['perimeter_' side]) = NaN(size(S.([side '_havg'])));
    S.(['Rh_' side]) = NaN(size(S.([side '_havg'])));
end
end
