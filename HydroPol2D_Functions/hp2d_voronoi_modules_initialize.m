function modules = hp2d_voronoi_modules_initialize(mesh, mesh_file, config, overlap_file, uses_channels, uses_subgrid)
%HP2D_VORONOI_MODULES_INITIALIZE The grid's optional modules on mesh cells.
%   Reservoir controls, water quality (mass-based wash-off) and the human
%   risk index, as HydroPol2D-Python runs them on a mesh (mesh_modules/).
%   Each is read from the prepared case and is off when its field is absent
%   or empty:
%
%   config.reservoirs     the grid's Reservoir_Data table: k1, h1, k2, k3, h2,
%       k4 and the map coordinates x_us, y_us, x_ds1, y_ds1, x_ds2, y_ds2
%       (destination 0, 0 leaves the model; NaN switches the outlet off).
%       Each point goes to the mesh cell that holds most of the raster cell
%       under it, so the case needs an overlap file. Mesh cell ids can be
%       given instead as source_cell, ds1_cell, ds2_cell (0 = leaves, NaN =
%       off).
%   config.water_quality  C_3 [kg/h per (m3/s)^C_4 per kg], C_4, min_Bt
%       [g/m2], and initial_mass_kg per cell, or C_1 [kg/ha], C_2 [1/day]
%       and ADD [days] for B0 = C_1 (1 - exp(-C_2 ADD)), as on the grid.
%   config.human_risk     mu, Cd, weight_person [kg], width1_person and
%       width2_person [m], ro_water [kg/m3], gravity [m/s2] (the grid's
%       flag_human_instability = 1 index). Mode 3 does not run on meshes.

modules = struct('active', false, 'reservoirs', [], 'water_quality', [], 'human_risk', []);
n = mesh.n_cells;
mapping = [];

if isfield(config, 'reservoirs') && ~isempty(config.reservoirs)
    % As on the grid and in HydroPol2D-Python: a reservoir cell stores only
    % its surface water.
    assert(~uses_channels && ~uses_subgrid, 'HydroPol2D:MeshModuleNotSupported', ...
        ['Reservoir controls are not coupled to mesh river channels or sub-grid ' ...
         'tables. Use a mesh without them, or the regular (D4) grid.']);
    R = config.reservoirs;
    n_controls = numel(R.k1);
    if isfield(R, 'source_cell')
        source = double(R.source_cell(:));
        ds1 = double(R.ds1_cell(:));
        ds2 = double(R.ds2_cell(:));
    else
        mapping = read_mapping(mapping, overlap_file, mesh_file);
        source = zeros(n_controls, 1);
        ds1 = zeros(n_controls, 1);
        ds2 = zeros(n_controls, 1);
        for ii = 1:n_controls
            source(ii) = host_of_point(mapping, R.x_us(ii), R.y_us(ii), ...
                sprintf('Reservoir %d', ii));
            ds1(ii) = destination_cell(mapping, R.x_ds1(ii), R.y_ds1(ii), ii);
            ds2(ii) = destination_cell(mapping, R.x_ds2(ii), R.y_ds2(ii), ii);
        end
    end
    assert(all(source >= 1 & source <= n & source == round(source)), ...
        'HydroPol2D:InvalidVoronoiCase', 'Reservoir source cells must be mesh cell ids.');
    for ii = 1:n_controls
        assert(ds1(ii) ~= source(ii) && ds2(ii) ~= source(ii), 'HydroPol2D:InvalidVoronoiCase', ...
            ['Reservoir %d: its storage cell and release cell fall in the same mesh cell, ' ...
             'so the release would return to the reservoir. Refine the mesh there, or ' ...
             'move the release cell.'], ii);
    end
    modules.reservoirs = struct('source', source, 'ds1', ds1, 'ds2', ds2, ...
        'k1', double(R.k1(:)), 'h1', double(R.h1(:)), 'k2', double(R.k2(:)), ...
        'k3', double(R.k3(:)), 'h2', double(R.h2(:)), 'k4', double(R.k4(:)), ...
        'discharge_m3_s', zeros(n_controls, 1), 'cumulative_outflow_m3', 0);
end

if isfield(config, 'water_quality') && ~isempty(config.water_quality)
    assert(~uses_channels, 'HydroPol2D:MeshModuleNotSupported', ...
        ['Water quality on a mesh does not carry pollutant along river channels. ' ...
         'Use a mesh without them, or the regular (D4) grid.']);
    Q = config.water_quality;
    assert(~(isfield(Q, 'flag_wq_model') && Q.flag_wq_model == 1), ...
        'HydroPol2D:MeshModuleNotSupported', ...
        'The rating-curve wash-off (flag_wq_model = 1) runs on the grid only.');
    area = mesh.cell_area(:);
    if isfield(Q, 'initial_mass_kg')
        mass = expand(Q.initial_mass_kg, n);
    else
        mass = expand(Q.C_1, n) .* (1 - exp(-expand(Q.C_2, n) .* expand(Q.ADD, n))) .* area / 1e4;
    end
    modules.water_quality = struct( ...
        'coefficient_kg_s', expand(Q.C_3, n) / 3600, ...
        'exponent', expand(Q.C_4, n), ...
        'threshold_kg', expand(Q.min_Bt, n) / 1000 .* area, ...
        'mass_kg', max(mass, 0), 'cumulative_export_kg', 0, ...
        'outlet_concentration_mg_l', 0, 'cell_concentration_mg_l', zeros(n, 1));
end

if isfield(config, 'human_risk') && ~isempty(config.human_risk)
    H = config.human_risk;
    assert(~(isfield(H, 'mode') && H.mode ~= 1), 'HydroPol2D:MeshModuleNotSupported', ...
        ['Human risk mode 3 (slide, topple, drowning) runs on the grid only; meshes ' ...
         'run the risk index (mode 1).']);
    [discharge_x, discharge_y] = cell_discharge_operators(mesh);
    maximum_velocity = inf;
    if isfield(H, 'maximum_velocity_m_s'), maximum_velocity = H.maximum_velocity_m_s; end
    modules.human_risk = struct('parameters', H, 'discharge_x', discharge_x, ...
        'discharge_y', discharge_y, 'maximum_velocity_m_s', maximum_velocity, ...
        'risk', zeros(n, 1), 'max_risk', zeros(n, 1));
end

modules.active = ~isempty(modules.reservoirs) || ~isempty(modules.water_quality) || ...
    ~isempty(modules.human_risk);
end

function [discharge_x, discharge_y] = cell_discharge_operators(mesh)
% (cell x edge) operators from edge discharge per width to the cell's hu, hv:
% (hu, hv)_c = (1/A_c) sum_e F_e L_e (x_e - x_c), F_e leaving the cell
% (Perot 2000). Exact for a uniform flow; the mean of opposite faces on a
% square cell.
edges = (1:mesh.n_edges)';
internal = mesh.edge_neighbor > 0;
cells = [mesh.edge_owner(:); mesh.edge_neighbor(internal)];
columns = [edges; edges(internal)];
leaving = [ones(mesh.n_edges, 1); -ones(nnz(internal), 1)];
weight = leaving .* mesh.edge_length(columns) ./ mesh.cell_area(cells);
lever_x = mesh.edge_midpoint_x(columns) - mesh.cell_x(cells);
lever_y = mesh.edge_midpoint_y(columns) - mesh.cell_y(cells);
discharge_x = sparse(cells, columns, weight .* lever_x, mesh.n_cells, mesh.n_edges);
discharge_y = sparse(cells, columns, weight .* lever_y, mesh.n_cells, mesh.n_edges);
end

function cell = destination_cell(mapping, x, y, index)
if isnan(x) || isnan(y)
    cell = NaN;
elseif x == 0 && y == 0
    cell = 0;
else
    cell = host_of_point(mapping, x, y, sprintf('The release cell of reservoir %d', index));
end
end

function mapping = read_mapping(mapping, overlap_file, mesh_file)
if ~isempty(mapping), return; end
assert(~isempty(overlap_file) && isfile(overlap_file), 'HydroPol2D:InvalidVoronoiCase', ...
    ['Reservoir map coordinates need the case overlap_file to find their mesh cells; ' ...
     'or give source_cell, ds1_cell and ds2_cell.']);
mapping = HydroPol2D_Read_Overlap(overlap_file, mesh_file);
end

function host = host_of_point(mapping, x, y, label)
% The mesh cell that holds most of the raster cell under (x, y).
col = edge_bin(mapping.x_edges, x);
row = edge_bin(mapping.y_edges, y);
assert(~isnan(col) && ~isnan(row), 'HydroPol2D:InvalidVoronoiCase', ...
    '%s (x %.3f, y %.3f) lies outside the grid.', label, x, y);
if mapping.y_edges(1) > mapping.y_edges(end)
    row = mapping.raster_rows - row + 1;   % count from the south
end
raster_cell = (row - 1) * mapping.raster_cols + col;   % south-up, row-major
[weight, host] = max(mapping.mesh_to_raster(raster_cell, :));
assert(full(weight) > 0, 'HydroPol2D:InvalidVoronoiCase', ...
    ['%s (x %.3f, y %.3f) lies outside the mesh. Move it inside the mesh, or run ' ...
     'this model on the regular (D4) grid.'], label, x, y);
end

function index = edge_bin(edges, value)
edges = edges(:);
low = min(edges(1:end-1), edges(2:end));
high = max(edges(1:end-1), edges(2:end));
index = find(value >= low & value < high, 1);
if isempty(index), index = NaN; end
end

function value = expand(value, n)
value = double(value(:));
if isscalar(value), value = repmat(value, n, 1); end
assert(numel(value) == n, 'HydroPol2D:InvalidVoronoiCase', ...
    'A mesh module field must be a scalar or one value per mesh cell.');
end
