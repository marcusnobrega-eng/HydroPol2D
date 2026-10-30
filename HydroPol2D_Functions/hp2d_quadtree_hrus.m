function hrus = hp2d_quadtree_hrus(raster_to_mesh, lulc, soil, continuous)
%HP2D_QUADTREE_HRUS Land-use x soil pieces of each mesh cell, with their shares.
%   raster_to_mesh is the overlap operator (n_cells x n_raster, weights
%   overlap area / mesh cell area). lulc and soil are the class codes of the
%   raster cells (south-up vectors); NaN marks cells outside the domain.
%   continuous is a struct of raster vectors (LAI, initial theta, impervious
%   fraction, ...) to average over each piece.
%
%   HydroPol2D-Python gives every mesh-cell x raster-cell overlap its own HRU
%   (voronoi.py, hru_host / hru_raster). The class parameters are the same in
%   every raster cell of one land-use x soil pair, so here the overlaps of a
%   mesh cell are grouped by that pair: one HRU per pair, its share of the
%   cell's active area, and the area mean of each continuous field over the
%   group. A cell made of one class gets one HRU; a mixed urban cell keeps its
%   impervious and pervious parts apart instead of taking the majority class.
%
%   hrus.fraction [n x h]: share of the cell's active area (rows sum to 1).
%   hrus.lulc, hrus.soil [n x h]: class codes of each HRU.
%   hrus.<field> [n x h]: area means of continuous.<field> (NaN where the
%       group has no finite value). Unused HRU slots have fraction 0 and copy
%       the cell's first HRU, so every value stays finite.

if nargin < 4, continuous = struct(); end
lulc = double(lulc(:)); soil = double(soil(:));
n = size(raster_to_mesh, 1);
n_raster = size(raster_to_mesh, 2);
valid = isfinite(lulc) & isfinite(soil);
[pairs, ~, pair_of_valid] = unique([lulc(valid) soil(valid)], 'rows');
valid_index = find(valid);
indicator = sparse(valid_index, pair_of_valid, 1, n_raster, size(pairs, 1));
share = raster_to_mesh * indicator;                       % n x pairs, area shares

[cell_id, pair_id, weight] = find(share);
keep = weight > 0;
cell_id = cell_id(keep); pair_id = pair_id(keep); weight = weight(keep);
[~, order] = sortrows([cell_id, -weight, pair_id]);
cell_id = cell_id(order); pair_id = pair_id(order); weight = weight(order);
first = [true; diff(cell_id) ~= 0];
start = find(first);
slot = (1:numel(cell_id))' - repelem(start, diff([start; numel(cell_id) + 1]));
slot = slot + 1;
h = max([slot; 1]);
covered = accumarray(cell_id, weight, [n 1]);
assert(all(covered > 0), 'HydroPol2D:IncompletePreparedRasterCoverage', ...
    'Every mesh cell must overlap at least one land-use and soil raster cell.');

linear = sub2ind([n h], cell_id, slot);
hrus.fraction = zeros(n, h);
hrus.fraction(linear) = weight ./ covered(cell_id);
hrus.lulc = fill_slots(n, h, linear, pairs(pair_id, 1), cell_id, slot);
hrus.soil = fill_slots(n, h, linear, pairs(pair_id, 2), cell_id, slot);
hrus.pair_index = fill_slots(n, h, linear, pair_id, cell_id, slot);

names = fieldnames(continuous);
for k = 1:numel(names)
    values = double(continuous.(names{k})(:));
    finite = isfinite(values);
    values(~finite) = 0;
    totals = raster_to_mesh * (indicator .* values);
    weights = raster_to_mesh * (indicator .* double(finite));
    group_total = full(totals(sub2ind(size(totals), cell_id, pair_id)));
    group_weight = full(weights(sub2ind(size(weights), cell_id, pair_id)));
    mean_value = group_total ./ group_weight;
    mean_value(group_weight <= 0) = NaN;
    hrus.(names{k}) = fill_slots(n, h, linear, mean_value, cell_id, slot);
end
end

function value = fill_slots(n, h, linear, data, cell_id, slot)
% Used slots take the data; the others copy the cell's first HRU.
value = nan(n, h);
value(linear) = data;
first_value = nan(n, 1);
first_value(cell_id(slot == 1)) = data(slot == 1);
unused = isnan(value) & ~ismember(reshape(1:n*h, n, h), linear);
fill = repmat(first_value, 1, h);
value(unused) = fill(unused);
end
