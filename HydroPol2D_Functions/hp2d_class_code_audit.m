function [unmapped_mask, fallback_values, Audit] = hp2d_class_code_audit( ...
    class_grid, valid_mask, class_index, parameter_values, parameter_names, class_label, cell_area_m2)
%HP2D_CLASS_CODE_AUDIT Validate class mappings and describe average fallbacks.
%
% Unmapped raster classes are assigned an area-weighted average of the
% properties represented by mapped cells. Invalid tables remain errors.

class_index = double(class_index(:));
parameter_values = double(parameter_values);
parameter_names = string(parameter_names(:))';

if size(parameter_values,1) ~= numel(class_index)
    error('%s class table has inconsistent index and parameter row counts.', class_label);
end
if size(parameter_values,2) ~= numel(parameter_names)
    error('%s class table has inconsistent parameter names and values.', class_label);
end
if isempty(class_index) || any(~isfinite(class_index)) || ...
        any(abs(class_index - round(class_index)) > 1e-9) || numel(unique(class_index)) ~= numel(class_index)
    error('%s class table must contain unique finite integer class indices.', class_label);
end
if ~isequal(size(class_grid), size(valid_mask))
    error('%s class raster and domain mask must have the same dimensions.', class_label);
end
if ~isscalar(cell_area_m2) || ~isfinite(cell_area_m2) || cell_area_m2 <= 0
    error('Cell area must be a positive scalar in m2.');
end

valid_mask = logical(valid_mask) & isfinite(class_grid);
grid_codes = unique(double(class_grid(valid_mask)));
if isempty(grid_codes)
    error('%s raster has no valid cells in the simulation domain.', class_label);
end
if any(abs(grid_codes - round(grid_codes)) > 1e-9)
    error('%s raster contains non-integer class codes.', class_label);
end

counts = zeros(numel(class_index), 1);
for i = 1:numel(class_index)
    counts(i) = nnz(valid_mask & class_grid == class_index(i));
end
if sum(counts) == 0
    error(['%s table has no code in common with the raster, so an average ', ...
        'fallback cannot be calculated.'], class_label);
end

fallback_values = nan(1, size(parameter_values,2));
for j = 1:size(parameter_values,2)
    values = parameter_values(:,j);
    use = counts > 0 & isfinite(values);
    if any(use)
        fallback_values(j) = sum(counts(use) .* values(use)) / sum(counts(use));
    end
end

mapped_mask = false(size(class_grid));
for i = 1:numel(class_index)
    mapped_mask = mapped_mask | (class_grid == class_index(i));
end
unmapped_mask = valid_mask & ~mapped_mask;

n_codes = numel(grid_codes);
map_name = repmat(string(class_label), n_codes, 1);
cell_count = zeros(n_codes, 1);
cell_area = zeros(n_codes, 1);
status = strings(n_codes, 1);
table_index = nan(n_codes, 1);
fallback_description = strings(n_codes, 1);
for i = 1:n_codes
    code = grid_codes(i);
    cell_count(i) = nnz(valid_mask & class_grid == code);
    cell_area(i) = cell_count(i) * cell_area_m2;
    table_row = find(class_index == code, 1);
    if isempty(table_row)
        status(i) = "fallback_area_weighted_mean";
        fallback_description(i) = format_fallback(parameter_names, fallback_values);
    else
        status(i) = "mapped";
        table_index(i) = class_index(table_row);
    end
end

Audit = table(map_name, grid_codes, cell_count, cell_area, status, table_index, fallback_description, ...
    'VariableNames', {'class_map','raster_code','cell_count','area_m2','mapping_status', ...
    'table_index','fallback_properties'});
end

function text_value = format_fallback(parameter_names, values)
parts = strings(1, numel(parameter_names));
for i = 1:numel(parameter_names)
    if isfinite(values(i))
        parts(i) = parameter_names(i) + "=" + string(sprintf('%.8g', values(i)));
    else
        parts(i) = parameter_names(i) + "=NaN";
    end
end
text_value = strjoin(parts, '; ');
end
