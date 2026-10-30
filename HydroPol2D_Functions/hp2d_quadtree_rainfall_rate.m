function rate = hp2d_quadtree_rainfall_rate(time_s, files, forcing_root, interval_s, latitude, longitude, mapping, valid)
%HP2D_QUADTREE_RAINFALL_RATE Rain rate of each mesh cell from the prepared rain maps [m/s].
%   Map i (files sorted by name) is the rate from (i-1)*interval_s to
%   i*interval_s; after the last map the rate is 0, as on the grid and in
%   HydroPol2D-Python. Each map is sampled (nearest) at the centre of every
%   land-use raster cell (latitude, longitude, south-up like the overlap)
%   and averaged over each mesh cell with the overlap weights, over the
%   active raster cells (valid), as Python samples its HRU pieces.

persistent cached_path cached_rate
index = floor(time_s/interval_s) + 1;
if index > numel(files)
    rate = zeros(size(mapping.raster_to_mesh, 1), 1);
    return
end
index = max(1, index);
path = fullfile(forcing_root, 'Rainfall', files(index).name);
if isequal(cached_path, path)
    rate = cached_rate;
    return
end
[raw, reference] = readgeoraster(path, 'OutputType', 'double');
for missing = double(georasterinfo(path).MissingDataIndicator(:))'
    raw(raw == missing) = NaN;
end
[column, row] = geographicToIntrinsic(reference, latitude, longitude);
% The map cell that contains each point. interp2 'nearest' treats a point a
% rounding error outside [1, n] as off the map, which zeroed edge cells.
column = round(column(:)); row = round(row(:));
inside = column >= 1 & column <= size(raw, 2) & row >= 1 & row <= size(raw, 1);
sampled = zeros(numel(column), 1);
sampled(inside) = raw(sub2ind(size(raw), row(inside), column(inside)));
sampled(~isfinite(sampled) | ~valid(:)) = 0;
rate = (mapping.raster_to_mesh*max(sampled(:), 0)) ./ ...
    max(mapping.raster_to_mesh*double(valid(:)), eps) ./ 1000 ./ 3600;
cached_path = path; cached_rate = rate;
end
