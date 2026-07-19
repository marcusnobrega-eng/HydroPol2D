function results = run_topotoolbox_lite_verification(model_root, source_root)
%RUN_TOPOTOOLBOX_LITE_VERIFICATION Compare full and curated terrain runtimes.
%   SOURCE_ROOT is a local TopoToolbox v2.4 snapshot used only for this
%   pre-release comparison. The released model never adds that folder.

arguments
    model_root (1,:) char
    source_root (1,:) char
end

lite_root = fullfile(model_root, 'third_party', 'topotoolbox_lite');
dem_path = fullfile(model_root, 'Static', 'DEM.tif');
release_root = fileparts(mfilename('fullpath'));

assert(isfolder(lite_root), 'Bundled runtime is missing: %s', lite_root);
assert(isfolder(source_root), 'Source snapshot is missing: %s', source_root);
assert(isfile(dem_path), 'Representative DEM is missing: %s', dem_path);

source = terrain_signature(model_root, source_root, dem_path);
lite = terrain_signature(model_root, lite_root, dem_path);

compare_exact(source.crop_z, lite.crop_z, 'Crop raster values');
compare_exact(source.resample_z, lite.resample_z, 'Resampled raster values');
compare_exact(source.filled_z, lite.filled_z, 'Filled DEM');
compare_exact(source.flowacc_z, lite.flowacc_z, 'D8 flow accumulation');
compare_exact(source.receiver_ix, lite.receiver_ix, 'D8 receiver indices');
compare_exact(source.receiver_ixc, lite.receiver_ixc, 'D8 donor indices');
compare_exact(source.stream_ixgrid, lite.stream_ixgrid, 'Stream mask indices');
compare_exact(source.imposemin_z, lite.imposemin_z, 'Minimum-slope DEM');
compare_close(source.smooth_z, lite.smooth_z, 1e-10, 'CRS-smoothed DEM');
compare_exact(source.stream_x, lite.stream_x, 'Stream coordinates');
compare_exact(source.stream_y, lite.stream_y, 'Stream coordinates');
compare_exact(source.mapstruct_x, lite.mapstruct_x, 'Stream shapefile coordinates');
compare_exact(source.mapstruct_y, lite.mapstruct_y, 'Stream shapefile coordinates');
compare_exact(source.grid_size, lite.grid_size, 'Raster dimensions');
compare_exact(source.grid_refmat, lite.grid_refmat, 'Raster reference matrix');
compare_exact(source.grid_cellsize, lite.grid_cellsize, 'Raster cell size');

% The release must work with no source snapshot on the MATLAB path.
restoredefaultpath;
addpath(fullfile(model_root, 'HydroPol2D_Functions'));
addpath(release_root);
runtime = hydropol2d_add_runtime_paths(model_root);
symbols = {'GRIDobj', 'FLOWobj', 'STREAMobj', 'fillsinks'};
for i = 1:numel(symbols)
    resolved = which(symbols{i});
    assert(startsWith(resolved, runtime.topotoolbox_lite_root), ...
        '%s resolved outside the bundled runtime: %s', symbols{i}, resolved);
    assert(~contains(resolved, source_root), ...
        '%s resolved from the source snapshot: %s', symbols{i}, resolved);
end

trace_files = trace_topotoolbox_lite_runtime(lite_root, dem_path);
assert(~isempty(trace_files), 'The curated terrain trace did not reach any files.');
assert(all(startsWith(trace_files, lite_root)), ...
    'A traced terrain dependency resolved outside the curated runtime.');

results = struct();
results.status = "passed";
results.model_root = model_root;
results.source_runtime = source_root;
results.bundled_runtime = lite_root;
results.traced_file_count = numel(trace_files);
results.trace_files = trace_files;
results.checked_at = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss Z');
disp(results);
end

function signature = terrain_signature(model_root, runtime_root, dem_path)
setappdata(groot, 'HydroPol2DReleaseModelRoot', model_root);
setappdata(groot, 'HydroPol2DReleaseRuntimeRoot', runtime_root);
setappdata(groot, 'HydroPol2DReleaseDemPath', dem_path);
restoredefaultpath;
clear classes;
model_root = getappdata(groot, 'HydroPol2DReleaseModelRoot');
runtime_root = getappdata(groot, 'HydroPol2DReleaseRuntimeRoot');
dem_path = getappdata(groot, 'HydroPol2DReleaseDemPath');
addpath(fullfile(model_root, 'HydroPol2D_Functions'));
addpath(runtime_root);

DEM = GRIDobj(dem_path);
DEM_nan_border = DEM;
DEM_nan_border.Z(1, :) = nan;
DEM_crop = crop(DEM_nan_border);
DEM_coarse = resample(DEM_crop, DEM_crop.cellsize * 2, 'bilinear');
DEM_resampled = resample(DEM_coarse, DEM_crop, 'nearest');
mask = DEM_crop;
mask.Z = ~isnan(mask.Z);
DEM_clipped = clip(DEM_resampled, mask);

DEM_filled = fillsinks(DEM_clipped);
FD = FLOWobj(DEM_filled, 'preprocess', 'none');
A = flowacc(FD);
S = STREAMobj(FD, 'minarea', 1);
S = klargestconncomps(S, 1);
S = trunk(S);
DEM_min = imposemin(S, DEM_filled, 1e-4);
[x, y] = STREAMobj2XY(S);
M = STREAMobj2mapstruct(S);

% This follows the constrained regularized smoothing calculation used by
% DEM_smoothening, with split=0 to keep the release check deterministic and
% avoid creating a parallel pool for this small reference DEM.
FD_smooth = FLOWobj(DEM_filled, 'preprocess', 'fill');
S_smooth = STREAMobj(FD_smooth, 'minarea', 1);
S_smooth = klargestconncomps(S_smooth, 1);
S_smooth = trunk(S_smooth);
zs = crs(S_smooth, DEM_filled, 'K', 1, 'tau', 0.2, 'split', 0);
DEM_smooth = DEM_filled;
DEM_smooth.Z(S_smooth.IXgrid) = zs;

signature = struct();
signature.crop_z = DEM_crop.Z;
signature.resample_z = DEM_resampled.Z;
signature.filled_z = DEM_filled.Z;
signature.flowacc_z = A.Z;
signature.receiver_ix = FD.ix;
signature.receiver_ixc = FD.ixc;
signature.stream_ixgrid = S.IXgrid;
signature.imposemin_z = DEM_min.Z;
signature.smooth_z = DEM_smooth.Z;
signature.stream_x = x;
signature.stream_y = y;
signature.mapstruct_x = {M.X};
signature.mapstruct_y = {M.Y};
signature.grid_size = DEM.size;
signature.grid_refmat = DEM.refmat;
signature.grid_cellsize = DEM.cellsize;
end

function compare_exact(actual, expected, label)
assert(isequaln(actual, expected), '%s changed between runtimes.', label);
end

function compare_close(actual, expected, tolerance, label)
assert(isequal(size(actual), size(expected)), '%s changed array size.', label);
max_error = max(abs(actual(:) - expected(:)), [], 'omitnan');
assert(isempty(max_error) || max_error <= tolerance, ...
    '%s exceeded tolerance: %.16g > %.16g.', label, max_error, tolerance);
end
