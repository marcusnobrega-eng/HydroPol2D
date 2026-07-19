function manifest = build_topotoolbox_lite_manifest(model_root, source_root, output_file)
%BUILD_TOPOTOOLBOX_LITE_MANIFEST Build the vendored TopoToolbox file list.
%   This release-preparation utility combines a runtime trace of the active
%   HydroPol2D terrain workflow with a small set of flag-dependent branches.
%   MATLAB's generic class dependency analysis expands a class reference into
%   unrelated TopoToolbox suites, so it is retained only as a review aid and
%   not used as the release file list.

arguments
    model_root (1,:) char = fileparts(fileparts(fileparts(mfilename('fullpath'))))
    source_root (1,:) char = ''
    output_file (1,:) char = fullfile(model_root, 'third_party', ...
        'topotoolbox_lite', 'MANIFEST.txt')
end

if isempty(source_root)
    error('HydroPol2D:TopoToolboxLite:SourceRequired', ...
        'Provide a local TopoToolbox v2.4 source snapshot to rebuild the manifest.');
end
if ~isfolder(source_root)
    error('HydroPol2D:TopoToolboxLite:MissingSource', ...
        'TopoToolbox source snapshot not found: %s', source_root);
end

% Flag-dependent methods are retained even when their branch is not taken by
% the representative DEM. The list is limited to the HydroPol2D terrain API
% and its immediate dynamic fallbacks.
topotoolbox_roots = {
    '@GRIDobj/GRIDobj.m'
    '@GRIDobj/crop.m'
    '@GRIDobj/resample.m'
    '@GRIDobj/fillsinks.m'
    '@STREAMobj/imposemin.m'
    '@GRIDobj/arcslope.m'
    '@GRIDobj/hillshade.m'
    '@GRIDobj/clip.m'
    '@GRIDobj/imageschs.m'
    '@GRIDobj/shufflelabel.m'
    '@GRIDobj/GRIDobj2geotiff.m'
    '@FLOWobj/FLOWobj.m'
    '@FLOWobj/flowacc.m'
    '@FLOWobj/drainagebasins.m'
    '@STREAMobj/STREAMobj.m'
    '@STREAMobj/STREAMobj2mapstruct.m'
    '@STREAMobj/STREAMobj2XY.m'
    '@STREAMobj/klargestconncomps.m'
    '@STREAMobj/trunk.m'
    '@STREAMobj/crs.m'
    '@STREAMobj/STREAMobj2cell.m'
    '@STREAMobj/quantcarve.m'
    '@STREAMobj/modify.m'
    '@FLOWobj/validatealignment.m'
    '@GRIDobj/minus.m'
    '@GRIDobj/plus.m'
    '@GRIDobj/mtimes.m'
    '@GRIDobj/mrdivide.m'
    '@GRIDobj/rdivide.m'
    };

root_files = fullfile(source_root, topotoolbox_roots);
missing = root_files(~cellfun(@isfile, root_files));
if ~isempty(missing)
    error('HydroPol2D:TopoToolboxLite:MissingMethod', ...
        'Manifest source file not found: %s', strjoin(missing, ', '));
end

dem_path = fullfile(model_root, 'Static', 'DEM.tif');
trace_files = trace_topotoolbox_lite_runtime(source_root, dem_path);
source_root = string(source_root);

% MATLAB's static analysis closes over class methods that a representative
% terrain trace cannot always reach. Analyse the currently bundled runtime,
% then express those paths relative to the supplied source snapshot. This
% prevents an older source snapshot from hiding methods required by the
% released HydroPol2D workflows.
runtime = hydropol2d_add_runtime_paths(model_root);
bundled_root = string(runtime.topotoolbox_lite_root);
static_targets = fullfile(model_root, 'HydroPol2D_Functions', {
    'HydroPol2D_preprocessing.m'
    'DEM_smoothening.m'
    'Plot_Initial_Maps.m'
    'post_processing.m'
    'NWP_rainfall_processing.m'
    'Satellite_rainfall_processing.m'
    'Automatic_Calibrator_HydroPol2D.m'});
[static_files, ~] = matlab.codetools.requiredFilesAndProducts(static_targets);
static_files = string(static_files(:));
static_files = static_files(startsWith(static_files, bundled_root + filesep));
static_files = fullfile(source_root, erase(static_files, bundled_root + filesep));
missing_static = static_files(~isfile(static_files));
if ~isempty(missing_static)
    error('HydroPol2D:TopoToolboxLite:StaticClosureMissing', ...
        'Source snapshot is missing a required static dependency: %s', ...
        strjoin(missing_static, ', '));
end

root_paths = strings(numel(topotoolbox_roots), 1);
for i = 1:numel(topotoolbox_roots)
    root_paths(i) = string(fullfile(source_root, topotoolbox_roots{i}));
end
manifest = unique([trace_files(:); root_paths; static_files]);
manifest = sort(manifest);

relative_files = erase(manifest, source_root + filesep);
output_folder = fileparts(output_file);
if ~isfolder(output_folder)
    mkdir(output_folder);
end

fid = fopen(output_file, 'w');
if fid < 0
    error('HydroPol2D:TopoToolboxLite:ManifestWriteFailed', ...
        'Cannot write manifest: %s', output_file);
end
cleanup_file = onCleanup(@() fclose(fid));

fprintf(fid, '# HydroPol2D TopoToolbox Lite manifest\n');
fprintf(fid, '# Upstream: TopoToolbox v2.4 (14-Jun-2022)\n');
fprintf(fid, '# Base source: https://github.com/wschwanghart/topotoolbox/tree/2.4\n');
fprintf(fid, '# Base tag commit: cdc21040d86bdebc4723bc28f689344db471410b\n');
fprintf(fid, '# Compatibility files: crop.m=c67012b, modify.m=bf442e3, quantcarve.m=164b441\n');
fprintf(fid, '# Generated: %s\n', char(datetime('now', ...
    'Format', 'yyyy-MM-dd HH:mm:ss Z')));
fprintf(fid, '# Files: %d\n\n', numel(relative_files));
fprintf(fid, '%s\n', relative_files);
end
