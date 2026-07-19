function review = review_topotoolbox_lite_dependencies(model_root)
%REVIEW_TOPOTOOLBOX_LITE_DEPENDENCIES Review MATLAB's static dependencies.
%   The runtime trace is the release authority because class dependency
%   analysis can conservatively include unused methods. This utility records
%   that broader MATLAB view so a release review can inspect any extra files.

arguments
    model_root (1,:) char = fileparts(fileparts(fileparts(mfilename('fullpath'))))
end

functions_root = fullfile(model_root, 'HydroPol2D_Functions');
runtime = hydropol2d_add_runtime_paths(model_root);

targets = fullfile(functions_root, {
    'HydroPol2D_preprocessing.m'
    'DEM_smoothening.m'
    'Plot_Initial_Maps.m'
    'post_processing.m'
    'NWP_rainfall_processing.m'
    'Satellite_rainfall_processing.m'
    'Automatic_Calibrator_HydroPol2D.m'});

[required_files, products] = matlab.codetools.requiredFilesAndProducts(targets);
required_files = string(required_files(:));
bundled_files = sort(required_files(startsWith(required_files, ...
    string(runtime.topotoolbox_lite_root) + filesep)));

manifest_path = fullfile(runtime.topotoolbox_lite_root, 'MANIFEST.txt');
manifest_lines = readlines(manifest_path);
manifest_rel = manifest_lines(~startsWith(manifest_lines, "#") & manifest_lines ~= "");
manifest_files = sort(string(fullfile(runtime.topotoolbox_lite_root, manifest_rel)));

review = struct();
review.targets = string(targets(:));
review.required_file_count = numel(required_files);
review.bundled_static_file_count = numel(bundled_files);
review.bundled_static_files = bundled_files;
review.not_listed_in_manifest = setdiff(bundled_files, manifest_files);
review.manifest_files_not_found_by_static_analysis = setdiff(manifest_files, bundled_files);
review.products = products;
review.note = "Static analysis is conservative for MATLAB class methods. " + ...
    "Use trace_topotoolbox_lite_runtime for the authoritative runtime closure.";

disp(review);
end
