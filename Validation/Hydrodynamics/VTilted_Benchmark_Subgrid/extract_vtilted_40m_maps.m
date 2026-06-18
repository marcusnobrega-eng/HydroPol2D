clear; clc;

case_dir = fileparts(mfilename('fullpath'));
out_root = fullfile(case_dir, 'Outputs', 'VTilted_40m_Comparison');
static_root = fullfile(case_dir, 'Static_40m_Comparison');
map_dir = fullfile(out_root, 'Maps');
if ~exist(map_dir, 'dir'); mkdir(map_dir); end

scenarios = ["reference20m", "baseline40m", "subgrid40m"];
for s = scenarios
    temp_dir = fullfile(out_root, char(s), 'HydroPol2D_Output', 'Temporary_Files');
    files = dir(fullfile(temp_dir, 'save_map_hydro_*.mat'));
    if isempty(files)
        error('No map files found for %s in %s', s, temp_dir);
    end
    max_depth_m = [];
    for k = 1:numel(files)
        f = fullfile(files(k).folder, files(k).name);
        D = load(f, 'Maps');
        if ~isfield(D, 'Maps') || ~isfield(D.Maps, 'Hydro') || ~isfield(D.Maps.Hydro, 'd')
            error('Maps.Hydro.d not found in %s', f);
        end
        d = double(gather(D.Maps.Hydro.d)) ./ 1000;
        d(~isfinite(d)) = NaN;
        dmax = max(d, [], 3, 'omitnan');
        if isempty(max_depth_m)
            max_depth_m = dmax;
        else
            max_depth_m = max(max_depth_m, dmax);
        end
    end
    writematrix(max_depth_m, fullfile(map_dir, sprintf('%s_max_depth_native_m.csv', s)));
end

dem20 = double(readgeoraster(fullfile(static_root, 'Static_20m_Crop', 'DEM.tif')));
dem40 = double(readgeoraster(fullfile(static_root, 'Static_40m_Resampled', 'DEM.tif')));
dem20(~isfinite(dem20)) = NaN;
dem40(~isfinite(dem40)) = NaN;

ref20 = readmatrix(fullfile(map_dir, 'reference20m_max_depth_native_m.csv'));
base40 = readmatrix(fullfile(map_dir, 'baseline40m_max_depth_native_m.csv'));
sub40 = readmatrix(fullfile(map_dir, 'subgrid40m_max_depth_native_m.csv'));

base40_projected = project_coarse_wse_to_fine(dem40 + base40, dem20, 2);
invert40 = block_min(dem20, 2);
sub40_projected = project_coarse_wse_to_fine(invert40 + sub40, dem20, 2);

writematrix(ref20, fullfile(map_dir, 'reference20m_max_depth_20m_m.csv'));
writematrix(base40_projected, fullfile(map_dir, 'baseline40m_projected_max_depth_20m_m.csv'));
writematrix(sub40_projected, fullfile(map_dir, 'subgrid40m_projected_max_depth_20m_m.csv'));
writematrix(base40_projected - ref20, fullfile(map_dir, 'baseline40m_projected_minus_reference_20m_m.csv'));
writematrix(sub40_projected - ref20, fullfile(map_dir, 'subgrid40m_projected_minus_reference_20m_m.csv'));

function projected = project_coarse_wse_to_fine(coarse_eta, fine_dem, ratio)
[nyc, nxc] = size(coarse_eta);
projected = zeros(size(fine_dem));
for i = 1:nyc
    rows = (i - 1) * ratio + (1:ratio);
    for j = 1:nxc
        cols = (j - 1) * ratio + (1:ratio);
        projected(rows, cols) = max(coarse_eta(i,j) - fine_dem(rows, cols), 0);
    end
end
end

function out = block_min(a, ratio)
out = NaN(size(a,1)/ratio, size(a,2)/ratio);
for i = 1:size(out,1)
    rows = (i - 1) * ratio + (1:ratio);
    for j = 1:size(out,2)
        cols = (j - 1) * ratio + (1:ratio);
        block = a(rows, cols);
        out(i,j) = min(block(:), [], 'omitnan');
    end
end
end
