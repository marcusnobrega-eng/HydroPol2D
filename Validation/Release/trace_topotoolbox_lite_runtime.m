function trace_files = trace_topotoolbox_lite_runtime(topotoolbox_root, dem_path)
%TRACE_TOPOTOOLBOX_LITE_RUNTIME Trace the active HydroPol2D terrain API.
%   The trace covers the TopoToolbox calls made by preprocessing, optional
%   DEM smoothing, spatial-forcing alignment, and map/shapefile export.
%   It returns only TopoToolbox files reached by this representative run.

arguments
    topotoolbox_root (1,:) char
    dem_path (1,:) char
end

if ~isfolder(topotoolbox_root)
    error('HydroPol2D:TopoToolboxLite:MissingRuntime', ...
        'TopoToolbox runtime not found: %s', topotoolbox_root);
end
if ~isfile(dem_path)
    error('HydroPol2D:TopoToolboxLite:MissingDEM', ...
        'Representative DEM not found: %s', dem_path);
end

old_path = path;
restore_path = onCleanup(@() path(old_path));
addpath(topotoolbox_root);

old_visibility = get(groot, 'DefaultFigureVisible');
restore_visibility = onCleanup(@() set(groot, ...
    'DefaultFigureVisible', old_visibility));
set(groot, 'DefaultFigureVisible', 'off');

profile off;
profile clear;
profile on -history;

DEM = GRIDobj(dem_path);
DEM_nan_border = DEM;
DEM_nan_border.Z(1, :) = nan;
DEM = crop(DEM_nan_border);
DEM_resampled = resample(DEM, DEM.cellsize * 2, 'bilinear');
mask = DEM;
mask.Z = ~isnan(mask.Z);
DEM_resampled = resample(DEM_resampled, DEM, 'nearest');
DEM_clipped = clip(DEM_resampled, mask);

DEM_filled = fillsinks(DEM_clipped);
FD = FLOWobj(DEM_filled, 'preprocess', 'none');
A = flowacc(FD); %#ok<NASGU>
S = STREAMobj(FD, 'minarea', 1);
D = drainagebasins(FD);

figure('Visible', 'off');
imageschs(DEM_filled, shufflelabel(D));
close(gcf);

if ~isempty(S.IXgrid)
    S = klargestconncomps(S, 1);
    S = trunk(S);
    DEM_min = imposemin(S, DEM_filled, 1e-4); %#ok<NASGU>
    [~, ~] = STREAMobj2XY(S);
    STREAMobj2mapstruct(S);
    figure('Visible', 'off');
    plot(S);
    close(gcf);
    if numel(S.IXgrid) >= 3
        try
            crs(S, DEM_filled, 'K', 1, 'tau', 0.5, 'split', 0);
        catch ME
            % Older upstream releases use a linprog signature that current
            % MATLAB no longer accepts. The manifest's explicit root list
            % still retains crs and quantcarve for review.
            warning('HydroPol2D:TopoToolboxLite:TraceCRS', '%s', ME.message);
        end
    end
end

arcslope(DEM_filled);
hillshade(DEM_filled);
GRIDobj2geotiff(DEM_filled, fullfile(tempdir, ...
    'hydropol2d_topotoolbox_lite_trace.tif'));

profile off;
profile_data = profile('info');
files = string({profile_data.FunctionTable.FileName});
root = string(topotoolbox_root);
trace_files = sort(unique(files(startsWith(files, root + filesep))));
end
