function Inundation_Maps_Groundwater_Only(Paths, flags, running_control, max_GW_depth, ...
    Soil_Properties, DEM_raster, Elevation_Properties, GIS_data, Wshed_Properties, ...
    saver_memory_maps, idx_nan)
%INUNDATION_MAPS_GROUNDWATER_ONLY Export HydroPol-style groundwater-depth animation.

close all;

[Spectrum, ~, ~, ~, ~, ~, ~, ~, ~, ~] = coloramps();

OUT = struct();
OUT.ROOT = Paths.Root;
OUT.VIDEOS = Paths.Anim;
OUT.GIFS = Paths.Anim;
OUT.STATIC = Paths.FigPDF;
OUT.FIG = Paths.FigFIG;

mk(OUT.ROOT);
mk(OUT.VIDEOS);
mk(OUT.GIFS);
mk(OUT.STATIC);
mk(OUT.FIG);

GLOBAL_VIDEO = struct();
GLOBAL_VIDEO.VISIBLE = false;
GLOBAL_VIDEO.FPS = 4;
GLOBAL_VIDEO.TARGET_HEIGHT_PX = 1024;
GLOBAL_VIDEO.AVI_QUALITY = 95;
GLOBAL_VIDEO.CONVERT_TO_MP4 = true;
GLOBAL_VIDEO.MP4_CRF = 23;
GLOBAL_VIDEO.MP4_PRESET = 'medium';
GLOBAL_VIDEO.DELETE_AVI_AFTER_MP4 = false;
GLOBAL_VIDEO.AVI_PROFILE = 'Motion JPEG AVI';
GLOBAL_VIDEO.FORCE_CONST_FRAME_SIZE = false;

DEM_maps = gather_if_needed(Elevation_Properties.elevation_cell);
grid_resolution = Wshed_Properties.Resolution;
xmax = size(DEM_maps, 2);
ymax = size(DEM_maps, 1);
x_grid = GIS_data.xulcorner + grid_resolution * (1:xmax);
y_grid = GIS_data.yulcorner - grid_resolution * (1:ymax);

flag_export_groundwater_maps = isfield(flags, 'flag_export_groundwater_maps') && ...
    flags.flag_export_groundwater_maps == 1;
if ~(flags.flag_groundwater_modeling == 1 && flag_export_groundwater_maps)
    error('Groundwater animation requested, but groundwater export flag is disabled.');
end

zmax = max(max_GW_depth(:));
zmin = min(max_GW_depth(:));
if ~isfinite(zmin)
    zmin = 0;
end
if ~isfinite(zmax) || zmax <= zmin
    zmax = zmin + 1;
end

make_static_max_depth_figure(OUT, Spectrum, x_grid, y_grid, max_GW_depth, idx_nan, zmin, zmax);

baseName = 'GW_Depths';
[video, aviPath, mp4Path, fig] = start_global_video(OUT.VIDEOS, baseName, GLOBAL_VIDEO);
clf(fig);
set(fig, 'DefaultTextInterpreter', 'latex');
set(fig, 'Color', 'w');

ax = axes('Parent', fig);
hold(ax, 'on');
cb = colorbar(ax);
cb.Label.String = 'GW Depth [m]';
cb.Label.Interpreter = 'latex';
cb.FontName = 'Garamond';
cb.FontSize = 12;
cb.TickDirection = 'out';

store = 1;
flag_loader = true;
targetH = [];
targetW = [];

for t = 1:length(running_control.time_records)
    cla(ax);

    if t > saver_memory_maps * store
        store = store + 1;
        load(fullfile('Temporary_Files', sprintf('save_map_hydro_%d', store)), 'Maps');
        flag_loader = false;
    elseif flag_loader
        load(fullfile('Temporary_Files', sprintf('save_map_hydro_%d', store)), 'Maps');
        flag_loader = false;
    end

    local_t = t - (store - 1) * saver_memory_maps;
    if ~isfield(Maps, 'Hydro') || ~isfield(Maps.Hydro, 'GWdepth_save') || ...
            isempty(Maps.Hydro.GWdepth_save) || local_t > size(Maps.Hydro.GWdepth_save, 3)
        continue;
    end

    GW = local_get_groundwater_depth(Maps.Hydro, local_t, Soil_Properties);
    GW(GW < 0) = NaN;
    GW(idx_nan) = NaN;

    surf(ax, x_grid, y_grid, GW, 'EdgeColor', 'none');
    shading(ax, 'interp');
    view(ax, 0, 90);
    axis(ax, [min(x_grid) max(x_grid) min(y_grid) max(y_grid) zmin zmax]);
    colormap(ax, Spectrum);
    caxis(ax, [zmin zmax]);

    t_title = running_control.time_records(t);
    if isa(t_title, 'datetime')
        time_str = sprintf('Time = %s', datestr(t_title, 'dd-mmm-yyyy HH:MM'));
    else
        if isfield(flags, 'flag_elapsed_time') && flags.flag_elapsed_time
            time_str = sprintf('Time [h] = %.2f', t_title / 60);
        else
            time_str = sprintf('Time [min] = %.2f', t_title);
        end
    end

    title(ax, time_str, 'Interpreter', 'latex', 'FontSize', 14);
    xlabel(ax, 'Easting [m]', 'Interpreter', 'latex');
    ylabel(ax, 'Northing [m]', 'Interpreter', 'latex');
    zlabel(ax, 'GW Depth [m]', 'Interpreter', 'latex');
    set(ax, 'FontName', 'Garamond', 'FontSize', 12, 'LineWidth', 2, 'TickDir', 'out');
    xtickformat(ax, '%.0f');
    ytickformat(ax, '%.0f');
    box(ax, 'on');

    drawnow;
    fr = getframe(fig);
    img = fr.cdata;
    if GLOBAL_VIDEO.FORCE_CONST_FRAME_SIZE
        if isempty(targetH)
            [targetH, targetW, ~] = size(img);
        elseif size(img,1) ~= targetH || size(img,2) ~= targetW
            img = safe_imresize(img, [targetH targetW]);
        end
    end
    writeVideo(video, img);
end

finish_global_video(video, aviPath, mp4Path, GLOBAL_VIDEO);
close(fig);
close all;
end

function make_static_max_depth_figure(OUT, Spectrum, x_grid, y_grid, max_GW_depth, idx_nan, zmin, zmax)
F = gather_if_needed(max_GW_depth);
F(F < 0) = NaN;
F(idx_nan) = NaN;

fig = figure('Units', 'inches', 'Position', [2, 2, 8, 6], 'Visible', 'off', 'Color', 'w');
surf(x_grid, y_grid, F, 'EdgeColor', 'none');
shading interp;
axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid) zmin zmax]);
view(0, 90);
title('Maximum Groundwater Depth', 'Interpreter', 'latex', 'FontSize', 12);
colormap(Spectrum);
k = colorbar;
k.FontName = 'Garamond';
k.FontSize = 12;
k.TickDirection = 'out';
ylabel(k, 'GW Depth [m]', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Easting [m]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Northing [m]', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('GW Depth [m]', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontName', 'Garamond', 'FontSize', 12, 'LineWidth', 2);
xtickformat('%.0f');
ytickformat('%.0f');
box on;
exportgraphics(fig, fullfile(OUT.STATIC, 'Max_GW_Depth.png'), ...
    'ContentType', 'image', 'Colorspace', 'rgb', 'Resolution', 300);
saveas(fig, fullfile(OUT.FIG, 'Max_GW_Depth.fig'));
close(fig);
end

function [video, aviPath, mp4Path, fig] = start_global_video(folderOut, baseName, GV)
defaultAspect = 7 / 12;
targetH = GV.TARGET_HEIGHT_PX;
targetW = max(2, round(targetH * defaultAspect)); %#ok<NASGU>

fig = figure( ...
    'Color', 'w', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Resize', 'off', ...
    'Renderer', 'opengl', ...
    'Visible', on_off(GV.VISIBLE));

aviPath = fullfile(folderOut, [baseName '.avi']);
mp4Path = fullfile(folderOut, [baseName '.mp4']);
safe_delete(aviPath);
safe_delete(mp4Path);

video = VideoWriter(aviPath, GV.AVI_PROFILE);
video.FrameRate = GV.FPS;
try
    video.Quality = GV.AVI_QUALITY;
catch
end
open(video);
end

function finish_global_video(video, aviPath, mp4Path, GV)
close(video);
if GV.CONVERT_TO_MP4 && ffmpeg_exists()
    ok = convert_avi_to_mp4_ffmpeg(aviPath, mp4Path, GV.FPS, GV.MP4_CRF, GV.MP4_PRESET);
    if ok && GV.DELETE_AVI_AFTER_MP4
        safe_delete(aviPath);
    end
end
end

function tf = ffmpeg_exists()
[status, ~] = system('ffmpeg -version');
tf = (status == 0);
end

function ok = convert_avi_to_mp4_ffmpeg(aviPath, mp4Path, fps, crf, preset)
ok = false;
if ~exist(aviPath, 'file')
    return;
end
cmd = sprintf('ffmpeg -y -hide_banner -loglevel error -r %g -i "%s" -c:v libx264 -pix_fmt yuv420p -crf %d -preset %s "%s"', ...
    fps, aviPath, crf, preset, mp4Path);
[status, ~] = system(cmd);
ok = (status == 0) && exist(mp4Path, 'file') == 2;
end

function v = gather_if_needed(v)
if isa(v, 'gpuArray')
    v = gather(v);
end
end

function s = on_off(tf)
if tf
    s = 'on';
else
    s = 'off';
end
end

function mk(folderPath)
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
end

function safe_delete(filePath)
if exist(filePath, 'file') == 2
    delete(filePath);
end
end

function img = safe_imresize(img, targetSize)
img = imresize(img, targetSize);
end

function GW = local_get_groundwater_depth(HydroMaps, local_t, Soil_Properties)
GW = gather_if_needed(HydroMaps.GWdepth_save(:, :, local_t));
is_new_zwt = isfield(HydroMaps, 'GWdepth_is_zwt') && isequal(HydroMaps.GWdepth_is_zwt, 1);
if ~is_new_zwt
    GW = Soil_Properties.Soil_Depth - GW;
end
GW = max(GW, 0);
if isfield(Soil_Properties, 'Soil_Depth')
    GW = min(GW, Soil_Properties.Soil_Depth);
end
end
