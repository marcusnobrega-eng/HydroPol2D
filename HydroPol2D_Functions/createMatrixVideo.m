function createMatrixVideo(flags, running_control, saver_memory_maps, folderName, FileName_String, Colormap_matrix, no_plot, A, RA, S_p, x_grid, y_grid, xbegin, xend, ybegin, yend, Elevation_Properties, Matrix_3D)
    % createSpatialRainfallVideo generates a video showing water level profile over time.
    %
    % Inputs:
    %   flags - Structure with control flags
    %   running_control - Structure with time records
    %   saver_memory_maps - Memory saving parameter
    %   folderName - Folder to save video
    %   FileName_String - String for filename
    %   Colormap_matrix - Colormap
    %   no_plot - Flag to disable plotting
    %   A, RA, S_p - Map data for overlays
    %   x_grid, y_grid - Grid coordinates
    %   xbegin, xend, ybegin, yend - Plot limits
    %   Elevation_Properties - Structure with elevation data
    %   Matrix_3D - Structure with rainfall parameters
    %
    close all;
    Video_Name = 'Rainfall_Intensities.mp4';
    video = VideoWriter(fullfile(folderName, Video_Name), 'MPEG-4');
    video.FrameRate = 2;
    open(video);
    h = figure;
    set(gcf, 'units', 'inches', 'position', [0, 0, 7, 12]);
    set(gcf, 'DefaultTextInterpreter', 'latex');
    
    store = 1;
    flag_loader = 1;
    Max_rains = -inf;
    Min_rains = inf;
    
    for i = 1:length(running_control.time_records)
        if i > saver_memory_maps * store
            store = store + 1;
            load(fullfile('Temporary_Files', sprintf('save_map_hydro_%d', store)), 'Maps');
        elseif flag_loader == 1
            load(fullfile('Temporary_Files', sprintf('save_map_hydro_%d', store)), 'Maps');
            flag_loader = 0;
        end
        Max_rains = max([Max_rains; max(Maps.Hydro.spatial_rainfall_maps, [], 'all')]);
        Min_rains = min([Min_rains; min(Maps.Hydro.spatial_rainfall_maps, [], 'all')]);
    end
    
    zmax = Max_rains;
    zmin = Min_rains;
    if zmin == zmax
        zmax = zmin + 10;
    end
    
    store = 1;
    flag_loader = 1;
    
    for t = 1:length(running_control.time_records)
        clf;
        if t > saver_memory_maps * store
            store = store + 1;
            load(fullfile('Temporary_Files', sprintf('save_map_hydro_%d', store)), 'Maps');
        elseif flag_loader == 1
            load(fullfile('Temporary_Files', sprintf('save_map_hydro_%d', store)), 'Maps');
            flag_loader = 0;
        end
        
        rain = Maps.Hydro.spatial_rainfall_maps(:, :, t - (store - 1) * saver_memory_maps);
        idx = isnan(Elevation_Properties.elevation_cell);
        rain(rain <= 0) = nan;
        rain(idx) = nan;
        
        if t > length(Matrix_3D.rainfall_spatial_duration)
            t_title = Matrix_3D.rainfall_spatial_duration(end) + Matrix_3D.rainfall_spatial_duration_agg(2);
            z = zeros(size(Elevation_Properties.elevation_cell));
        else
            t_title = Matrix_3D.rainfall_spatial_duration(t);
            z = rain;
        end
        
        hold on;
        F = z(ybegin:yend, xbegin:xend);
        
        if no_plot == 0
            try
                mapshow(A, RA, "AlphaData", 0.45); hold on;
                mapshow(S_p, 'FaceColor', 'n'); hold on;
            catch
            end
        end
        
        surf(x_grid, y_grid, F);
        axis([min(x_grid(:)) max(x_grid(:)) min(y_grid(:)) max(y_grid(:)) zmin zmax]);
        shading INTERP;
        
        if flags.flag_elapsed_time == 1
            title(sprintf('Time [h] = %4.2f', t_title / 60), 'Interpreter', 'Latex', 'FontSize', 12);
        else
            title(sprintf(string(t_title)), 'Interpreter', 'Latex', 'FontSize', 12);
        end
        
        view(0, 90);
        colorbar;
        caxis([zmin zmax]);
        colormap(Colormap_matrix);
        k = colorbar;
        k.FontName = 'Garamond';
        k.FontSize = 12;
        k.TickDirection = 'out';
        ylabel(k, 'Rainfall (mm/h)', 'Interpreter', 'Latex', 'FontSize', 12);
        xlabel('Easting [m]', 'Interpreter', 'Latex', 'FontSize', 12);
        ylabel('Northing [m]', 'Interpreter', 'Latex', 'FontSize', 12);
        zlabel('Rainfall (mm/h)', 'Interpreter', 'Latex', 'FontSize', 12);
        
        set(gca, 'FontName', 'Garamond', 'FontSize', 12, 'LineWidth', 2, 'TickDir', 'out', 'TickLength', [0.02 0.01]);
        xtickformat('%.0f');
        ytickformat('%.0f');
        
        frame = getframe(gcf);
        writeVideo(video, frame);
        hold off;
        clf;
    end
    
    close(video);
    close all;
end
