function files = HydroPol2D_Export_Voronoi_Rasters(mesh_file,overlap_file,results,output_directory,prefix)
%HYDROPOL2D_EXPORT_VORONOI_RASTERS Conservatively remap saved polygon states.

if exist('geotiffwrite','file') ~= 2 || exist('maprefcells','file') ~= 2
    error('HydroPol2D:MissingMappingToolbox', ...
        'Voronoi GeoTIFF export requires geotiffwrite and maprefcells.');
end
if exist(output_directory,'dir') ~= 7, mkdir(output_directory); end
mapping = HydroPol2D_Read_Overlap(overlap_file, mesh_file);
mesh = HydroPol2D_Read_UGRID(mesh_file);
assert(size(mapping.mesh_to_raster,2)==mesh.n_cells, ...
    'HydroPol2D:InvalidVoronoiCase','Overlap weights do not match the mesh.');

native = struct();
native.peak_surface_depth_m = max(results.surface_depth_m,[],2);
native.final_surface_depth_m = results.surface_depth_m(:,end);
native.peak_surface_velocity_m_s = max(results.surface_velocity_m_s,[],2);
native.final_surface_velocity_m_s = results.surface_velocity_m_s(:,end);
if ~isempty(results.surface_velocity_x_m_s)
    native.final_surface_velocity_x_m_s = results.surface_velocity_x_m_s(:,end);
    native.final_surface_velocity_y_m_s = results.surface_velocity_y_m_s(:,end);
end
if ~isempty(fieldnames(results.hydrology))
    native.cumulative_infiltration_m = results.hydrology.cumulative_infiltration_m(:,end);
    native.soil_water_m = results.hydrology.soil_water_m(:,end);
    native.cumulative_actual_et_m = results.hydrology.cumulative_actual_et_m(:,end);
    native.cumulative_recharge_m = results.hydrology.cumulative_recharge_m(:,end);
end
if ~isempty(results.groundwater_head_m)
    native.groundwater_head_m = results.groundwater_head_m(:,end);
    native.depth_to_groundwater_m = max(mesh.surface_bed(:)-results.groundwater_head_m(:,end),0);
    native.groundwater_exchange_m_s = results.groundwater_river_exchange_m_s(:,end);
end

R = maprefcells([mapping.x_edges(1) mapping.x_edges(end)], ...
    [mapping.y_edges(1) mapping.y_edges(end)], ...
    [mapping.raster_rows mapping.raster_cols],'ColumnsStartFrom','north');
geotiff_options={};
if ~isempty(mesh.crs_wkt) && exist('projcrs','file')==2
    try
        projected=projcrs(mesh.crs_wkt);
        if ~isempty(projected.AuthorityCode) && projected.AuthorityCode>0
            geotiff_options={'CoordRefSysCode',projected.AuthorityCode};
        end
    catch
        warning('HydroPol2D:UnknownVoronoiCRS', ...
            'GeoTIFF transform was written, but the UGRID CRS has no recognized EPSG code.');
    end
end
names = fieldnames(native); files = strings(numel(names),1); mapped = struct();
for k=1:numel(names)
    % Every field in `native` is INTENSIVE -- depths, velocities, heads,
    % cumulative metres -- so the conservative remap would scale each partly
    % covered raster cell down by its coverage fraction.  See
    % hp2d_remap_intensive for the measured symptom.
    flat = hp2d_remap_intensive(mapping, native.(names{k})(:), 0.5, NaN);
    map = reshape(flat,[mapping.raster_cols mapping.raster_rows])';
    mapped.(names{k}) = flipud(map);
    files(k)=fullfile(output_directory,[prefix '-' strrep(names{k},'_','-') '.tif']);
    geotiffwrite(char(files(k)),single(mapped.(names{k})),R,geotiff_options{:});
end
save(fullfile(output_directory,[prefix '-remapped-maps.mat']),'native','mapped','files');

panel_names = intersect({'peak_surface_depth_m','peak_surface_velocity_m_s', ...
    'cumulative_infiltration_m','cumulative_actual_et_m', ...
    'groundwater_head_m','depth_to_groundwater_m'},names,'stable');
if ~isempty(panel_names)
    figure_handle=figure('Visible','off','Color','w','Position',[50 50 1200 850]);
    layout=tiledlayout(figure_handle,3,2,'TileSpacing','compact','Padding','compact');
    for k=1:numel(panel_names)
        axes_handle=nexttile(layout); imagesc(axes_handle,mapped.(panel_names{k}));
        axis(axes_handle,'image','off'); colormap(axes_handle,parula(256));
        colorbar(axes_handle); title(axes_handle,strrep(panel_names{k},'_',' '), ...
            'FontName','Helvetica','FontWeight','normal','FontSize',10);
        set(axes_handle,'FontName','Helvetica','LineWidth',1.5);
    end
    exportgraphics(figure_handle,fullfile(output_directory,[prefix '-state-summary.png']),'Resolution',300);
    exportgraphics(figure_handle,fullfile(output_directory,[prefix '-state-summary.pdf']),'ContentType','vector');
    close(figure_handle);
end
end
