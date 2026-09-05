function HydroPol2D_Write_Voronoi_Output(output_file, mesh, results, overwrite)
%HYDROPOL2D_WRITE_VORONOI_OUTPUT Write the native unstructured-output archive.

arguments
    output_file (1,:) char
    mesh struct
    results struct
    overwrite (1,1) logical = false
end
if exist(output_file,'file')==2
    if ~overwrite
        error('HydroPol2D:OutputExists','Output exists and overwrite_output is false: %s',output_file);
    end
    delete(output_file);
end
parent=fileparts(output_file); if exist(parent,'dir')~=7, mkdir(parent); end
nt=numel(results.time_s);

write_vector(output_file,'time_s',results.time_s,'time_s','s');
write_vector(output_file,'map_time_s',results.map_time_s,'map_time_s','s');
write_vector(output_file,'face',(0:mesh.n_cells-1)','face','1');
write_vector(output_file,'edge',(0:mesh.n_edges-1)','edge','1');
if mesh.channel.n_nodes>0
    write_vector(output_file,'channel_node',(0:mesh.channel.n_nodes-1)', ...
        'channel_node','1');
end
if mesh.channel.n_links>0
    write_vector(output_file,'channel_link',(0:mesh.channel.n_links-1)', ...
        'channel_link','1');
end
write_matrix(output_file,'surface_depth_face_m',results.surface_depth_m,'face',mesh.n_cells,nt,'m');
write_matrix(output_file,'surface_velocity_face_m_s',results.surface_velocity_m_s,'face',mesh.n_cells,nt,'m s-1');
write_matrix(output_file,'edge_discharge_per_width_m2_s', ...
    results.edge_discharge_per_width_history_m2_s,'edge',mesh.n_edges,nt,'m2 s-1');
write_matrix(output_file,'surface_source_rate_face_m_s', ...
    results.surface_source_rate_m_s,'face',mesh.n_cells,nt,'m s-1');
write_matrix(output_file,'cumulative_surface_source_face_m', ...
    results.cumulative_surface_source_m,'face',mesh.n_cells,nt,'m');
write_matrix(output_file,'potential_et_rate_face_m_s', ...
    results.potential_et_rate_m_s,'face',mesh.n_cells,nt,'m s-1');
write_vector(output_file,'outlet_discharge_m3_s',results.outlet_discharge_m3_s,'time_s','m3 s-1');

if ~isempty(results.surface_momentum_x_m2_s)
    write_matrix(output_file,'surface_momentum_x_face_m2_s',results.surface_momentum_x_m2_s,'face',mesh.n_cells,nt,'m2 s-1');
    write_matrix(output_file,'surface_momentum_y_face_m2_s',results.surface_momentum_y_m2_s,'face',mesh.n_cells,nt,'m2 s-1');
    write_matrix(output_file,'surface_velocity_x_face_m_s',results.surface_velocity_x_m_s,'face',mesh.n_cells,nt,'m s-1');
    write_matrix(output_file,'surface_velocity_y_face_m_s',results.surface_velocity_y_m_s,'face',mesh.n_cells,nt,'m s-1');
end
write_vector(output_file,'final_surface_volume_face_m3',results.final_surface_volume_m3,'face','m3');
write_vector(output_file,'final_edge_discharge_per_width_m2_s', ...
    results.edge_discharge_per_width_m2_s,'edge','m2 s-1');

if mesh.channel.n_nodes>0
    write_matrix(output_file,'channel_depth_node_m',results.channel_depth_m, ...
        'channel_node',mesh.channel.n_nodes,nt,'m');
    write_vector(output_file,'final_channel_volume_node_m3',results.final_channel_volume_m3, ...
        'channel_node','m3');
end
if mesh.channel.n_links>0
    write_matrix(output_file,'channel_discharge_link_m3_s',results.channel_discharge_history_m3_s, ...
        'channel_link',mesh.channel.n_links,nt,'m3 s-1');
end
if mesh.channel.n_transitions>0
    write_matrix(output_file,'channel_transition_discharge_m3_s', ...
        results.transition_discharge_history_m3_s,'channel_transition', ...
        mesh.channel.n_transitions,nt,'m3 s-1');
end

if ~isempty(results.groundwater_head_m)
    write_matrix(output_file,'groundwater_head_face_m',results.groundwater_head_m,'face',mesh.n_cells,nt,'m');
    write_matrix(output_file,'groundwater_depth_face_m', ...
        max(mesh.surface_bed(:)-results.groundwater_head_m,0),'face',mesh.n_cells,nt,'m');
    write_matrix(output_file,'groundwater_river_exchange_face_m_s', ...
        results.groundwater_river_exchange_m_s,'face',mesh.n_cells,nt,'m s-1');
    write_matrix(output_file,'groundwater_seepage_face_m_s', ...
        results.groundwater_seepage_m_s,'face',mesh.n_cells,nt,'m s-1');
end
if ~isempty(fieldnames(results.hydrology))
    fields={'soil_water_m','canopy_storage_m','cumulative_infiltration_m', ...
        'cumulative_recharge_m','cumulative_actual_et_m','infiltration_rate_m_s', ...
        'recharge_rate_m_s','capillary_rate_m_s','actual_et_rate_m_s'};
    units={'m','m','m','m','m','m s-1','m s-1','m s-1','m s-1'};
    for k=1:numel(fields)
        write_matrix(output_file,[fields{k} '_face'],results.hydrology.(fields{k}), ...
            'face',mesh.n_cells,nt,units{k});
    end
    n_hru=size(results.hydrology.hru_fraction,2);
    hru_fields={'hru_fraction','final_canopy_storage_m','final_near_surface_storage_m', ...
        'final_root_zone_storage_m','final_transmission_storage_m'};
    for k=1:numel(hru_fields)
        name=hru_fields{k}; values=results.hydrology.(name);
        nccreate(output_file,name,'Dimensions',{'face',mesh.n_cells,'hru',n_hru},'Datatype','double');
        ncwrite(output_file,name,values);
    end
end

n_gauges=numel(results.gauge_names);
if n_gauges>0
    write_vector(output_file,'gauge',(0:n_gauges-1)','gauge','1');
    write_matrix(output_file,'gauge_discharge_m3_s',results.gauge_discharge_m3_s, ...
        'gauge',n_gauges,nt,'m3 s-1');
    write_matrix(output_file,'gauge_surface_depth_m',results.gauge_surface_depth_m, ...
        'gauge',n_gauges,nt,'m');
    write_matrix(output_file,'gauge_water_surface_elevation_m', ...
        results.gauge_water_surface_elevation_m,'gauge',n_gauges,nt,'m');
    ncwriteatt(output_file,'/','gauge_names',strjoin(cellstr(results.gauge_names),','));
end

if ~isempty(results.diagnostics)
    write_vector(output_file,'diagnostic_time_s',[results.diagnostics.time_s]','step','s');
    fields={'dt_s','mass_m3','step_mass_residual_m3','max_surface_depth_m', ...
        'max_surface_velocity_m_s','max_channel_depth_m','max_channel_velocity_m_s', ...
        'max_transition_velocity_m_s','max_groundwater_flux_m3_s','surface_channel_exchange_m3', ...
        'boundary_net_inflow_volume_m3','groundwater_substep_count', ...
        'groundwater_river_exchange_m3','groundwater_seepage_volume_m3', ...
        'precipitation_volume_m3','surface_source_volume_m3','infiltration_volume_m3', ...
        'actual_et_volume_m3','recharge_volume_m3','capillary_volume_m3', ...
        'hydrology_mass_residual_m3'};
    for k=1:numel(fields)
        write_vector(output_file,fields{k},[results.diagnostics.(fields{k})]','step','');
    end
end

% Store raster-time histories as a portable convenience. Polygon histories
% remain authoritative, and temporal maxima are still calculated only after
% remapping every saved native timestep during post-processing.
if ~isempty(results.overlap_file)
    mapping=HydroPol2D_Read_Overlap(results.overlap_file,results.mesh_file);
    assert(size(mapping.mesh_to_raster,2)==mesh.n_cells, ...
        'HydroPol2D:InvalidVoronoiOutput','Overlap weights do not match the mesh.');
    map_indices=nearest_indices(results.time_s(:),results.map_time_s(:));
    write_vector(output_file,'row',(0:mapping.raster_rows-1)','row','1');
    write_vector(output_file,'column',(0:mapping.raster_cols-1)','column','1');

    write_face_map('surface_depth_map_m',results.surface_depth_m,'m');
    write_face_map('surface_velocity_map_m_s',results.surface_velocity_m_s,'m s-1');
    write_face_map('water_surface_elevation_map_m', ...
        results.surface_depth_m+mesh.surface_bed(:),'m');
    write_face_map('rainfall_rate_map_m_s',results.surface_source_rate_m_s,'m s-1');
    write_face_map('cumulative_rainfall_map_m',results.cumulative_surface_source_m,'m');
    write_face_map('potential_et_rate_map_m_s',results.potential_et_rate_m_s,'m s-1');

    if ~isempty(results.groundwater_head_m)
        write_face_map('groundwater_head_map_m',results.groundwater_head_m,'m');
        write_face_map('groundwater_depth_map_m', ...
            max(mesh.surface_bed(:)-results.groundwater_head_m,0),'m');
        write_face_map('groundwater_river_exchange_map_m_s', ...
            results.groundwater_river_exchange_m_s,'m s-1');
        write_face_map('groundwater_seepage_map_m_s',results.groundwater_seepage_m_s,'m s-1');
    end
    if ~isempty(fieldnames(results.hydrology))
        map_fields={'soil_water_m','cumulative_infiltration_m','cumulative_recharge_m', ...
            'cumulative_actual_et_m','infiltration_rate_m_s','recharge_rate_m_s', ...
            'actual_et_rate_m_s'};
        map_names={'soil_water_map_m','cumulative_infiltration_map_m', ...
            'cumulative_recharge_map_m','cumulative_actual_et_map_m', ...
            'infiltration_rate_map_m_s','recharge_rate_map_m_s','actual_et_rate_map_m_s'};
        map_units={'m','m','m','m','m s-1','m s-1','m s-1'};
        for k=1:numel(map_fields)
            write_face_map(map_names{k},results.hydrology.(map_fields{k}),map_units{k});
        end
    end
    if mesh.channel.n_nodes>0 && ~isempty(results.channel_depth_m)
        face_values=zeros(mesh.n_cells,nt);
        face_values(mesh.channel.host_cell,:)=results.channel_depth_m;
        write_face_map('channel_depth_map_m',face_values,'m');
    end
    if mesh.channel.n_links>0 && ~isempty(results.channel_discharge_history_m3_s)
        face_values=zeros(mesh.n_cells,nt);
        for link=1:mesh.channel.n_links
            host=mesh.channel.host_cell(mesh.channel.link_down(link));
            face_values(host,:)=face_values(host,:)+results.channel_discharge_history_m3_s(link,:);
        end
        write_face_map('channel_discharge_map_m3_s',face_values,'m3 s-1');
    end
end

ncwriteatt(output_file,'/','Conventions','CF-1.13, UGRID-1.0');
ncwriteatt(output_file,'/','output_schema','hydropol2d-unstructured-output-1.0');
ncwriteatt(output_file,'/','mesh_file',results.mesh_file);
ncwriteatt(output_file,'/','overlap_file',results.overlap_file);
ncwriteatt(output_file,'/','mesh_file_sha256',file_sha256(results.mesh_file));
if ~isempty(results.overlap_file)
    ncwriteatt(output_file,'/','overlap_file_sha256',file_sha256(results.overlap_file));
end
ncwriteatt(output_file,'/','crs_wkt',mesh.crs_wkt);
ncwriteatt(output_file,'/','title','HydroPol2D native Voronoi finite-volume output');

    function write_face_map(name,face_history,units)
        if isempty(face_history), return; end
        maps=map_history(mapping,double(face_history(:,map_indices)));
        write_map_cube(output_file,name,maps,units);
    end
end

function write_matrix(path,name,values,space_dimension,n_space,n_time,units)
if isempty(values), return; end
assert(isequal(size(values),[n_space n_time]),'HydroPol2D:InvalidVoronoiOutput', ...
    '%s has inconsistent dimensions.',name);
% Native histories can be large; time-slice chunks keep appends and later
% post-processing efficient without changing numerical precision.
nccreate_arguments={'Dimensions',{space_dimension,n_space,'time_s',n_time}, ...
    'Datatype','double','ChunkSize',[min(n_space,65536) 1], ...
    'DeflateLevel',4,'Shuffle',true};
nccreate(path,name,nccreate_arguments{:});
ncwrite(path,name,values); if ~isempty(units), ncwriteatt(path,name,'units',units); end
end

function write_vector(path,name,values,dimension,units)
values=double(values(:));
nccreate(path,name,'Dimensions',{dimension,numel(values)},'Datatype','double');
ncwrite(path,name,values); if ~isempty(units), ncwriteatt(path,name,'units',units); end
end

function write_map_cube(path,name,maps,units)
% MATLAB reverses dimension order in the NetCDF file. This declaration and
% permutation expose the portable order (map_time_s, row, column).
values=permute(double(maps),[2 1 3]);
nccreate(path,name,'Dimensions',{'column',size(values,1),'row',size(values,2), ...
    'map_time_s',size(values,3)},'Datatype','double', ...
    'ChunkSize',[min(size(values,1),256) min(size(values,2),256) 1], ...
    'DeflateLevel',4,'Shuffle',true);
ncwrite(path,name,values); if ~isempty(units), ncwriteatt(path,name,'units',units); end
end

function indices=nearest_indices(native_times,map_times)
indices=zeros(numel(map_times),1);
for k=1:numel(map_times)
    [difference,indices(k)]=min(abs(native_times-map_times(k)));
    assert(difference<=1e-8*max(abs(map_times(k)),1), ...
        'HydroPol2D:InvalidVoronoiOutput','map_time_s is not on the native output schedule.');
end
end

function maps=map_history(mapping,face_history)
flat=mapping.mesh_to_raster*double(face_history); nt=size(flat,2);
maps=zeros(mapping.raster_rows,mapping.raster_cols,nt);
for k=1:nt
    maps(:,:,k)=flipud(reshape(flat(:,k),[mapping.raster_cols mapping.raster_rows])');
end
end

function digest=file_sha256(path)
quoted=['''' strrep(char(path),'''','''"''"''') ''''];
[status,text]=system(['sha256sum ' quoted]);
if status~=0, [status,text]=system(['shasum -a 256 ' quoted]); end
if status~=0, error('HydroPol2D:HashFailure','Could not hash %s.',path); end
parts=strsplit(strtrim(text)); digest=parts{1};
end
