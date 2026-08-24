function HydroPol2D_Write_Voronoi_Output(output_file, mesh, results, overwrite)
%HYDROPOL2D_WRITE_VORONOI_OUTPUT Write native polygon/channel states to NetCDF.

arguments
    output_file (1,:) char
    mesh struct
    results struct
    overwrite (1,1) logical = false
end
if exist(output_file, 'file') == 2
    if ~overwrite
        error('HydroPol2D:OutputExists', 'Output exists and overwrite_output is false: %s', output_file);
    end
    delete(output_file);
end

nccreate(output_file,'time_s','Dimensions',{'time',numel(results.time_s)},'Datatype','double');
ncwrite(output_file,'time_s',results.time_s);
nccreate(output_file,'surface_depth_m','Dimensions',{'cell',mesh.n_cells,'time',numel(results.time_s)},'Datatype','double');
ncwrite(output_file,'surface_depth_m',results.surface_depth_m);
nccreate(output_file,'final_surface_volume_m3','Dimensions',{'cell',mesh.n_cells},'Datatype','double');
ncwrite(output_file,'final_surface_volume_m3',results.final_surface_volume_m3);
nccreate(output_file,'final_edge_discharge_per_width_m2_s','Dimensions',{'edge',mesh.n_edges},'Datatype','double');
ncwrite(output_file,'final_edge_discharge_per_width_m2_s',results.edge_discharge_per_width_m2_s);
if mesh.channel.n_nodes > 0
    nccreate(output_file,'channel_depth_m','Dimensions',{'channel_node',mesh.channel.n_nodes,'time',numel(results.time_s)},'Datatype','double');
    ncwrite(output_file,'channel_depth_m',results.channel_depth_m);
    nccreate(output_file,'final_channel_volume_m3','Dimensions',{'channel_node',mesh.channel.n_nodes},'Datatype','double');
    ncwrite(output_file,'final_channel_volume_m3',results.final_channel_volume_m3);
end
if mesh.channel.n_links > 0
    nccreate(output_file,'final_channel_discharge_m3_s','Dimensions',{'channel_link',mesh.channel.n_links},'Datatype','double');
    ncwrite(output_file,'final_channel_discharge_m3_s',results.channel_discharge_m3_s);
end
if mesh.channel.n_transitions > 0
    nccreate(output_file,'final_transition_discharge_m3_s','Dimensions',{'channel_transition',mesh.channel.n_transitions},'Datatype','double');
    ncwrite(output_file,'final_transition_discharge_m3_s',results.channel_transition_discharge_m3_s);
end
if ~isempty(results.groundwater_head_m)
    nccreate(output_file,'groundwater_head_m','Dimensions',{'cell',mesh.n_cells,'time',numel(results.time_s)},'Datatype','double');
    ncwrite(output_file,'groundwater_head_m',results.groundwater_head_m);
    nccreate(output_file,'depth_to_groundwater_m','Dimensions',{'cell',mesh.n_cells,'time',numel(results.time_s)},'Datatype','double');
    ncwrite(output_file,'depth_to_groundwater_m',max(mesh.surface_bed(:)-results.groundwater_head_m,0));
    nccreate(output_file,'groundwater_river_exchange_m_s','Dimensions',{'cell',mesh.n_cells,'time',numel(results.time_s)},'Datatype','double');
    ncwrite(output_file,'groundwater_river_exchange_m_s',results.groundwater_river_exchange_m_s);
    nccreate(output_file,'groundwater_seepage_m_s','Dimensions',{'cell',mesh.n_cells,'time',numel(results.time_s)},'Datatype','double');
    ncwrite(output_file,'groundwater_seepage_m_s',results.groundwater_seepage_m_s);
end
if ~isempty(fieldnames(results.hydrology))
    cell_fields={'soil_water_m','canopy_storage_m','cumulative_infiltration_m', ...
        'cumulative_recharge_m','cumulative_actual_et_m','infiltration_rate_m_s', ...
        'recharge_rate_m_s','capillary_rate_m_s','actual_et_rate_m_s'};
    for k=1:numel(cell_fields)
        name=cell_fields{k};
        nccreate(output_file,name,'Dimensions',{'cell',mesh.n_cells,'time',numel(results.time_s)},'Datatype','double');
        ncwrite(output_file,name,results.hydrology.(name));
    end
    n_hru=size(results.hydrology.hru_fraction,2);
    hru_fields={'hru_fraction','final_canopy_storage_m','final_near_surface_storage_m', ...
        'final_root_zone_storage_m','final_transmission_storage_m'};
    for k=1:numel(hru_fields)
        name=hru_fields{k};
        nccreate(output_file,name,'Dimensions',{'cell',mesh.n_cells,'hru',n_hru},'Datatype','double');
        ncwrite(output_file,name,results.hydrology.(name));
    end
end
if ~isempty(results.diagnostics)
    values = [results.diagnostics.time_s]';
    nccreate(output_file,'diagnostic_time_s','Dimensions',{'step',numel(values)},'Datatype','double');
    ncwrite(output_file,'diagnostic_time_s',values);
    fields = {'dt_s','mass_m3','step_mass_residual_m3','max_surface_depth_m', ...
        'max_surface_velocity_m_s','max_channel_depth_m','max_channel_velocity_m_s', ...
        'max_transition_velocity_m_s','max_groundwater_flux_m3_s','surface_channel_exchange_m3', ...
        'boundary_net_inflow_volume_m3','groundwater_substep_count', ...
        'groundwater_river_exchange_m3','groundwater_seepage_volume_m3', ...
        'precipitation_volume_m3','infiltration_volume_m3','actual_et_volume_m3', ...
        'recharge_volume_m3','capillary_volume_m3','hydrology_mass_residual_m3'};
    for k = 1:numel(fields)
        name = fields{k}; values = [results.diagnostics.(name)]';
        nccreate(output_file,name,'Dimensions',{'step',numel(values)},'Datatype','double');
        ncwrite(output_file,name,values);
    end
end
ncwriteatt(output_file,'/','Conventions','CF-1.13, UGRID-1.0');
ncwriteatt(output_file,'/','mesh_file',results.mesh_file);
ncwriteatt(output_file,'/','title','HydroPol2D native Voronoi finite-volume output');
end
