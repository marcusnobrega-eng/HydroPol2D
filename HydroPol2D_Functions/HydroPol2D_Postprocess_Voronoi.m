function report = HydroPol2D_Postprocess_Voronoi(native_file,mesh_file,overlap_file,output_root,options)
%HYDROPOL2D_POSTPROCESS_VORONOI Conventional raster products from a native archive.
% The native polygon/channel histories remain authoritative. Temporal maxima
% are calculated only after every saved native state is conservatively remapped.

arguments
    native_file (1,:) char
    mesh_file (1,:) char
    overlap_file (1,:) char
    output_root (1,:) char
    options struct = struct()
end
options=defaults(options);
postprocess_clock=tic;
validate_inputs(native_file,mesh_file,overlap_file,options);
mesh=HydroPol2D_Read_UGRID(mesh_file); mapping=HydroPol2D_Read_Overlap(overlap_file,mesh_file);
assert(size(mapping.mesh_to_raster,2)==mesh.n_cells,'HydroPol2D:InvalidVoronoiOutput', ...
    'Overlap weights do not match the UGRID mesh.');
folders=create_folders(output_root);
time=read_required(native_file,'time_s'); time=time(:);
map_time=read_optional(native_file,'map_time_s',time); map_time=map_time(:);
map_indices=nearest_indices(time,map_time);
if ~isempty(options.raster_stack_interval_s)
    native_interval=median(diff(time)); ratio=options.raster_stack_interval_s/native_interval;
    assert(ratio>=1 && abs(ratio-round(ratio))<=1e-9*max(ratio,1), ...
        'HydroPol2D:InvalidVoronoiOutput', ...
        'raster_stack_interval_s must be an integer multiple of the native schedule.');
    map_indices=find(abs(time/options.raster_stack_interval_s- ...
        round(time/options.raster_stack_interval_s))<=1e-9); map_indices=unique([map_indices;numel(time)]);
    map_time=time(map_indices);
end

surface_depth=read_face_history(native_file,{'surface_depth_face_m','surface_depth_m'},mesh.n_cells);
velocity=read_face_history(native_file,{'surface_velocity_face_m_s','surface_velocity_m_s'},mesh.n_cells);
histories=struct(); metadata=struct();
[histories.surface_depth,metadata.surface_depth]=field(map_history(mapping,surface_depth), ...
    folders.water_depth,'Flood_Depth','m',false);
[histories.velocity,metadata.velocity]=field(map_history(mapping,velocity), ...
    folders.velocity,'Velocity','m_s',false);
wse=map_history(mapping,surface_depth+mesh.surface_bed(:));
[histories.wse,metadata.wse]=field(wse,folders.wse,'Water_Surface_Elevation','m',false);
[histories.depth_velocity,metadata.depth_velocity]=field(histories.surface_depth.*histories.velocity, ...
    folders.hazard,'Depth_Velocity_Index','m2_s',false);
instability=zeros(size(histories.depth_velocity));
instability(histories.depth_velocity>=0.4)=1; instability(histories.depth_velocity>=0.8)=2;
instability(histories.depth_velocity>=1.2)=3;
[histories.instability,metadata.instability]=field(instability,folders.hazard, ...
    'Instability_Class','class',false);

add_face('rainfall_rate',{'surface_source_rate_face_m_s'},folders.rainfall, ...
    'Rainfall_Intensity','mm_h',1000*3600,false);
add_face('rainfall_depth',{'cumulative_surface_source_face_m'},folders.rainfall, ...
    'Cumulative_Rainfall','mm',1000,false);
add_face('infiltration_rate',{'infiltration_rate_m_s_face','infiltration_rate_m_s'}, ...
    folders.infiltration,'Infiltration_Rate','mm_h',1000*3600,false);
add_face('infiltration_depth',{'cumulative_infiltration_m_face','cumulative_infiltration_m'}, ...
    folders.infiltration,'Cumulative_Infiltration','mm',1000,false);
add_face('soil_water',{'soil_water_m_face','soil_water_m'},folders.infiltration, ...
    'Soil_Water','mm',1000,false);
add_face('recharge_rate',{'recharge_rate_m_s_face','recharge_rate_m_s'},folders.recharge, ...
    'Recharge_Rate','mm_h',1000*3600,false);
add_face('recharge_depth',{'cumulative_recharge_m_face','cumulative_recharge_m'}, ...
    folders.recharge,'Cumulative_Recharge','mm',1000,false);
add_face('pet_rate',{'potential_et_rate_face_m_s'},folders.et,'Potential_ET','mm_day',1000*86400,false);
add_face('actual_et_rate',{'actual_et_rate_m_s_face','actual_et_rate_m_s'}, ...
    folders.et,'Actual_ET_Rate','mm_day',1000*86400,false);
add_face('actual_et_depth',{'cumulative_actual_et_m_face','cumulative_actual_et_m'}, ...
    folders.et,'Cumulative_Actual_ET','mm',1000,false);
add_face('groundwater_head',{'groundwater_head_face_m','groundwater_head_m'}, ...
    folders.groundwater,'Groundwater_Head','m',1,false);
add_face('groundwater_depth',{'groundwater_depth_face_m','depth_to_groundwater_m'}, ...
    folders.groundwater,'Depth_to_Groundwater','m',1,false);
add_face('groundwater_exchange',{'groundwater_river_exchange_face_m_s', ...
    'groundwater_river_exchange_m_s'},folders.groundwater,'Groundwater_River_Exchange', ...
    'mm_day',1000*86400,true);
add_face('groundwater_seepage',{'groundwater_seepage_face_m_s','groundwater_seepage_m_s'}, ...
    folders.groundwater,'Groundwater_Seepage','mm_day',1000*86400,false);

channel_depth=read_any(native_file,{'channel_depth_node_m','channel_depth_m'});
if ~isempty(channel_depth) && mesh.channel.n_nodes>0
    face_values=zeros(mesh.n_cells,size(channel_depth,2));
    face_values(mesh.channel.host_cell,:)=channel_depth;
    [histories.channel_depth,metadata.channel_depth]=field(map_history(mapping,face_values), ...
        folders.channel,'Channel_Depth','m',false);
end
channel_q=read_any(native_file,{'channel_discharge_link_m3_s'});
if ~isempty(channel_q) && mesh.channel.n_links>0
    face_values=zeros(mesh.n_cells,size(channel_q,2));
    for link=1:mesh.channel.n_links
        host=mesh.channel.host_cell(mesh.channel.link_down(link));
        face_values(host,:)=face_values(host,:)+channel_q(link,:);
    end
    [histories.channel_discharge,metadata.channel_discharge]=field( ...
        map_history(mapping,face_values),folders.channel,'Channel_Discharge','m3_s',true);
end

manifest=table(strings(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{'Variable','Time_s','Path'}); raster_files=strings(0,1);
names=fieldnames(histories);
for k=1:numel(names)
    name=names{k}; maps=histories.(name); meta=metadata.(name);
    if options.write_final_geotiffs
        raster_files(end+1)=write_tiff(maps(:,:,end),meta.folder, ...
            ['Final_' meta.stem '.tif'],mapping,mesh); %#ok<AGROW>
        if meta.absolute_max, maximum=max(abs(maps),[],3); else, maximum=max(maps,[],3); end
        raster_files(end+1)=write_tiff(maximum,meta.folder, ...
            ['Maximum_' meta.stem '.tif'],mapping,mesh); %#ok<AGROW>
    end
    if options.write_temporal_geotiffs
        for j=1:numel(map_indices)
            label=sprintf('%010.0fs',map_time(j)); filename=[meta.stem '_' label '.tif'];
            path=write_tiff(maps(:,:,map_indices(j)),meta.folder,filename,mapping,mesh);
            raster_files(end+1)=path; %#ok<AGROW>
            manifest(end+1,:)={string(name),map_time(j),relative_path(path,output_root)}; %#ok<AGROW>
        end
    end
end
writetable(manifest,fullfile(folders.manifests,'temporal_rasters.csv'));

table_files=write_tables(native_file,time,folders.tables,options);
figure_files=strings(0,1); video_files=strings(0,1);
if options.write_figures, figure_files=write_figures(histories,metadata,time,native_file,folders); end
if options.write_videos, video_files=write_videos(histories,metadata,map_indices,map_time,folders.videos,options.video_fps); end
entries=dir(fullfile(output_root,'**','*')); entries=entries(~[entries.isdir]);
workspace=whos('histories','mapping','surface_depth','velocity','wse');
performance=table(toc(postprocess_clock),numel(entries),sum([entries.bytes]),sum([workspace.bytes]), ...
    'VariableNames',{'Postprocessing_Runtime_s','Output_File_Count','Output_Disk_Bytes', ...
    'Estimated_Working_Memory_Bytes'});
performance_file=fullfile(folders.tables,'Postprocessing_Performance.csv');
writetable(performance,performance_file); table_files(end+1)=performance_file;
report=struct('native_file',string(native_file),'raster_files',raster_files, ...
    'table_files',table_files,'figure_files',figure_files,'video_files',video_files, ...
    'native_time_count',numel(time),'raster_time_count',numel(map_indices), ...
    'output_root',string(output_root),'postprocessing_runtime_s',performance.Postprocessing_Runtime_s, ...
    'output_file_count',performance.Output_File_Count,'output_disk_bytes',performance.Output_Disk_Bytes, ...
    'estimated_working_memory_bytes',performance.Estimated_Working_Memory_Bytes);

    function add_face(name,candidates,folder_name,stem,units,scale,absolute_max)
        values=read_any(native_file,candidates); if isempty(values), return; end
        [histories.(name),metadata.(name)]=field(map_history(mapping,values.*scale), ...
            folder_name,stem,units,absolute_max);
    end
end

function options=defaults(options)
values=struct('write_final_geotiffs',true,'write_temporal_geotiffs',true, ...
    'raster_stack_interval_s',[],'write_figures',true,'write_videos',true, ...
    'write_gauge_hydrographs',true,'video_fps',8);
names=fieldnames(values); for k=1:numel(names), if ~isfield(options,names{k}), options.(names{k})=values.(names{k}); end, end
end

function validate_inputs(native_file,mesh_file,overlap_file,options)
for path={native_file,mesh_file,overlap_file}
    assert(exist(path{1},'file')==2,'HydroPol2D:MissingVoronoiOutput','Missing required file: %s',path{1});
end
if options.write_final_geotiffs || options.write_temporal_geotiffs
    assert(exist('geotiffwrite','file')==2 && exist('maprefcells','file')==2, ...
        'HydroPol2D:MissingMappingToolbox','Voronoi raster output requires Mapping Toolbox.');
end
info=ncinfo(native_file); attributes=string({info.Attributes.Name});
if any(attributes=="output_schema")
    schema=string(ncreadatt(native_file,'/','output_schema'));
    assert(schema=="hydropol2d-unstructured-output-1.0",'HydroPol2D:InvalidVoronoiOutput', ...
        'Unsupported native output schema: %s',schema);
end
variables=string({info.Variables.Name});
assert(any(variables=="time_s") && any(variables=="surface_depth_face_m" | variables=="surface_depth_m"), ...
    'HydroPol2D:InvalidVoronoiOutput','Native archive lacks time or surface-depth histories.');
end

function folders=create_folders(root)
names={'Native','Tables_CSV','Figures_PNG','Figures_PDF','Figures_SVG','GIFs_MP4', ...
    'Rasters_Water_Depths','Rasters_Velocity','Rasters_WSE','Rasters_Hazard', ...
    'Rasters_Rainfall','Rasters_Infiltration','Rasters_ET','Rasters_Recharge', ...
    'Rasters_Groundwater','Rasters_Channel','Manifests'};
keys={'native','tables','figures_png','figures_pdf','figures_svg','videos', ...
    'water_depth','velocity','wse','hazard','rainfall','infiltration','et', ...
    'recharge','groundwater','channel','manifests'};
for k=1:numel(names), folders.(keys{k})=fullfile(root,names{k}); if exist(folders.(keys{k}),'dir')~=7, mkdir(folders.(keys{k})); end, end
end

function values=read_required(path,name)
values=double(ncread(path,name));
end

function values=read_optional(path,name,default_value)
values=read_any(path,{name}); if isempty(values), values=default_value; end
end

function values=read_face_history(path,names,n_faces)
values=read_any(path,names); assert(~isempty(values) && size(values,1)==n_faces, ...
    'HydroPol2D:InvalidVoronoiOutput','Native face history has inconsistent dimensions.');
end

function values=read_any(path,names)
info=ncinfo(path); available=string({info.Variables.Name}); values=[];
for k=1:numel(names), if any(available==string(names{k})), values=double(ncread(path,names{k})); return; end, end
end

function indices=nearest_indices(native_times,map_times)
indices=zeros(numel(map_times),1);
for k=1:numel(map_times), [difference,indices(k)]=min(abs(native_times-map_times(k))); assert(difference<=1e-8*max(abs(map_times(k)),1)); end
indices=unique(indices); if indices(end)~=numel(native_times), indices(end+1)=numel(native_times); end
end

function maps=map_history(mapping,face_history)
flat=mapping.mesh_to_raster*double(face_history); nt=size(flat,2);
maps=zeros(mapping.raster_rows,mapping.raster_cols,nt);
for k=1:nt, maps(:,:,k)=flipud(reshape(flat(:,k),[mapping.raster_cols mapping.raster_rows])'); end
end

function [maps,meta]=field(maps,folder,stem,units,absolute_max)
meta=struct('folder',folder,'stem',stem,'units',units,'absolute_max',absolute_max);
end

function path=write_tiff(values,folder,filename,mapping,mesh)
R=maprefcells([mapping.x_edges(1) mapping.x_edges(end)], ...
    [mapping.y_edges(1) mapping.y_edges(end)],[mapping.raster_rows mapping.raster_cols], ...
    'ColumnsStartFrom','north'); path=string(fullfile(folder,filename)); epsg=[];
if ~isempty(mesh.crs_wkt) && exist('projcrs','file')==2
    try
        crs=projcrs(mesh.crs_wkt);
        if ~isempty(crs.AuthorityCode) && crs.AuthorityCode>0
            epsg=crs.AuthorityCode;
        end
    catch
    end
end
if ~isempty(epsg)
    geotiffwrite(char(path),single(values),R,'CoordRefSysCode',epsg);
else
    write_affine_tiff(char(path),single(values),mapping);
end
end

function write_affine_tiff(path,values,mapping)
% GeoTIFF affine tags without a false CRS for synthetic/local validation grids.
tiff=Tiff(path,'w'); tags=struct(); tags.ImageLength=size(values,1); tags.ImageWidth=size(values,2);
tags.Photometric=Tiff.Photometric.MinIsBlack; tags.BitsPerSample=32;
tags.SamplesPerPixel=1; tags.SampleFormat=Tiff.SampleFormat.IEEEFP;
tags.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky; tags.Compression=Tiff.Compression.LZW;
tags.RowsPerStrip=min(size(values,1),256);
tags.ModelPixelScaleTag=[abs(diff(mapping.x_edges(1:2))) abs(diff(mapping.y_edges(1:2))) 0];
tags.ModelTiepointTag=[0 0 0 mapping.x_edges(1) mapping.y_edges(end) 0];
tiff.setTag(tags); tiff.write(values); tiff.close();
end

function text=relative_path(path,root)
text=erase(string(path),string(root)+filesep);
end

function files=write_tables(path,time,folder,options)
files=strings(0,1); outlet=read_any(path,{'outlet_discharge_m3_s'});
if ~isempty(outlet), file=fullfile(folder,'Outlet_Hydrograph.csv'); writetable(table(time,outlet(:),'VariableNames',{'Time_s','Discharge_m3_s'}),file); files(end+1)=file; end
diagnostic_time=read_any(path,{'diagnostic_time_s'});
if ~isempty(diagnostic_time)
    info=ncinfo(path); names=string({info.Variables.Name}); columns={'dt_s','mass_m3','step_mass_residual_m3', ...
        'max_surface_depth_m','max_surface_velocity_m_s','max_channel_depth_m', ...
        'max_channel_velocity_m_s','boundary_net_inflow_volume_m3','precipitation_volume_m3', ...
        'infiltration_volume_m3','actual_et_volume_m3','recharge_volume_m3'};
    T=table(diagnostic_time(:),'VariableNames',{'Time_s'});
    for k=1:numel(columns), if any(names==columns{k}), T.(columns{k})=double(ncread(path,columns{k})); end, end
    file=fullfile(folder,'Timestep_Diagnostics.csv'); writetable(T,file); files(end+1)=file;
    balance=table(diagnostic_time(:),'VariableNames',{'Time_s'});
    cumulative={'precipitation_volume_m3','boundary_net_inflow_volume_m3','infiltration_volume_m3', ...
        'actual_et_volume_m3','recharge_volume_m3'};
    for k=1:numel(cumulative), if any(names==cumulative{k}), balance.(['Cumulative_' cumulative{k}])=cumsum(double(ncread(path,cumulative{k}))); end, end
    if any(names=="step_mass_residual_m3"), balance.Water_Balance_Residual_m3=double(ncread(path,'step_mass_residual_m3')); end
    file=fullfile(folder,'Water_Balance.csv'); writetable(balance,file); files(end+1)=file;
end
if options.write_gauge_hydrographs
    gauge_q=read_any(path,{'gauge_discharge_m3_s'});
    if ~isempty(gauge_q)
        T=array2table(gauge_q','VariableNames',gauge_names(path,size(gauge_q,1))); T=addvars(T,time,'Before',1,'NewVariableNames','Time_s');
        file=fullfile(folder,'Gauge_Hydrographs.csv'); writetable(T,file); files(end+1)=file;
    end
end
end

function names=gauge_names(path,count)
try
    names=split(string(ncreadatt(path,'/','gauge_names')),',')';
catch
    names="Gauge_"+(1:count);
end
names=matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(cellstr(names)));
end

function files=write_figures(histories,metadata,time,native_file,folders)
files=strings(0,1); names=fieldnames(histories); preferred={'surface_depth','velocity','infiltration_depth','actual_et_depth','groundwater_head','groundwater_depth'};
selected=intersect(preferred,names,'stable'); if isempty(selected), return; end
figure_handle=figure('Visible','off','Color','w','Position',[50 50 1200 800]); layout=tiledlayout(figure_handle,2,3,'TileSpacing','compact','Padding','compact');
for k=1:numel(selected), name=selected{k}; ax=nexttile(layout); imagesc(ax,histories.(name)(:,:,end)); axis(ax,'image','off'); colormap(ax,parula(256)); colorbar(ax); title(ax,strrep(metadata.(name).stem,'_',' '),'FontName','Helvetica','FontWeight','normal'); set(ax,'FontName','Helvetica','LineWidth',1.5); end
files=[files; save_figure(figure_handle,folders,'State_Summary')]; close(figure_handle);
outlet=read_any(native_file,{'outlet_discharge_m3_s'});
if ~isempty(outlet)
    figure_handle=figure('Visible','off','Color','w','Position',[50 50 720 420]); plot(time/3600,outlet,'LineWidth',1.5,'Color',[0 0.447 0.741]); grid on; xlabel('Time (h)'); ylabel('Discharge (m^3/s)'); set(gca,'FontName','Helvetica','LineWidth',1.5); title('Outlet hydrograph','FontWeight','normal');
    files=[files;save_figure(figure_handle,folders,'Outlet_Hydrograph')]; close(figure_handle);
end
gauge_q=read_any(native_file,{'gauge_discharge_m3_s'});
if ~isempty(gauge_q)
    figure_handle=figure('Visible','off','Color','w','Position',[50 50 720 420]);
    plot(time/3600,gauge_q','LineWidth',1.5); grid on; xlabel('Time (h)'); ylabel('Discharge (m^3/s)');
    set(gca,'FontName','Helvetica','LineWidth',1.5); title('Configured-gauge hydrographs','FontWeight','normal');
    legend(gauge_names(native_file,size(gauge_q,1)),'Location','best','Interpreter','none');
    files=[files;save_figure(figure_handle,folders,'Gauge_Hydrographs')]; close(figure_handle);
end
end

function files=save_figure(handle,folders,stem)
files=[string(fullfile(folders.figures_png,[stem '.png']));string(fullfile(folders.figures_pdf,[stem '.pdf']));string(fullfile(folders.figures_svg,[stem '.svg']))];
exportgraphics(handle,files(1),'Resolution',300); exportgraphics(handle,files(2),'ContentType','vector'); exportgraphics(handle,files(3),'ContentType','vector');
end

function files=write_videos(histories,metadata,indices,map_time,folder,fps)
files=strings(0,1); names=fieldnames(histories);
for k=1:numel(names), name=names{k}; maps=histories.(name); if size(maps,3)<2, continue; end
    file=fullfile(folder,[metadata.(name).stem '.mp4']); video=VideoWriter(file,'MPEG-4'); video.FrameRate=fps; open(video);
    limits=[min(maps(:)) max(maps(:))]; if limits(2)<=limits(1), limits(2)=limits(1)+1; end
    figure_handle=figure('Visible','off','Color','w','Position',[50 50 720 520]);
    for j=1:numel(indices), imagesc(maps(:,:,indices(j)),limits); axis image off; colormap(parula(256)); colorbar; title(sprintf('%s | t = %.2f h',strrep(metadata.(name).stem,'_',' '),map_time(j)/3600),'FontName','Helvetica','FontWeight','normal'); drawnow; writeVideo(video,getframe(figure_handle)); end
    close(video); close(figure_handle); files(end+1)=file; %#ok<AGROW>
end
preferred={'surface_depth','velocity','infiltration_depth','actual_et_depth', ...
    'groundwater_head','groundwater_depth'};
selected=intersect(preferred,names,'stable');
if numel(selected)>=2
    file=fullfile(folder,'Standard_Multipanel.mp4'); video=VideoWriter(file,'MPEG-4');
    video.FrameRate=fps; open(video); figure_handle=figure('Visible','off','Color','w','Position',[50 50 1200 760]);
    for j=1:numel(indices)
        clf(figure_handle); layout=tiledlayout(figure_handle,2,3,'TileSpacing','compact','Padding','compact');
        for k=1:numel(selected)
            name=selected{k}; maps=histories.(name); ax=nexttile(layout);
            imagesc(ax,maps(:,:,indices(j))); axis(ax,'image','off'); colormap(ax,parula(256)); colorbar(ax);
            title(ax,strrep(metadata.(name).stem,'_',' '),'FontName','Helvetica','FontWeight','normal');
        end
        title(layout,sprintf('HydroPol2D | t = %.2f h',map_time(j)/3600), ...
            'FontName','Helvetica','FontWeight','normal'); drawnow; writeVideo(video,getframe(figure_handle));
    end
    close(video); close(figure_handle); files(end+1)=file;
end
end
