function results = HydroPol2D_Run_Prepared_Quadtree(case_root)
%HYDROPOL2D_RUN_PREPARED_QUADTREE Run a Case Builder quadtree package.

arguments
    case_root (1,:) char
end
config_json=jsondecode(fileread(fullfile(case_root,'hydropol_config.json')));
run_progress_path=strtrim(getenv('HYDROPOL_PROGRESS_FILE'));
run_metrics_path=strtrim(getenv('HYDROPOL_METRICS_FILE'));
hp2d_write_run_monitor(run_progress_path,run_metrics_path, ...
    struct('stage','preprocessing','percent',NaN),false);
static_root=fullfile(case_root,'processed_inputs','Static');
forcing_root=fullfile(case_root,'processed_inputs','Forcing');
mesh_file=fullfile(case_root,'processed_inputs','Mesh','hydropol_quadtree_mesh.nc');
overlap_file=fullfile(case_root,'processed_inputs','Mesh','hydropol_quadtree_overlap.nc');
mesh=HydroPol2D_Read_UGRID(mesh_file);
mapping=HydroPol2D_Read_Overlap(overlap_file,mesh_file);

lulc=map_category(fullfile(static_root,'LULC.tif'),mapping);
soil=max(1,min(12,map_category(fullfile(static_root,'SOIL.tif'),mapping)));
lai=map_mean(fullfile(static_root,'LAI.tif'),mapping,0);
initial_theta=map_mean(fullfile(static_root,'Initial_Soil_Moisture_Fraction.tif'),mapping,0.25);
dtb=max(map_mean(fullfile(static_root,'DTB.tif'),mapping,2),0.1);
groundwater_head=map_mean(fullfile(static_root,'GW_table.tif'),mapping,mesh.cell_bed-2);

roughness=lookup(lulc,[10 20 30 40 50 60 70 80 90 95 100], ...
    [0.100 0.080 0.060 0.050 0.030 0.035 0.020 0.035 0.120 0.150 0.080],0.05);
theta_r=[0.068 0.089 0.075 0.095 0.089 0.065 0.067 0.078 0.100 0.034 0.049 0.045];
theta_sat=[0.385 0.423 0.321 0.309 0.432 0.330 0.432 0.399 0.387 0.481 0.390 0.430];
alpha_vg=[1.09 1.23 1.31 1.31 1.23 1.48 1.68 1.56 1.89 1.37 1.75 2.68];
n_vg=[1.8 2.0 2.5 2.9 2.5 4.0 3.0 4.6 7.0 2.6 12.0 15.5];
ksat_mm_h=[0.3 0.5 0.6 1.0 1.0 1.5 7.6 3.4 10.9 3.4 29.9 117.8];
soil_index=round(soil);
tr=theta_r(soil_index)'; ts=theta_sat(soil_index)';

flags=config_json.flags;
hydrology=struct( ...
    'pervious_fraction',double(lulc~=50 & lulc~=80), ...
    'lai',lai.*double(flags.canopy_interception), ...
    'soil_depth_m',dtb,'root_depth_m',min(dtb,1), ...
    'theta_r',tr,'theta_sat',ts, ...
    'alpha_vg_1_m',alpha_vg(soil_index)', ...
    'n_vg',n_vg(soil_index)', ...
    'ksat_m_s',ksat_mm_h(soil_index)'./1000./3600, ...
    'initial_soil_saturation',min(max((initial_theta-tr)./max(ts-tr,eps),0),1));
if ~flags.infiltration, hydrology.pervious_fraction(:)=0; end

simulation=config_json.simulation;
solver_config=struct( ...
    'duration_s',double(simulation.duration_seconds), ...
    'output_interval_s',3600,'raster_output_interval_s',3600, ...
    'routing_solver','local_inertial','surface_roughness',roughness, ...
    'courant',0.2,'min_dt_s',double(simulation.minimum_timestep_s), ...
    'max_dt_s',double(simulation.maximum_timestep_s),'critical_flow',true, ...
    'hydrology_enabled',logical(flags.canopy_interception || flags.infiltration || flags.internal_etp), ...
    'hydrology',hydrology,'groundwater_enabled',logical(flags.groundwater), ...
    'initial_groundwater_head_m',min(max(groundwater_head,mesh.cell_bed-dtb),mesh.cell_bed), ...
    'aquifer_bottom_m',mesh.cell_bed-dtb,'hydraulic_conductivity_m_s',1e-5, ...
    'specific_yield',0.2,'groundwater_update_interval_s',3600, ...
    'forcing_interval_s',double(config_json.forcing.rainfall_interval_minutes)*60, ...
    'progress_interval_s',max(double(simulation.duration_seconds)/100,60), ...
    'progress_checkpoint_file',run_progress_path);
subgrid_path=fullfile(case_root,'processed_inputs','Mesh','hydropol_quadtree_subgrid.nc');
if exist(subgrid_path,'file')==2
    solver_config.voronoi_subgrid_enabled=true;
    solver_config.subgrid_table_path=subgrid_path;
end

rainfall_files=dir(fullfile(forcing_root,'Rainfall','*.tif'));
[~,order]=sort({rainfall_files.name}); rainfall_files=rainfall_files(order);
assert(~isempty(rainfall_files),'HydroPol2D:MissingRainfall','No prepared rainfall maps were found.');
projected=projcrs(mesh.crs_wkt);
[mesh_latitude,mesh_longitude]=projinv(projected,mesh.cell_x,mesh.cell_y);
forcing=struct();
forcing.surface_source_m_s=@(time_s,state,current_mesh) rainfall_at( ...
    time_s,rainfall_files,forcing_root,double(config_json.forcing.rainfall_interval_minutes)*60, ...
    mesh_latitude,mesh_longitude);
if flags.internal_etp
    etp_path=fullfile(forcing_root,'Evapotranspiration','ETP_input_data.xlsx');
    if exist(etp_path,'file')==2
        meteorology=read_etp_workbook(etp_path,mesh);
        forcing.meteorology=@(time_s,state,current_mesh) meteorology_at(time_s,meteorology);
    else
        forcing.potential_et_m_s=0;
        warning('HydroPol2D:QuadtreeETFallback', ...
            'The prepared ETP workbook is absent; evapotranspiration is zero.');
    end
end
boundary=find(mesh.edge_neighbor==0);
forcing.surface_boundary=struct('edge_id',boundary,'type',"critical_flow",'value',0);

case_definition=struct( ...
    'mesh_file',mesh_file,'overlap_file',overlap_file, ...
    'fine_resolution_m',double(simulation.resolution_m), ...
    'coarse_resolution_m',double(simulation.coarse_resolution_m), ...
    'urban_refinement_buffer_m',double(simulation.urban_refinement_buffer_m), ...
    'config',solver_config,'forcing',forcing,'output_prefix','quadtree', ...
    'output_controls',struct('write_native_archive',true, ...
        'native_output_interval_s',3600,'write_final_geotiffs',false, ...
        'write_temporal_geotiffs',false,'write_figures',false, ...
        'write_videos',false,'write_gauge_hydrographs',false));
if exist(subgrid_path,'file')==2
    case_definition.options=struct('voronoi_subgrid_enabled',true, ...
        'subgrid_table_path',subgrid_path);
end
hp2d_write_run_monitor(run_progress_path,run_metrics_path, ...
    struct('stage','simulation','percent',0),false);
results=HydroPol2D_Run_Quadtree_Case(case_definition,fullfile(case_root,'Outputs'),false);
hp2d_write_run_monitor(run_progress_path,run_metrics_path, ...
    struct('stage','complete','percent',100),false);
end

function value=map_mean(path,mapping,fallback)
raw=read_south_up(path); valid=isfinite(raw);
raw(~valid)=0;
numerator=mapping.raster_to_mesh*raw;
denominator=mapping.raster_to_mesh*double(valid);
value=numerator./max(denominator,eps);
if isscalar(fallback), fallback=repmat(fallback,size(value)); end
value(denominator<=0)=fallback(denominator<=0);
end

function value=map_category(path,mapping)
raw=read_south_up(path); classes=unique(raw(isfinite(raw)));
score=zeros(size(mapping.raster_to_mesh,1),numel(classes));
for k=1:numel(classes)
    score(:,k)=mapping.raster_to_mesh*double(raw==classes(k));
end
[~,selected]=max(score,[],2); value=classes(selected);
end

function raw=read_south_up(path)
[array,~]=readgeoraster(path,'OutputType','double');
info=georasterinfo(path);
for missing=double(info.MissingDataIndicator(:))'
    array(array==missing)=NaN;
end
raw=flipud(array); raw=raw(:);
end

function value=lookup(code,keys,values,fallback)
value=repmat(fallback,size(code));
for k=1:numel(keys), value(code==keys(k))=values(k); end
end

function rate=rainfall_at(time_s,files,forcing_root,interval_s,latitude,longitude)
index=min(numel(files),max(1,floor(time_s/interval_s)+1));
path=fullfile(forcing_root,'Rainfall',files(index).name);
[raw,reference]=readgeoraster(path,'OutputType','double');
for missing=double(georasterinfo(path).MissingDataIndicator(:))'
    raw(raw==missing)=NaN;
end
[column,row]=geographicToIntrinsic(reference,latitude,longitude);
sampled=interp2(raw,column,row,'nearest',0);
sampled(~isfinite(sampled))=0;
rate=max(sampled(:),0)./1000/3600;
end

function data=read_etp_workbook(path,mesh)
raw=readcell(path,'Sheet','ETP');
station_count=floor((size(raw,2)-2)/6);
assert(station_count>0,'HydroPol2D:InvalidETP','The ETP workbook has no ERA5 stations.');
x=zeros(station_count,1); y=zeros(station_count,1);
for station=1:station_count
    first=3+6*(station-1);
    x(station)=double(raw{1,first+3}); y(station)=double(raw{1,first+5});
end
nearest=zeros(mesh.n_cells,1);
for cell_id=1:mesh.n_cells
    [~,nearest(cell_id)]=min((x-mesh.cell_x(cell_id)).^2+(y-mesh.cell_y(cell_id)).^2);
end
rows=3:size(raw,1); count=numel(rows);
fields={'maximum_temperature_c','minimum_temperature_c','temperature_c', ...
    'wind_speed_m_s','relative_humidity_pct','ground_heat_flux_mj_m2_day'};
for name=1:numel(fields), data.(fields{name})=zeros(count,mesh.n_cells); end
data.times=NaT(count,1,'TimeZone','UTC');
for day_index=1:count
    timestamp=raw{rows(day_index),2};
    if isdatetime(timestamp), data.times(day_index)=datetime(timestamp,'TimeZone','UTC');
    else, data.times(day_index)=datetime(timestamp,'ConvertFrom','excel','TimeZone','UTC'); end
    for cell_id=1:mesh.n_cells
        first=3+6*(nearest(cell_id)-1);
        for name=1:numel(fields)
            data.(fields{name})(day_index,cell_id)=double(raw{rows(day_index),first+name-1});
        end
    end
end
try
    projected=projcrs(mesh.crs_wkt);
    [data.latitude_deg,~]=projinv(projected,mesh.cell_x,mesh.cell_y);
catch
    data.latitude_deg=zeros(mesh.n_cells,1);
    warning('HydroPol2D:QuadtreeLatitudeFallback', ...
        'Could not transform mesh centres to latitude; ET uses 0 degrees latitude.');
end
end

function result=meteorology_at(time_s,data)
index=min(numel(data.times),max(1,floor(time_s/86400)+1));
result=struct( ...
    'temperature_c',data.temperature_c(index,:)', ...
    'maximum_temperature_c',data.maximum_temperature_c(index,:)', ...
    'minimum_temperature_c',data.minimum_temperature_c(index,:)', ...
    'day_of_year',day(data.times(index),'dayofyear'), ...
    'latitude_deg',data.latitude_deg, ...
    'wind_speed_m_s',data.wind_speed_m_s(index,:)', ...
    'relative_humidity_pct',data.relative_humidity_pct(index,:)', ...
    'ground_heat_flux_mj_m2_day',data.ground_heat_flux_mj_m2_day(index,:)', ...
    'albedo',0.25,'krs',0.16);
end
