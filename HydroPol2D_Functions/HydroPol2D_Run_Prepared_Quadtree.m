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
groundwater_head=map_mean_required(fullfile(static_root,'GW_table.tif'),mapping,'initial groundwater head');
[groundwater_fixed_mask,groundwater_fixed_head]=map_fixed_head( ...
    fullfile(static_root,'GW_Dirichlet_Head.tif'),mapping,groundwater_head);

land_cover_parameters=config_json.class_parameters.land_cover;
soil_parameters=config_json.class_parameters.soil;
roughness=parameter_lookup(land_cover_parameters,lulc,'manning_n', ...
    double(config_json.parameters.surface_manning_fallback_n));
root_depth=parameter_lookup_required(land_cover_parameters,lulc,'root_depth_m');
impervious=logical(parameter_lookup_required(land_cover_parameters,lulc,'impervious'));
tr=parameter_lookup_required(soil_parameters,soil,'theta_r');
ts=parameter_lookup_required(soil_parameters,soil,'theta_s');
alpha_vg=parameter_lookup_required(soil_parameters,soil,'alpha_vg_m_inverse');
n_vg=parameter_lookup_required(soil_parameters,soil,'n_vg');
ksat_mm_h=parameter_lookup_required(soil_parameters,soil,'ksat_mm_h');
specific_yield=parameter_lookup_required(soil_parameters,soil,'specific_yield');
groundwater_ksat_mm_h=parameter_lookup_required(soil_parameters,soil,'groundwater_ksat_mm_h');

flags=config_json.flags;
hydrology=struct( ...
    'pervious_fraction',double(~impervious), ...
    'lai',lai.*double(flags.canopy_interception), ...
    'soil_depth_m',dtb,'root_depth_m',min(dtb,root_depth), ...
    'theta_r',tr,'theta_sat',ts, ...
    'alpha_vg_1_m',alpha_vg,'n_vg',n_vg, ...
    'ksat_m_s',ksat_mm_h./1000./3600, ...
    'surface_layer_depth_m',double(config_json.parameters.surface_layer_depth_m), ...
    'infiltration_head_cap_m',double(config_json.parameters.infiltration_head_cap_m), ...
    'surface_conductivity_weight',double(config_json.parameters.surface_conductivity_weight), ...
    'wilting_saturation',double(config_json.parameters.wilting_saturation), ...
    'capillary_field_capacity_head_m',double(config_json.parameters.capillary_field_capacity_head_m), ...
    'initial_soil_saturation',min(max((initial_theta-tr)./max(ts-tr,eps),0),1));
if ~flags.infiltration, hydrology.pervious_fraction(:)=0; end

simulation=config_json.simulation;
outputs=struct('profile','standard','map_interval_minutes',60, ...
    'diagnostic_interval_minutes',10,'map_variables',{{'surface_depth','surface_velocity'}});
if isfield(config_json,'outputs')
    names=fieldnames(config_json.outputs);
    for k=1:numel(names), outputs.(names{k})=config_json.outputs.(names{k}); end
end
map_interval_s=double(outputs.map_interval_minutes)*60;
diagnostic_interval_s=double(outputs.diagnostic_interval_minutes)*60;
profile=lower(string(outputs.profile));
write_full=profile=="full";
write_native=profile=="standard" || write_full;
write_media=profile=="standard" || write_full;
groundwater_relative_tolerance=1e-3;
groundwater_head_tolerance_m=1e-3;
if isfield(config_json.parameters,'groundwater_relative_tolerance')
    groundwater_relative_tolerance=double(config_json.parameters.groundwater_relative_tolerance);
end
if isfield(config_json.parameters,'groundwater_head_tolerance_m')
    groundwater_head_tolerance_m=double(config_json.parameters.groundwater_head_tolerance_m);
end
solver_config=struct( ...
    'duration_s',double(simulation.duration_seconds), ...
    'output_interval_s',map_interval_s,'raster_output_interval_s',map_interval_s, ...
    'routing_solver','local_inertial','surface_roughness',roughness, ...
    'courant',double(config_json.parameters.routing_courant), ...
    'min_dt_s',double(simulation.minimum_timestep_s), ...
    'timestep_floor_action','clamp', ...
    'max_dt_s',double(simulation.maximum_timestep_s),'critical_flow',true, ...
    'hydrology_enabled',logical(flags.canopy_interception || flags.infiltration || flags.internal_etp), ...
    'hydrology',hydrology,'groundwater_enabled',logical(flags.groundwater), ...
    'river_groundwater_exchange',logical(flags.river_groundwater_exchange), ...
    'riverbed_thickness_m',double(config_json.parameters.riverbed_thickness_m), ...
    'groundwater_lateral_flow',~strcmpi(string(config_json.parameters.groundwater_boundary),'local'), ...
    'initial_groundwater_head_m',min(max(groundwater_head,mesh.cell_bed-dtb),mesh.cell_bed), ...
    'aquifer_bottom_m',mesh.cell_bed-dtb, ...
    'hydraulic_conductivity_m_s',groundwater_ksat_mm_h./1000./3600, ...
    'specific_yield',specific_yield, ...
    'groundwater_relative_tolerance',groundwater_relative_tolerance, ...
    'groundwater_head_tolerance_m',groundwater_head_tolerance_m, ...
    'groundwater_courant',double(config_json.parameters.groundwater_courant), ...
    'groundwater_fixed_head_mask',groundwater_fixed_mask, ...
    'groundwater_fixed_head_m',groundwater_fixed_head, ...
    'map_variables',{cellstr(string(outputs.map_variables))}, ...
    'forcing_interval_s',double(config_json.forcing.rainfall_interval_minutes)*60, ...
    'progress_interval_s',diagnostic_interval_s, ...
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
        meteorology=read_etp_workbook(etp_path,mesh, ...
            double(config_json.parameters.reference_et_krs), ...
            double(config_json.parameters.reference_crop_albedo));
        forcing.meteorology=@(time_s,state,current_mesh) meteorology_at(time_s,meteorology);
    else
        forcing.potential_et_m_s=0;
        warning('HydroPol2D:QuadtreeETFallback', ...
            'The prepared ETP workbook is absent; evapotranspiration is zero.');
    end
end
forcing.surface_boundary=prepared_surface_boundary(case_root,mesh,simulation.start_utc);

case_definition=struct( ...
    'mesh_file',mesh_file,'overlap_file',overlap_file, ...
    'fine_resolution_m',double(simulation.resolution_m), ...
    'coarse_resolution_m',double(simulation.coarse_resolution_m), ...
    'urban_refinement_buffer_m',double(simulation.urban_refinement_buffer_m), ...
    'config',solver_config,'forcing',forcing,'output_prefix','quadtree', ...
    'output_controls',struct('write_native_archive',write_native, ...
        'native_output_interval_s',map_interval_s,'write_final_geotiffs',true, ...
        'write_temporal_geotiffs',write_full,'raster_stack_interval_s',map_interval_s, ...
        'write_figures',write_media,'write_videos',write_media, ...
        'write_gauge_hydrographs',true, ...
        'map_variables',{cellstr(string(outputs.map_variables))}));
if exist(subgrid_path,'file')==2
    case_definition.options=struct('voronoi_subgrid_enabled',true, ...
        'subgrid_table_path',subgrid_path);
end

hp2d_write_run_monitor(run_progress_path,run_metrics_path, ...
    struct('stage','simulation','percent',0),false);
results=HydroPol2D_Run_Quadtree_Case(case_definition,fullfile(case_root,'Outputs'),true);
HydroPol2D_Write_Run_Manifest(case_root,results,fullfile(case_root,'Outputs','Modeling_Results'));
hp2d_write_run_monitor(run_progress_path,run_metrics_path, ...
    struct('stage','complete','percent',100),false);
end

function forcing=prepared_surface_boundary(case_root,mesh,start_utc)
plan_path=fullfile(case_root,'processed_inputs','Forcing','Boundary','surface_boundary_plan.json');
series_path=fullfile(case_root,'processed_inputs','Forcing','Boundary','surface_boundary.csv');
if exist(plan_path,'file')~=2
    boundary=find(mesh.edge_neighbor==0);
    forcing=struct('edge_id',boundary,'type',"critical_flow",'value',0);
    return
end
plan=jsondecode(fileread(plan_path));
critical=double(plan.critical_flow_edge_ids(:))+1;
if exist(series_path,'file')~=2
    forcing=struct('edge_id',critical,'type',"critical_flow",'value',0);
    return
end
table=readtable(series_path,'TextType','string');
if isempty(table)
    forcing=struct('edge_id',critical,'type',"critical_flow",'value',0);
    return
end
origin=datetime(start_utc,'TimeZone','UTC','InputFormat','yyyy-MM-dd''T''HH:mm:ssXXX');
timestamps=datetime(table.time_utc,'TimeZone','UTC','InputFormat','yyyy-MM-dd''T''HH:mm:ssXXX');
data=struct('time_s',seconds(timestamps-origin),'edge_id',double(table.edge_id)+1, ...
    'type',string(table.type),'value',double(table.value),'critical',critical, ...
    'glofas_interpolation',string(plan.inflows_interpolation),'tide_interpolation',string(plan.tide.interpolation));
forcing=@(time_s,state,current_mesh) boundary_at_time(data,time_s);
end

function boundary=boundary_at_time(data,time_s)
critical=data.critical(:);
dynamic=unique(data.edge_id(:));
n_critical=numel(critical);
edge_id=zeros(n_critical+numel(dynamic),1);
types=strings(size(edge_id));
values=zeros(size(edge_id));
edge_id(1:n_critical)=critical;
types(1:n_critical)="critical_flow";
for k=1:numel(dynamic)
    selected=data.edge_id==dynamic(k);
    times=data.time_s(selected); series=data.value(selected); kind=data.type(find(selected,1));
    [times,order]=sort(times); series=series(order);
    method=data.glofas_interpolation;
    if kind=="stage", method=data.tide_interpolation; end
    if method=="linear"
        value=interp1(times,series,time_s,'linear','extrap');
    else
        index=find(times<=time_s,1,'last');
        if isempty(index), index=1; end
        value=series(index);
    end
    output_index=n_critical+k;
    edge_id(output_index)=dynamic(k);
    types(output_index)=kind;
    values(output_index)=value;
end
boundary=struct('edge_id',edge_id,'type',types,'value',values);
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

function value=map_mean_required(path,mapping,label)
assert(isfile(path),'HydroPol2D:MissingPreparedRaster', ...
    'The current solver requires the prepared %s raster: %s.',label,path);
raw=read_south_up(path); valid=isfinite(raw);
raw(~valid)=0;
numerator=mapping.raster_to_mesh*raw;
denominator=mapping.raster_to_mesh*double(valid);
assert(all(denominator>0),'HydroPol2D:IncompletePreparedRasterCoverage', ...
    ['The prepared %s raster does not cover every native mesh cell. ' ...
     'Rebuild the case instead of using an implicit groundwater fallback.'],label);
value=numerator./denominator;
assert(all(isfinite(value)),'HydroPol2D:InvalidPreparedRaster', ...
    'The prepared %s raster produced non-finite native values.',label);
end

function [mask,value]=map_fixed_head(path,mapping,fallback)
if exist(path,'file')~=2
    mask=false(size(fallback)); value=fallback; return
end
raw=read_south_up(path); valid=isfinite(raw); raw(~valid)=0;
weight=mapping.raster_to_mesh*double(valid);
value=(mapping.raster_to_mesh*raw)./max(weight,eps);
mask=weight>0; value(~mask)=fallback(~mask);
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
% Overlap raster indices are C/NetCDF row-major. MATLAB's ``(:)`` is
% column-major, so transpose before vectorising to preserve the shared
% south-up row-major contract used by the Builder and Python solver.
raw=flipud(array).'; raw=raw(:);
end

function value=parameter_lookup(records,code,field_name,fallback)
value=repmat(double(fallback),size(code));
for k=1:numel(records)
    selected=code==double(records(k).index);
    raw=records(k).(field_name);
    value(selected)=double(raw);
end
end

function value=parameter_lookup_required(records,code,field_name)
value=nan(size(code));
for k=1:numel(records)
    assert(isfield(records(k),field_name), ...
        'HydroPol2D:MissingPhysicsV2Parameter', ...
        'HydroPol2D class %g is missing %s.',double(records(k).index),field_name);
    selected=code==double(records(k).index);
    value(selected)=double(records(k).(field_name));
end
assert(all(isfinite(value)), ...
    'HydroPol2D:MissingPhysicsV2Class', ...
    'HydroPol2D has no %s value for one or more mapped classes.',field_name);
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

function data=read_etp_workbook(path,mesh,reference_et_krs,reference_crop_albedo)
raw=readcell(path,'Sheet','ETP');
station_count=floor((size(raw,2)-2)/6);
assert(station_count>0,'HydroPol2D:InvalidETP','The ETP workbook has no ERA5 stations.');
x=zeros(station_count,1); y=zeros(station_count,1);
for station=1:station_count
    first=3+6*(station-1);
    x(station)=double(raw{1,first+3}); y(station)=double(raw{1,first+5});
end
projected=projcrs(mesh.crs_wkt);
[station_latitude,station_longitude]=projinv(projected,x,y);
[mesh_latitude,mesh_longitude]=projinv(projected,mesh.cell_x,mesh.cell_y);
[longitude,~,longitude_index]=unique(round(station_longitude,10));
[latitude,~,latitude_index]=unique(round(station_latitude,10));
use_linear=numel(longitude)>=2 && numel(latitude)>=2 && ...
    numel(longitude)*numel(latitude)==station_count;
if use_linear
    target_longitude=min(max(mesh_longitude,min(longitude)),max(longitude));
    target_latitude=min(max(mesh_latitude,min(latitude)),max(latitude));
else
    nearest=zeros(mesh.n_cells,1);
    for cell_id=1:mesh.n_cells
        [~,nearest(cell_id)]=min((x-mesh.cell_x(cell_id)).^2+(y-mesh.cell_y(cell_id)).^2);
    end
end
rows=3:size(raw,1); count=numel(rows);
fields={'maximum_temperature_c','minimum_temperature_c','temperature_c', ...
    'wind_speed_m_s','relative_humidity_pct','ground_heat_flux_mj_m2_day'};
for name=1:numel(fields), data.(fields{name})=zeros(count,mesh.n_cells); end
data.times=NaT(count,1,'TimeZone','UTC');
data.reference_et_krs=reference_et_krs;
data.reference_crop_albedo=reference_crop_albedo;
for day_index=1:count
    timestamp=raw{rows(day_index),2};
    if isdatetime(timestamp), data.times(day_index)=datetime(timestamp,'TimeZone','UTC');
    else, data.times(day_index)=datetime(timestamp,'ConvertFrom','excel','TimeZone','UTC'); end
    for name=1:numel(fields)
        station_values=zeros(station_count,1);
        for station=1:station_count
            first=3+6*(station-1);
            station_values(station)=double(raw{rows(day_index),first+name-1});
        end
        if use_linear
            grid=nan(numel(latitude),numel(longitude));
            grid(sub2ind(size(grid),latitude_index,longitude_index))=station_values;
            data.(fields{name})(day_index,:)=interp2(longitude,latitude,grid, ...
                target_longitude,target_latitude,'linear')';
        else
            data.(fields{name})(day_index,:)=station_values(nearest)';
        end
    end
end
data.latitude_deg=mesh_latitude;
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
    'albedo',data.reference_crop_albedo,'krs',data.reference_et_krs);
end
