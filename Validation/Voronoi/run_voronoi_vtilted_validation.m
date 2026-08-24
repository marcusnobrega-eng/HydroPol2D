function Results = run_voronoi_vtilted_validation(mesh_directory, output_directory)
%RUN_VORONOI_VTILTED_VALIDATION Compare D4 and Voronoi V-tilted routing.

arguments
    mesh_directory (1,:) char
    output_directory (1,:) char
end
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(fullfile(root, 'HydroPol2D_Functions'));
if ~exist(output_directory, 'dir'), mkdir(output_directory); end

reference = run_d4_reference();
names = ["fine"; "variable"];
Summary = table();
all_series = cell(numel(names),1); all_maps = cell(numel(names),1);
for k = 1:numel(names)
    name = names(k);
    mesh_file = fullfile(mesh_directory, "vtilted-" + name + "-mesh.nc");
    overlap_file = fullfile(mesh_directory, "vtilted-" + name + "-overlap.nc");
    [series, final_map, metrics] = run_voronoi_case(mesh_file, overlap_file, reference, ...
        fullfile(output_directory, "vtilted-" + name + "-results.nc"));
    writetable(series, fullfile(output_directory, "vtilted-" + name + "-timeseries.csv"));
    Summary = [Summary; struct2table(metrics)]; %#ok<AGROW>
    all_series{k}=series; all_maps{k}=final_map;
    if name == "fine", fine_map = final_map; else, variable_map = final_map; end
end
writetable(reference.series, fullfile(output_directory, 'vtilted-d4-timeseries.csv'));
writetable(Summary, fullfile(output_directory, 'vtilted-comparison-summary.csv'));
Convergence = compare_voronoi_meshes(all_series{1},all_series{2},all_maps{1},all_maps{2}, ...
    Summary.cell_count(1),Summary.cell_count(2));
writetable(Convergence, fullfile(output_directory, 'vtilted-voronoi-convergence-summary.csv'));

[rows, cols] = size(reference.final_depth_m);
[cc, rr] = meshgrid(1:cols, 1:rows);
x_m = (cc - 0.5) * 20;
y_m = (rr - 0.5) * 20;
Maps = table(rr(:),cc(:),x_m(:),y_m(:),reference.final_depth_m(:),fine_map(:),variable_map(:), ...
    'VariableNames',{'row','column','x_m','y_m','d4_depth_m','fine_voronoi_depth_m','variable_voronoi_depth_m'});
writetable(Maps, fullfile(output_directory, 'vtilted-final-depth-maps.csv'));
save(fullfile(output_directory, 'vtilted-comparison.mat'), 'reference', 'fine_map', 'variable_map', 'Summary', 'Convergence');
Results = struct('d4_comparison',Summary,'mesh_convergence',Convergence);
disp(Summary);
disp(Convergence);
end

function Reference = run_d4_reference()
dx = 20; nx = 41; ny = 31; dt = 1; duration = 90*60; record_dt = 60;
side_slope = 0.02; downslope = 0.004; roughness = 0.04; rain = 30e-3/3600;
[cc, rr] = meshgrid(1:nx, 1:ny); center = ceil(nx/2);
z = side_slope * abs(cc-center) * dx + downslope * (ny-rr) * dx;
h = zeros(ny,nx); n = roughness * ones(ny,nx); outflow = zeros(ny,nx,3);
outlet = false(ny,nx); outlet(ny,center-1:center+1) = true;
[outlet_row,outlet_col] = find(outlet);
area = dx^2; records = (0:record_dt:duration)'; count = numel(records);
q = zeros(count,1); storage = zeros(count,1); cumulative_out = zeros(count,1);
mass_error_pct = zeros(count,1); symmetry = zeros(count,1); cumulative = 0; rec = 1;
for step = 1:duration/dt
    h = h + rain*dt;
    [~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,~] = Local_Inertial_Model_D4( ...
        1,[],[],[],[],[],[],[],[],[],[],[],[],0,z,h*1000,h*1000,n,n.^2,area, ...
        dt/60,dx,outlet,1,downslope,outlet_row,outlet_col,1e-6,outflow,false(ny,nx), ...
        0,0,[],[],zeros(ny,nx),zeros(ny,nx),0,0,0,0,area,[],0,0,[]);
    h = max(double(d_t)/1000,0);
    q_now = sum(double(outlet_flow),'all')/1000/3600*area;
    cumulative = cumulative + q_now*dt;
    if mod(step,record_dt) == 0
        rec = rec+1; q(rec)=q_now; cumulative_out(rec)=cumulative;
        storage(rec)=sum(h,'all')*area;
        input = rain*step*dt*nx*ny*area;
        mass_error_pct(rec)=100*(storage(rec)+cumulative-input)/input;
        symmetry(rec)=sqrt(mean((h(:,1:center-1)-fliplr(h(:,center+1:end))).^2,'all'));
    end
end
Reference.series = table(records/60,q,storage,cumulative_out,mass_error_pct,symmetry, ...
    'VariableNames',{'time_min','outlet_discharge_m3_s','storage_m3','outlet_volume_m3','mass_error_pct','symmetry_rmse_m'});
Reference.final_depth_m = flipud(h);
Reference.rainfall_volume_m3 = rain*duration*nx*ny*area;
Reference.outlet_volume_m3 = cumulative;
Reference.mass_error_pct = mass_error_pct(end);
Reference.cell_count = nx*ny;
end

function [Series, final_map, metrics] = run_voronoi_case(mesh_file, overlap_file, reference, output_file)
mesh = HydroPol2D_Read_UGRID(mesh_file);
mapping = HydroPol2D_Read_Overlap(overlap_file);
outlet = find(mesh.edge_neighbor == 0 & abs(mesh.edge_midpoint_y) < 1e-8 & ...
    mesh.edge_midpoint_x >= 380 & mesh.edge_midpoint_x <= 440);
assert(abs(sum(mesh.edge_length(outlet))-60) <= 1e-8, 'V-tilted outlet width must be exactly 60 m.');
config = struct('duration_s',90*60,'min_dt_s',0.01,'max_dt_s',1, ...
    'output_interval_s',60,'forcing_interval_s',inf,'courant',0.6, ...
    'surface_roughness',0.04,'critical_flow',false,'compute_backend','cpu', ...
    'output_netcdf',output_file,'overwrite_output',true);
forcing = struct('surface_source_m_s',30e-3/3600, ...
    'surface_boundary',struct('edge_id',outlet,'type',"normal_flow",'value',0.004));
run = HydroPol2D_Voronoi_Run(mesh_file,config,forcing);
d = run.diagnostics; diagnostic_time = [d.time_s]'; dt = [d.dt_s]';
diagnostic_q = -[d.boundary_net_inflow_volume_m3]'./dt;
diagnostic_outlet_volume = cumsum(diagnostic_q.*dt);
target_s = reference.series.time_min*60;
q = [0; interp1(diagnostic_time,diagnostic_q,target_s(2:end),'linear')];
outlet_volume = [0; interp1(diagnostic_time,diagnostic_outlet_volume,target_s(2:end),'linear')];
output_storage = sum(run.surface_depth_m.*mesh.surface_area,1)';
storage = interp1([0;run.time_s],[0;output_storage],target_s,'linear');
input = 30e-3/3600.*target_s.*sum(mesh.cell_area);
mass_error_pct = zeros(size(input));
mass_error_pct(2:end) = 100*(storage(2:end)+outlet_volume(2:end)-input(2:end))./input(2:end);

mapped = mapping.mesh_to_raster * run.surface_depth_m;
mapped = reshape(mapped,[41,31,numel(run.time_s)]);
mapped = permute(mapped,[2,1,3]);
symmetry_output = zeros(numel(run.time_s),1);
for i = 1:numel(run.time_s)
    left = mapped(:,1:20,i); right = fliplr(mapped(:,22:end,i));
    symmetry_output(i)=sqrt(mean((left-right).^2,'all'));
end
symmetry = interp1([0;run.time_s],[0;symmetry_output],target_s,'linear');
final_map = mapped(:,:,end);
Series = table(target_s/60,q,storage,outlet_volume,mass_error_pct,symmetry, ...
    'VariableNames',{'time_min','outlet_discharge_m3_s','storage_m3','outlet_volume_m3','mass_error_pct','symmetry_rmse_m'});

q_ref = reference.series.outlet_discharge_m3_s;
q_error = q-q_ref;
depth_error = final_map-reference.final_depth_m;
[peak_q,peak_i]=max(q); [peak_ref,peak_ref_i]=max(q_ref);
threshold = 1e-4;
metrics = struct( ...
    'mesh',string(erase(string(mesh_file),[string(fileparts(mesh_file))+filesep,"vtilted-","-mesh.nc"])), ...
    'cell_count',mesh.n_cells,'cell_reduction_fraction',1-mesh.n_cells/reference.cell_count, ...
    'hydrograph_rmse_m3_s',sqrt(mean(q_error.^2)), ...
    'hydrograph_relative_l2',norm(q_error)/max(norm(q_ref),eps), ...
    'hydrograph_nse',nse(q,q_ref), ...
    'peak_discharge_error_pct',100*abs(peak_q-peak_ref)/max(peak_ref,eps), ...
    'peak_time_error_min',abs(peak_i-peak_ref_i), ...
    'outlet_volume_error_pct',100*abs(outlet_volume(end)-reference.outlet_volume_m3)/reference.outlet_volume_m3, ...
    'final_depth_rmse_m',sqrt(mean(depth_error.^2,'all')), ...
    'final_depth_relative_l2',norm(depth_error(:))/max(norm(reference.final_depth_m(:)),eps), ...
    'inundated_area_error_pct',100*abs(sum(final_map>threshold,'all')-sum(reference.final_depth_m>threshold,'all'))/max(sum(reference.final_depth_m>threshold,'all'),1), ...
    'maximum_symmetry_rmse_m',max(symmetry), ...
    'mass_error_pct',max(abs(mass_error_pct)), ...
    'minimum_dt_s',min(dt),'maximum_dt_s',max(dt));
metrics.passed = metrics.mass_error_pct < 0.1 && metrics.maximum_symmetry_rmse_m < 1e-3 && ...
    metrics.hydrograph_nse >= 0.95 && metrics.peak_discharge_error_pct <= 5;
end

function Comparison = compare_voronoi_meshes(fine,variable,fine_map,variable_map,fine_cells,variable_cells)
q_error = variable.outlet_discharge_m3_s-fine.outlet_discharge_m3_s;
depth_error = variable_map-fine_map;
[peak_variable,peak_variable_i]=max(variable.outlet_discharge_m3_s);
[peak_fine,peak_fine_i]=max(fine.outlet_discharge_m3_s);
correlation = corrcoef(fine_map(:),variable_map(:));
threshold=1e-4;
Comparison = table( ...
    fine_cells,variable_cells,1-variable_cells/fine_cells, ...
    sqrt(mean(q_error.^2)),norm(q_error)/max(norm(fine.outlet_discharge_m3_s),eps), ...
    nse(variable.outlet_discharge_m3_s,fine.outlet_discharge_m3_s), ...
    100*abs(peak_variable-peak_fine)/max(peak_fine,eps),abs(peak_variable_i-peak_fine_i), ...
    100*abs(variable.outlet_volume_m3(end)-fine.outlet_volume_m3(end))/fine.outlet_volume_m3(end), ...
    sqrt(mean(depth_error.^2,'all')),norm(depth_error(:))/max(norm(fine_map(:)),eps),correlation(1,2), ...
    100*abs(sum(variable_map>threshold,'all')-sum(fine_map>threshold,'all'))/max(sum(fine_map>threshold,'all'),1), ...
    'VariableNames',{'fine_cell_count','variable_cell_count','cell_reduction_fraction', ...
    'hydrograph_rmse_m3_s','hydrograph_relative_l2','hydrograph_nse', ...
    'peak_discharge_error_pct','peak_time_error_min','outlet_volume_error_pct', ...
    'final_depth_rmse_m','final_depth_relative_l2','final_depth_correlation','inundated_area_error_pct'});
Comparison.passed = Comparison.outlet_volume_error_pct <= 2 && ...
    Comparison.peak_discharge_error_pct <= 5 && Comparison.final_depth_relative_l2 <= 0.05 && ...
    Comparison.inundated_area_error_pct <= 5;
end

function value = nse(model, reference)
denominator = sum((reference-mean(reference)).^2);
value = 1-sum((model-reference).^2)/max(denominator,eps);
end
