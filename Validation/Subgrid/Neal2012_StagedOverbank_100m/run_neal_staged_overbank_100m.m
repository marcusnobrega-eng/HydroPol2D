function run_neal_staged_overbank_100m()
%RUN_NEAL_STAGED_OVERBANK_100M Test channel, bench, and floodplain activation.

case_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(fileparts(fileparts(case_dir)));
addpath(fullfile(model_dir, 'HydroPol2D_Functions'));

Cfg = default_config();
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
fig_dir = fullfile(case_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

G = build_grids(Cfg);
[T, Mass, Final, Profiles, MaxEta] = run_models(Cfg, G);
Metrics = score_models(Cfg, T, Mass);
MapMetrics = score_max_depth_maps(Cfg, G, MaxEta);
Transitions = transition_table(Cfg, T);
Rating = analytical_rating_curve(Cfg);
Pass = make_pass_table(Cfg, Metrics, MapMetrics, Transitions);

writetable(T, fullfile(out_dir, 'TimeSeries.csv'));
writetable(Mass, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(Metrics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(MapMetrics, fullfile(out_dir, 'MaxDepth_MapMetrics.csv'));
writetable(Transitions, fullfile(out_dir, 'Transition_Timing.csv'));
writetable(Rating, fullfile(out_dir, 'Analytical_Rating_Curve.csv'));
writetable(Pass, fullfile(out_dir, 'Pass_Fail.csv'));

plot_geometry(Cfg, G, Rating, fig_dir);
plot_hydrographs(Cfg, T, fig_dir);
plot_depths_and_wetting(Cfg, T, fig_dir);
plot_storage(Cfg, T, fig_dir);
plot_cross_sections(Cfg, G, T, Profiles, fig_dir);
plot_max_depth_maps(Cfg, G, MaxEta, fig_dir);
plot_stage_discharge(Cfg, T, Rating, fig_dir);

save(fullfile(out_dir, 'Results.mat'), 'Cfg', 'G', 'T', 'Mass', 'Final', ...
    'Profiles', 'MaxEta', 'Metrics', 'MapMetrics', 'Transitions', 'Rating', 'Pass');

fprintf('\n%s complete. Neal improvement pass = %d.\n', Cfg.case_id, Pass.Pass);
disp(Metrics);
disp(MapMetrics);
disp(Transitions);
end

function Cfg = default_config()
Cfg = struct();
Cfg.case_id = "P1-SUBGRID-NEAL-STAGED-OVERBANK-001";
Cfg.length_m = 1000;
Cfg.domain_width_m = 300;
Cfg.fine_dx_m = 5;
Cfg.coarse_dx_m = 100;
Cfg.channel_width_m = 50;
Cfg.bench_width_each_m = 25;
Cfg.channel_bank_height_m = 0.5;
Cfg.external_floodplain_elevation_m = 1.0;
Cfg.slope_mpm = 0.005;
Cfg.channel_manning_n = 0.035;
Cfg.floodplain_manning_n = 0.10;
Cfg.peak_inflow_m3s = 200;
Cfg.nash_time_to_peak_min = 45;
Cfg.nash_shape = 4;
Cfg.duration_min = 240;
Cfg.dt_min = 0.01;
Cfg.record_dt_min = 1;
Cfg.gauge_x_m = 750;
Cfg.midpoint_x_m = 500;
Cfg.depth_tolerance_m = 1e-6;
Cfg.wet_threshold_m = 0.01;
Cfg.map_thresholds_m = [0.01, 0.10];
Cfg.qfun = @(t) nash_hydrograph(t, Cfg.peak_inflow_m3s, ...
    Cfg.nash_time_to_peak_min, Cfg.nash_shape);
end

function G = build_grids(Cfg)
ratio = Cfg.coarse_dx_m / Cfg.fine_dx_m;
assert(ratio == round(ratio), 'Fine-to-coarse ratio must be an integer.');
assert(mod(Cfg.length_m, Cfg.coarse_dx_m) == 0, ...
    'Reach length must be divisible by the coarse resolution.');
assert(mod(Cfg.domain_width_m, Cfg.coarse_dx_m) == 0, ...
    'Domain width must be divisible by the coarse resolution.');

dx = Cfg.fine_dx_m;
x = (dx/2:dx:Cfg.length_m-dx/2)';
y = (-Cfg.domain_width_m/2+dx/2:dx:Cfg.domain_width_m/2-dx/2)';
[X,Y] = meshgrid(x, y);
base = -Cfg.slope_mpm .* X;

channel = abs(Y) < Cfg.channel_width_m / 2;
bench = abs(Y) >= Cfg.channel_width_m / 2 & ...
    abs(Y) < Cfg.channel_width_m / 2 + Cfg.bench_width_each_m;
external = ~(channel | bench);

cross_elevation = zeros(size(Y));
cross_elevation(bench) = Cfg.channel_bank_height_m;
cross_elevation(external) = Cfg.external_floodplain_elevation_m;
z_fine = base + cross_elevation;
roughness_fine = Cfg.floodplain_manning_n .* ones(size(z_fine));
roughness_fine(channel) = Cfg.channel_manning_n;

G.fine = make_grid("fine", dx, x, y, z_fine, roughness_fine, Cfg);
G.fine.channel_rows = find(channel(:,1));
G.fine.bench_rows = find(bench(:,1));
G.fine.external_rows = find(external(:,1));
G.fine.inlet_rows = G.fine.channel_rows;

dx = Cfg.coarse_dx_m;
xc = (dx/2:dx:Cfg.length_m-dx/2)';
yc = (-Cfg.domain_width_m/2+dx/2:dx:Cfg.domain_width_m/2-dx/2)';
[Xc,Yc] = meshgrid(xc, yc);

% Bilinear raster sampling at the coarse-cell centers.
z_coarse = interp2(X, Y, z_fine, Xc, Yc, 'linear');
roughness_coarse = interp2(X, Y, roughness_fine, Xc, Yc, 'linear');
G.coarse = make_grid("coarse", dx, xc, yc, z_coarse, roughness_coarse, Cfg);
[~, center_row] = min(abs(yc));
G.coarse.channel_rows = center_row;
G.coarse.inlet_rows = center_row;

% Neal cell elevation is the internal-bench elevation. The embedded
% rectangular channel occupies the remaining 50 m width below that level.
z_neal = -Cfg.slope_mpm .* Xc + Cfg.external_floodplain_elevation_m;
z_neal(center_row,:) = -Cfg.slope_mpm .* Xc(center_row,:) + ...
    Cfg.channel_bank_height_m;
roughness_neal = Cfg.floodplain_manning_n .* ones(size(z_neal));
G.neal = make_grid("neal", dx, xc, yc, z_neal, roughness_neal, Cfg);
G.neal.channel_rows = center_row;
G.neal.inlet_rows = center_row;
G.neal.River_Width(center_row,:) = Cfg.channel_width_m;
G.neal.River_Depth(center_row,:) = Cfg.channel_bank_height_m;
G.neal.n_channel(center_row,:) = Cfg.channel_manning_n;
end

function grid = make_grid(mode, dx, x, y, z, roughness, Cfg)
grid = struct();
grid.mode = mode;
grid.dx = dx;
grid.area = dx^2;
grid.nx = size(z,2);
grid.ny = size(z,1);
grid.x = x(:);
grid.y = y(:);
grid.z = z;
grid.roughness = roughness;
grid.n_channel = roughness;
grid.n_flood = Cfg.floodplain_manning_n .* ones(size(z));
grid.River_Width = zeros(size(z));
grid.River_Depth = zeros(size(z));
grid.outlet_index = false(size(z));
grid.outlet_type = 1;
grid.slope_outlet = Cfg.slope_mpm;
grid.row_outlet = (1:grid.ny)';
grid.col_outlet = repmat(grid.nx, grid.ny, 1);
grid.channel_rows = ceil(grid.ny/2);
grid.inlet_rows = grid.channel_rows;
end

function [T, Mass, Final, Profiles, MaxEta] = run_models(Cfg, G)
dt_s = Cfg.dt_min * 60;
nsteps = round(Cfg.duration_min / Cfg.dt_min);
record_every = round(Cfg.record_dt_min / Cfg.dt_min);
nrecords = floor(nsteps / record_every) + 1;

S.fine = initialize_state(G.fine);
S.coarse = initialize_state(G.coarse);
S.neal = initialize_state(G.neal);
models = {'fine','coarse','neal'};

initial_storage = zeros(3,1);
cum_out = zeros(3,1);
for k = 1:3
    initial_storage(k) = compute_storage(S.(models{k}), G.(models{k}));
end
cum_input = 0;

R = initialize_records(nrecords);
Profiles.fine = nan(G.fine.ny, nrecords);
Profiles.coarse = nan(G.fine.ny, nrecords);
Profiles.neal = nan(G.fine.ny, nrecords);
MaxEta.fine = -inf(size(G.fine.z));
MaxEta.coarse = -inf(size(G.coarse.z));
MaxEta.neal = -inf(size(G.neal.z));

[R,Profiles,MaxEta] = store_record(R, Profiles, MaxEta, 1, 0, S, G, Cfg);
rec = 1;

for it = 1:nsteps
    qin = Cfg.qfun((it - 1) * Cfg.dt_min);
    for k = 1:3
        name = models{k};
        S.(name) = apply_inflow(S.(name), G.(name), qin, dt_s);
        S.(name) = step_model(S.(name), G.(name), Cfg.dt_min);
        qout = outlet_discharge(S.(name), G.(name));
        cum_out(k) = cum_out(k) + qout * dt_s;
    end
    cum_input = cum_input + qin * dt_s;

    if mod(it, record_every) == 0 || it == nsteps
        rec = rec + 1;
        [R,Profiles,MaxEta] = store_record(R, Profiles, MaxEta, rec, ...
            it * Cfg.dt_min, S, G, Cfg);
    end
end

T = struct2table(R);
final_storage = [T.Storage_Fine_m3(end); T.Storage_Coarse_m3(end); ...
    T.Storage_Neal_m3(end)];
available = initial_storage + cum_input;
accounted = cum_out + final_storage;
residual = available - accounted;

Mass = table(repmat(Cfg.case_id,3,1), ...
    ["FineExplicit5m"; "CoarseOrdinary100m"; "NealSubgrid100m"], ...
    repmat(cum_input,3,1), initial_storage, available, cum_out, ...
    final_storage, accounted, residual, 100 .* residual ./ max(cum_input,1), ...
    'VariableNames', {'CaseID','Model','InputVolume_m3','InitialStorage_m3', ...
    'AvailableVolume_m3','OutletVolume_m3','FinalStorage_m3', ...
    'AccountedVolume_m3','Residual_m3','MassResidual_pct'});

Final = S;
end

function R = initialize_records(n)
names = {'Time_min','Inflow_m3s', ...
    'Outlet_Fine_m3s','Outlet_Coarse_m3s','Outlet_Neal_m3s', ...
    'Gauge_Fine_m3s','Gauge_Coarse_m3s','Gauge_Neal_m3s', ...
    'Storage_Fine_m3','Storage_Coarse_m3','Storage_Neal_m3', ...
    'ChannelDepth_Fine_m','ChannelDepth_Coarse_m','ChannelDepth_Neal_m', ...
    'BenchDepth_Fine_m','BenchDepth_Coarse_m','BenchDepth_Neal_m', ...
    'ExternalDepth_Fine_m','ExternalDepth_Coarse_m','ExternalDepth_Neal_m', ...
    'WetArea_Fine_m2','WetArea_Coarse_m2','WetArea_Neal_m2', ...
    'ExternalWetArea_Fine_m2','ExternalWetArea_Coarse_m2','ExternalWetArea_Neal_m2'};
for i = 1:numel(names)
    R.(names{i}) = zeros(n,1);
end
end

function state = initialize_state(grid)
state = struct();
state.d_tot = zeros(grid.ny,grid.nx);
state.d_p = zeros(grid.ny,grid.nx);
state.outflow = zeros(grid.ny,grid.nx,5);
state.outlet_flow = zeros(grid.ny,grid.nx);
state.Qc = zeros(grid.ny,grid.nx,2);
state.Qf = zeros(grid.ny,grid.nx,2);
state.Qci = zeros(grid.ny,grid.nx,2);
state.Qfi = zeros(grid.ny,grid.nx,2);
if grid.mode == "neal"
    state.C_a = hp2d_neal_cell_area(zeros(grid.ny,grid.nx), ...
        grid.River_Width, grid.River_Depth, grid.dx);
else
    state.C_a = grid.area .* ones(grid.ny,grid.nx);
end
end

function state = apply_inflow(state, grid, q_m3s, dt_s)
if q_m3s <= 0; return; end
rows = grid.inlet_rows(:);
volume = q_m3s * dt_s / numel(rows);
if grid.mode == "neal"
    h = state.d_tot(rows,1) ./ 1000;
    h = hp2d_neal_apply_volume_change(h, volume .* ones(size(h)), ...
        grid.River_Width(rows,1), grid.River_Depth(rows,1), grid.dx);
    state.d_tot(rows,1) = 1000 .* h;
else
    state.d_tot(rows,1) = state.d_tot(rows,1) + 1000 .* volume ./ grid.area;
end
end

function state = step_model(state, grid, dt_min)
idx_nan = false(grid.ny,grid.nx);
if grid.mode == "neal"
    [~,~,~,~,outlet_flow,d_t,~,outflow,~,Qc,Qf,Qci,Qfi,C_a] = ...
        Local_Inertial_Model_D4(1,[],[],[],[],[],[],[],[],[],[],[],[], ...
        0,grid.z,state.d_tot,state.d_p,grid.roughness,grid.roughness.^2, ...
        grid.area,dt_min,grid.dx,grid.outlet_index,grid.outlet_type, ...
        grid.slope_outlet,grid.row_outlet,grid.col_outlet,1e-6, ...
        state.outflow,idx_nan,0,1,grid.n_channel,grid.n_flood, ...
        grid.River_Width,grid.River_Depth,state.Qc,state.Qf,state.Qci, ...
        state.Qfi,state.C_a,[],1,0,[]);
    state.Qc = Qc; state.Qf = Qf; state.Qci = Qci; state.Qfi = Qfi;
else
    [~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = ...
        Local_Inertial_Model_D4(1,[],[],[],[],[],[],[],[],[],[],[],[], ...
        0,grid.z,state.d_tot,state.d_p,grid.roughness,grid.roughness.^2, ...
        grid.area,dt_min,grid.dx,grid.outlet_index,grid.outlet_type, ...
        grid.slope_outlet,grid.row_outlet,grid.col_outlet,1e-6, ...
        state.outflow,idx_nan,0,0,grid.n_channel,grid.n_flood, ...
        grid.River_Width,grid.River_Depth,[],[],[],[],state.C_a,[],0,0,[]);
end
state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.outlet_flow = outlet_flow;
state.C_a = C_a;
end

function [R,Profiles,MaxEta] = store_record(R, Profiles, MaxEta, i, time, S, G, Cfg)
R.Time_min(i) = time;
R.Inflow_m3s(i) = Cfg.qfun(time);

models = {'fine','coarse','neal'};
labels = {'Fine','Coarse','Neal'};
for k = 1:3
    name = models{k}; label = labels{k};
    R.(['Outlet_' label '_m3s'])(i) = outlet_discharge(S.(name),G.(name));
    R.(['Gauge_' label '_m3s'])(i) = section_discharge(S.(name),G.(name),Cfg.gauge_x_m);
    R.(['Storage_' label '_m3'])(i) = compute_storage(S.(name),G.(name));
    D = physical_diagnostics(S.(name),G.(name),G.fine,Cfg,Cfg.gauge_x_m);
    R.(['ChannelDepth_' label '_m'])(i) = D.channel_depth_m;
    R.(['BenchDepth_' label '_m'])(i) = D.bench_depth_m;
    R.(['ExternalDepth_' label '_m'])(i) = D.external_depth_m;
    R.(['WetArea_' label '_m2'])(i) = D.wet_area_m2;
    R.(['ExternalWetArea_' label '_m2'])(i) = D.external_wet_area_m2;
    Profiles.(name)(:,i) = D.profile_eta_m;

    eta = model_wse(S.(name),G.(name));
    wet = S.(name).d_tot > Cfg.depth_tolerance_m * 1000;
    MaxEta.(name)(wet) = max(MaxEta.(name)(wet),eta(wet));
end
end

function D = physical_diagnostics(state, grid, fine_grid, Cfg, x_gauge)
eta = project_eta_to_fine(state,grid,fine_grid,Cfg);
[~,col] = min(abs(fine_grid.x-x_gauge));
depth = max(eta-fine_grid.z,0);
depth(~isfinite(eta)) = 0;

y = fine_grid.y;
channel = abs(y) < Cfg.channel_width_m/2;
bench = abs(y) >= Cfg.channel_width_m/2 & ...
    abs(y) < Cfg.channel_width_m/2+Cfg.bench_width_each_m;
external = ~(channel|bench);

D.channel_depth_m = max(depth(channel,col),[],'omitnan');
D.bench_depth_m = mean(depth(bench,col),'omitnan');
D.external_depth_m = mean(depth(external,col),'omitnan');
D.wet_area_m2 = nnz(depth > Cfg.wet_threshold_m) * fine_grid.area;
D.external_wet_area_m2 = nnz(depth(external,:) > Cfg.wet_threshold_m) * fine_grid.area;
D.profile_eta_m = eta(:,col);
end

function eta_fine = project_eta_to_fine(state,grid,fine_grid,Cfg)
eta = model_wse(state,grid);
if grid.mode == "fine"
    eta(state.d_tot <= Cfg.depth_tolerance_m*1000) = NaN;
    eta_fine = eta;
else
    wet = state.d_tot > Cfg.depth_tolerance_m*1000;
    eta_fine = project_coarse_rows(eta,wet,grid,fine_grid);
end
end

function eta_fine = project_coarse_rows(eta,wet,grid,fine_grid)
% Reconstruct longitudinal water levels without smoothing cross-section steps.
ratio = grid.dx/fine_grid.dx;
bed = model_bed(grid);
eta_fill = eta;
eta_fill(~wet) = bed(~wet);
eta_x = zeros(grid.ny,fine_grid.nx);
for row = 1:grid.ny
    eta_x(row,:) = interp1(grid.x,eta_fill(row,:),fine_grid.x','linear','extrap');
end
eta_fine = repelem(eta_x,ratio,1);
wet_fine = repelem(wet,ratio,ratio);
eta_fine(~wet_fine) = NaN;
end

function bed = model_bed(grid)
if grid.mode == "neal"
    bed = grid.z-grid.River_Depth;
else
    bed = grid.z;
end
end

function eta = model_wse(state,grid)
h = max(state.d_tot./1000,0);
if grid.mode == "neal"
    eta = grid.z-grid.River_Depth+h;
else
    eta = grid.z+h;
end
end

function q = outlet_discharge(state,grid)
q = sum(state.outlet_flow(:),'omitnan') * grid.area / 1000 / 3600;
end

function q = section_discharge(state,grid,x_gauge)
[~,col] = min(abs(grid.x-x_gauge));
col = min(col,grid.nx-1);
q = sum(state.outflow(:,col,1),'omitnan') * grid.area / 1000 / 3600;
end

function storage = compute_storage(state,grid)
h = max(state.d_tot./1000,0);
if grid.mode == "neal"
    V = hp2d_neal_cell_volume(h,grid.River_Width,grid.River_Depth,grid.dx);
else
    V = h.*grid.area;
end
storage = sum(V(:),'omitnan');
end

function Metrics = score_models(Cfg,T,Mass)
models = ["FineExplicit5m";"CoarseOrdinary100m";"NealSubgrid100m"];
qout = {T.Outlet_Fine_m3s,T.Outlet_Coarse_m3s,T.Outlet_Neal_m3s};
qg = {T.Gauge_Fine_m3s,T.Gauge_Coarse_m3s,T.Gauge_Neal_m3s};
storage = {T.Storage_Fine_m3,T.Storage_Coarse_m3,T.Storage_Neal_m3};
channel = {T.ChannelDepth_Fine_m,T.ChannelDepth_Coarse_m,T.ChannelDepth_Neal_m};
bench = {T.BenchDepth_Fine_m,T.BenchDepth_Coarse_m,T.BenchDepth_Neal_m};
external = {T.ExternalDepth_Fine_m,T.ExternalDepth_Coarse_m,T.ExternalDepth_Neal_m};

n = 3;
out_rmse = zeros(n,1); out_nse = ones(n,1); gauge_rmse = zeros(n,1);
gauge_nse = ones(n,1); peak_error = zeros(n,1); peak_time = zeros(n,1);
out_vol_error = zeros(n,1); storage_rmse = zeros(n,1);
channel_rmse = zeros(n,1); bench_rmse = zeros(n,1); external_rmse = zeros(n,1);
[ref_peak,ref_peak_idx] = max(qout{1});
ref_out_vol = trapz(T.Time_min*60,qout{1});
for i = 1:n
    out_rmse(i) = rmse(qout{1},qout{i});
    out_nse(i) = nse(qout{1},qout{i});
    gauge_rmse(i) = rmse(qg{1},qg{i});
    gauge_nse(i) = nse(qg{1},qg{i});
    [peak,idx] = max(qout{i});
    peak_error(i) = 100*(peak-ref_peak)/max(ref_peak,1);
    peak_time(i) = T.Time_min(idx)-T.Time_min(ref_peak_idx);
    out_vol_error(i) = 100*(trapz(T.Time_min*60,qout{i})-ref_out_vol)/max(ref_out_vol,1);
    storage_rmse(i) = rmse(storage{1},storage{i});
    channel_rmse(i) = rmse(channel{1},channel{i});
    bench_rmse(i) = rmse(bench{1},bench{i});
    external_rmse(i) = rmse(external{1},external{i});
end

Metrics = table(repmat(Cfg.case_id,n,1),models,out_rmse,out_nse,gauge_rmse, ...
    gauge_nse,peak_time,peak_error,out_vol_error,storage_rmse,channel_rmse, ...
    bench_rmse,external_rmse,Mass.MassResidual_pct, ...
    'VariableNames',{'CaseID','Model','OutletRMSE_m3s','OutletNSE', ...
    'InternalGaugeRMSE_m3s','InternalGaugeNSE','PeakTimingError_min', ...
    'PeakDischargeError_pct','OutletVolumeError_pct','StorageRMSE_m3', ...
    'ChannelDepthRMSE_m','BenchDepthRMSE_m','ExternalDepthRMSE_m', ...
    'MassResidual_pct'});
end

function MapMetrics = score_max_depth_maps(Cfg,G,MaxEta)
ref = max(MaxEta.fine-G.fine.z,0);
ref(~isfinite(MaxEta.fine)) = 0;
models = ["FineExplicit5m";"CoarseOrdinary100m";"NealSubgrid100m"];
depths = cell(3,1);
depths{1} = ref;
depths{2} = projected_max_depth(MaxEta.coarse,G.coarse,G.fine);
depths{3} = projected_max_depth(MaxEta.neal,G.neal,G.fine);

map_rmse = zeros(3,1); map_mae = zeros(3,1); bias = zeros(3,1);
wet_error_001 = zeros(3,1); wet_error_010 = zeros(3,1);
csi_001 = zeros(3,1); csi_010 = zeros(3,1); volume_error = zeros(3,1);
ref_volume = sum(ref(:))*G.fine.area;
for i = 1:3
    d = depths{i}; diff = d-ref;
    map_rmse(i) = sqrt(mean(diff(:).^2));
    map_mae(i) = mean(abs(diff(:)));
    bias(i) = mean(diff(:));
    [wet_error_001(i),csi_001(i)] = wet_scores(ref,d,Cfg.map_thresholds_m(1),G.fine.area);
    [wet_error_010(i),csi_010(i)] = wet_scores(ref,d,Cfg.map_thresholds_m(2),G.fine.area);
    volume_error(i) = 100*(sum(d(:))*G.fine.area-ref_volume)/max(ref_volume,1);
end
MapMetrics = table(repmat(Cfg.case_id,3,1),models,map_rmse,map_mae,bias, ...
    wet_error_001,csi_001,wet_error_010,csi_010,volume_error, ...
    'VariableNames',{'CaseID','Model','MaxDepthRMSE_m','MaxDepthMAE_m', ...
    'MaxDepthBias_m','WetAreaError_001m_pct','CSI_001m', ...
    'WetAreaError_010m_pct','CSI_010m','DepthVolumeError_pct'});
end

function d = projected_max_depth(max_eta,grid,fine_grid)
wet = isfinite(max_eta);
eta = project_coarse_rows(max_eta,wet,grid,fine_grid);
d = max(eta-fine_grid.z,0);
d(~isfinite(eta)) = 0;
end

function [area_error,csi] = wet_scores(ref,sim,threshold,area)
r = ref>threshold; s = sim>threshold;
area_error = 100*(nnz(s)*area-nnz(r)*area)/max(nnz(r)*area,1);
csi = nnz(r&s)/max(nnz(r|s),1);
end

function Transitions = transition_table(Cfg,T)
models = ["FineExplicit5m";"CoarseOrdinary100m";"NealSubgrid100m"];
bench = {T.BenchDepth_Fine_m,T.BenchDepth_Coarse_m,T.BenchDepth_Neal_m};
external = {T.ExternalDepth_Fine_m,T.ExternalDepth_Coarse_m,T.ExternalDepth_Neal_m};
bench_wet = zeros(3,1); external_wet = zeros(3,1);
bench_dry = zeros(3,1); external_dry = zeros(3,1);
for i=1:3
    bench_wet(i) = first_crossing(T.Time_min,bench{i},Cfg.wet_threshold_m);
    external_wet(i) = first_crossing(T.Time_min,external{i},Cfg.wet_threshold_m);
    bench_dry(i) = recession_crossing(T.Time_min,bench{i},Cfg.wet_threshold_m);
    external_dry(i) = recession_crossing(T.Time_min,external{i},Cfg.wet_threshold_m);
end
Transitions = table(repmat(Cfg.case_id,3,1),models,bench_wet,external_wet, ...
    bench_dry,external_dry, ...
    'VariableNames',{'CaseID','Model','BenchFirstWet_min','ExternalFirstWet_min', ...
    'BenchRecessionDry_min','ExternalRecessionDry_min'});
end

function t = first_crossing(time,value,threshold)
idx = find(value>threshold,1,'first');
if isempty(idx); t=NaN; else; t=time(idx); end
end

function t = recession_crossing(time,value,threshold)
if max(value,[],'omitnan')<=threshold
    t=NaN;
    return;
end
[~,peak] = max(value);
idx = find(value(peak:end)<=threshold,1,'first');
if isempty(idx); t=NaN; else; t=time(peak+idx-1); end
end

function Rating = analytical_rating_curve(Cfg)
h = linspace(0,2,401)';
Ac = Cfg.channel_width_m.*h;
Pc = Cfg.channel_width_m+2.*min(h,Cfg.channel_bank_height_m);
Kc = Ac.*(Ac./max(Pc,eps)).^(2/3)./Cfg.channel_manning_n;
hb = max(h-Cfg.channel_bank_height_m,0);
Kb = (2*Cfg.bench_width_each_m).*hb.^(5/3)./Cfg.floodplain_manning_n;
he = max(h-Cfg.external_floodplain_elevation_m,0);
external_width = Cfg.domain_width_m-Cfg.channel_width_m-2*Cfg.bench_width_each_m;
Ke = external_width.*he.^(5/3)./Cfg.floodplain_manning_n;
Q = (Kc+Kb+Ke).*sqrt(Cfg.slope_mpm);
Rating = table(h,Q,Kc.*sqrt(Cfg.slope_mpm),Kb.*sqrt(Cfg.slope_mpm), ...
    Ke.*sqrt(Cfg.slope_mpm), ...
    'VariableNames',{'DepthAboveChannelBed_m','Discharge_m3s', ...
    'ChannelDischarge_m3s','BenchDischarge_m3s','ExternalFloodplainDischarge_m3s'});
end

function Pass = make_pass_table(Cfg,Metrics,MapMetrics,Transitions)
coarse = Metrics(Metrics.Model=="CoarseOrdinary100m",:);
neal = Metrics(Metrics.Model=="NealSubgrid100m",:);
coarse_map = MapMetrics(MapMetrics.Model=="CoarseOrdinary100m",:);
neal_map = MapMetrics(MapMetrics.Model=="NealSubgrid100m",:);
ref_transition = Transitions(Transitions.Model=="FineExplicit5m",:);
neal_transition = Transitions(Transitions.Model=="NealSubgrid100m",:);
mass_pass = abs(neal.MassResidual_pct)<0.1;
hydrograph_pass = neal.OutletRMSE_m3s<coarse.OutletRMSE_m3s && ...
    neal.InternalGaugeRMSE_m3s<coarse.InternalGaugeRMSE_m3s;
peak_extent_pass = neal.BenchDepthRMSE_m<coarse.BenchDepthRMSE_m && ...
    neal_map.MaxDepthRMSE_m<coarse_map.MaxDepthRMSE_m;
activation_pass = abs(neal_transition.BenchFirstWet_min-ref_transition.BenchFirstWet_min)<=5 && ...
    abs(neal_transition.ExternalFirstWet_min-ref_transition.ExternalFirstWet_min)<=5;
recession_pass = abs(neal_transition.BenchRecessionDry_min-ref_transition.BenchRecessionDry_min)<=15 && ...
    abs(neal_transition.ExternalRecessionDry_min-ref_transition.ExternalRecessionDry_min)<=15;
passed = mass_pass && hydrograph_pass && peak_extent_pass && activation_pass && recession_pass;
notes = sprintf(['Neal/ordinary outlet RMSE %.3f/%.3f m3/s; bench-depth RMSE ', ...
    '%.3f/%.3f m; max-depth RMSE %.3f/%.3f m; recession timing pass %d.'], ...
    neal.OutletRMSE_m3s,coarse.OutletRMSE_m3s,neal.BenchDepthRMSE_m, ...
    coarse.BenchDepthRMSE_m,neal_map.MaxDepthRMSE_m,coarse_map.MaxDepthRMSE_m,recession_pass);
Pass = table(Cfg.case_id,logical(mass_pass),logical(hydrograph_pass), ...
    logical(peak_extent_pass),logical(activation_pass),logical(recession_pass), ...
    logical(passed),string(notes), ...
    'VariableNames',{'CaseID','MassBalancePass','HydrographPass', ...
    'PeakExtentPass','ActivationTimingPass','RecessionTimingPass','Pass','Notes'});
end

function plot_geometry(Cfg,G,Rating,fig_dir)
colors = model_colors();
f = figure('Visible','off','Color','w','Position',[100 100 1100 720]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
[y,z] = step_profile(G.fine.y,G.fine.z(:,1)-min(G.fine.z(:,1)),G.fine.dx);
plot(y,z,'k-','LineWidth',2); hold on;
xline([-50 50],'Color',[0.4 0.4 0.4],'LineStyle','--');
xline([-25 25],'Color',[0.4 0.4 0.4],'LineStyle',':');
patch([-25 25 25 -25],[0 0 0.5 0.5],colors.channel,'FaceAlpha',0.18,'EdgeColor','none');
patch([-50 -25 -25 -50],[0.5 0.5 1 1],colors.neal,'FaceAlpha',0.14,'EdgeColor','none');
patch([25 50 50 25],[0.5 0.5 1 1],colors.neal,'FaceAlpha',0.14,'EdgeColor','none');
xlabel('Cross-channel distance (m)'); ylabel('Elevation above channel bed (m)');
title('Staged compound cross section'); grid on; xlim([-150 150]); ylim([-0.05 1.35]);
text(0,0.16,'50 m channel, n = 0.035','HorizontalAlignment','center');
text(-37.5,0.65,'25 m bench','HorizontalAlignment','center');
text(37.5,0.65,'25 m bench','HorizontalAlignment','center');
text(-100,1.10,'external floodplain, n = 0.10','HorizontalAlignment','center');
text(100,1.10,'external floodplain, n = 0.10','HorizontalAlignment','center');

nexttile;
plot(Rating.DepthAboveChannelBed_m,Rating.Discharge_m3s,'k-','LineWidth',2); hold on;
xline(Cfg.channel_bank_height_m,':','Internal benches activate','LineWidth',1.2);
xline(Cfg.external_floodplain_elevation_m,'--','External floodplain activates','LineWidth',1.2);
yline(Cfg.peak_inflow_m3s,'Color',[0.5 0.5 0.5],'LineStyle','-.');
grid on; xlim([0 1.6]);
xlabel('Depth above channel bed (m)'); ylabel('Uniform-flow discharge (m^3 s^{-1})');
title(sprintf('Analytical conveyance at slope %.3f',Cfg.slope_mpm));
exportgraphics(f,fullfile(fig_dir,'Geometry_and_AnalyticalRating.png'),'Resolution',240);
close(f);
end

function plot_hydrographs(Cfg,T,fig_dir)
c = model_colors();
f = figure('Visible','off','Color','w','Position',[100 100 1050 720]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(T.Time_min,T.Inflow_m3s,'k--','LineWidth',1.5); hold on;
plot_models(T.Time_min,T.Outlet_Fine_m3s,T.Outlet_Coarse_m3s,T.Outlet_Neal_m3s,c);
grid on; xlim([0 Cfg.duration_min]); ylabel('Discharge (m^3 s^{-1})');
title('Downstream boundary');
legend({'Inflow','Fine explicit, 5 m','Ordinary, 100 m','Neal, 100 m'},'Location','northeast');
nexttile;
plot_models(T.Time_min,T.Gauge_Fine_m3s,T.Gauge_Coarse_m3s,T.Gauge_Neal_m3s,c);
grid on; xlim([0 Cfg.duration_min]); xlabel('Time (min)'); ylabel('Discharge (m^3 s^{-1})');
title(sprintf('Internal section at x = %.0f m',Cfg.gauge_x_m));
legend({'Fine explicit, 5 m','Ordinary, 100 m','Neal, 100 m'},'Location','northeast');
exportgraphics(f,fullfile(fig_dir,'Hydrograph_Comparison.png'),'Resolution',240);
close(f);
end

function plot_depths_and_wetting(Cfg,T,fig_dir)
c = model_colors();
f = figure('Visible','off','Color','w','Position',[100 100 1100 850]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot_models(T.Time_min,T.ChannelDepth_Fine_m,T.ChannelDepth_Coarse_m,T.ChannelDepth_Neal_m,c); hold on;
yline(0.5,':','Internal benches','LineWidth',1.1); yline(1,'--','External floodplain','LineWidth',1.1);
grid on; xlim([0 Cfg.duration_min]); ylabel('Channel depth (m)'); title('Water depth above the channel bed');
legend({'Fine explicit, 5 m','Ordinary, 100 m','Neal, 100 m'},'Location','northeast');
nexttile;
plot_models(T.Time_min,T.BenchDepth_Fine_m,T.BenchDepth_Coarse_m,T.BenchDepth_Neal_m,c);
grid on; xlim([0 Cfg.duration_min]); ylabel('Bench depth (m)'); title('Mean depth on the two internal benches');
nexttile;
plot_models(T.Time_min,T.ExternalWetArea_Fine_m2/1e4,T.ExternalWetArea_Coarse_m2/1e4, ...
    T.ExternalWetArea_Neal_m2/1e4,c);
grid on; xlim([0 Cfg.duration_min]); xlabel('Time (min)'); ylabel('Wet area (ha)');
title('External floodplain area deeper than 0.01 m');
exportgraphics(f,fullfile(fig_dir,'Depth_and_FloodplainActivation.png'),'Resolution',240);
close(f);
end

function plot_storage(Cfg,T,fig_dir)
c = model_colors();
f = figure('Visible','off','Color','w','Position',[100 100 1050 530]);
plot_models(T.Time_min,T.Storage_Fine_m3,T.Storage_Coarse_m3,T.Storage_Neal_m3,c);
grid on; xlim([0 Cfg.duration_min]); xlabel('Time (min)'); ylabel('Stored water (m^3)');
title('Domain water storage');
legend({'Fine explicit, 5 m','Ordinary, 100 m','Neal, 100 m'},'Location','northeast');
exportgraphics(f,fullfile(fig_dir,'Storage_Comparison.png'),'Resolution',240);
close(f);
end

function plot_cross_sections(Cfg,G,T,Profiles,fig_dir)
c = model_colors();
idx = snapshot_indices(T);
f = figure('Visible','off','Color','w','Position',[100 100 1250 780]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
[~,gauge_col] = min(abs(G.fine.x-Cfg.gauge_x_m));
datum = min(G.fine.z(:,gauge_col));
bed = G.fine.z(:,gauge_col)-datum;
[yb,zb] = step_profile(G.fine.y,bed,G.fine.dx);
for k=1:3
    nexttile;
    plot(yb,zb,'k-','LineWidth',1.5); hold on;
    plot_wse_profile(G.fine.y,Profiles.fine(:,idx(k))-datum,G.fine.dx,c.fine,2.0);
    plot_wse_profile(G.fine.y,Profiles.coarse(:,idx(k))-datum,G.fine.dx,c.coarse,1.8);
    plot_wse_profile(G.fine.y,Profiles.neal(:,idx(k))-datum,G.fine.dx,c.neal,2.0);
    grid on; xlim([-150 150]); ylim([-0.05 1.8]); ylabel('Elevation (m)');
    title(sprintf('t = %.0f min',T.Time_min(idx(k))));
    if k==1
        legend({'Physical terrain','Fine explicit, 5 m','Ordinary, 100 m','Neal, 100 m'}, ...
            'Location','northoutside','NumColumns',4);
    end
end
xlabel('Cross-channel distance (m)');
exportgraphics(f,fullfile(fig_dir,'CrossSection_Evolution.png'),'Resolution',240);
close(f);
end

function idx = snapshot_indices(T)
[~,peak] = max(T.ChannelDepth_Fine_m);
rise = find(T.ChannelDepth_Fine_m>=0.75,1,'first');
recession = find((1:height(T))'>peak & T.ChannelDepth_Fine_m<=0.75 & ...
    T.ChannelDepth_Fine_m>0.05,1,'first');
if isempty(rise); rise=max(1,round(peak/2)); end
if isempty(recession); recession=min(height(T),peak+round((height(T)-peak)/2)); end
idx = [rise peak recession];
end

function plot_max_depth_maps(Cfg,G,MaxEta,fig_dir)
cmap = parula(256);
ref = max(MaxEta.fine-G.fine.z,0); ref(~isfinite(MaxEta.fine))=0;
ordinary = projected_max_depth(MaxEta.coarse,G.coarse,G.fine);
neal = projected_max_depth(MaxEta.neal,G.neal,G.fine);
limit = max([ref(:);ordinary(:);neal(:)]);
f = figure('Visible','off','Color','w','Position',[100 100 1300 430]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
maps = {ref,ordinary,neal}; names = {'Fine explicit, 5 m','Ordinary, 100 m','Neal, 100 m'};
for i=1:3
    nexttile; imagesc(G.fine.x,G.fine.y,maps{i}); axis xy image; clim([0 limit]); colormap(cmap);
    title(names{i}); if i==1; ylabel('Cross-channel distance (m)'); end
    xlabel('Distance downstream (m)');
end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='Maximum depth (m)';
exportgraphics(f,fullfile(fig_dir,'MaximumDepth_Maps.png'),'Resolution',240);
close(f);

f = figure('Visible','off','Color','w','Position',[100 100 1000 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
difflim = max(abs([ordinary(:)-ref(:);neal(:)-ref(:)]));
nexttile; imagesc(G.fine.x,G.fine.y,ordinary-ref); axis xy image; clim([-difflim difflim]);
colormap(gca,bluewhitered()); title('Ordinary minus fine'); xlabel('Distance downstream (m)');
ylabel('Cross-channel distance (m)');
nexttile; imagesc(G.fine.x,G.fine.y,neal-ref); axis xy image; clim([-difflim difflim]);
colormap(gca,bluewhitered()); title('Neal minus fine'); xlabel('Distance downstream (m)');
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='Maximum-depth difference (m)';
exportgraphics(f,fullfile(fig_dir,'MaximumDepth_Differences.png'),'Resolution',240);
close(f);
end

function plot_stage_discharge(Cfg,T,Rating,fig_dir)
c = model_colors();
f = figure('Visible','off','Color','w','Position',[100 100 800 620]);
plot(Rating.DepthAboveChannelBed_m,Rating.Discharge_m3s,'k-','LineWidth',2); hold on;
scatter(T.ChannelDepth_Fine_m,T.Gauge_Fine_m3s,16,T.Time_min,'filled','Marker','o');
plot(T.ChannelDepth_Neal_m,T.Gauge_Neal_m3s,'Color',c.neal,'LineWidth',1.4);
xline(0.5,':'); xline(1,'--'); grid on; xlim([0 1.6]);
xlabel('Depth above channel bed (m)'); ylabel('Discharge (m^3 s^{-1})');
title(sprintf('Stage-discharge response at x = %.0f m',Cfg.gauge_x_m));
legend({'Analytical uniform flow','Fine transient states','Neal transient path'},'Location','northwest');
cb=colorbar; cb.Label.String='Fine-model time (min)';
exportgraphics(f,fullfile(fig_dir,'Stage_Discharge_Response.png'),'Resolution',240);
close(f);
end

function plot_models(t,a,b,cval,c)
plot(t,a,'Color',c.fine,'LineWidth',2); hold on;
plot(t,b,'Color',c.coarse,'LineWidth',1.8);
plot(t,cval,'Color',c.neal,'LineWidth',2);
end

function plot_wse_profile(y,wse,dx,color,width)
[ys,ws] = step_profile(y,wse,dx);
plot(ys,ws,'Color',color,'LineWidth',width);
end

function c = model_colors()
c.fine = [0.00 0.35 0.70];
c.coarse = [0.85 0.40 0.05];
c.neal = [0.10 0.50 0.25];
c.channel = [0.16 0.45 0.70];
end

function [xstep,ystep] = step_profile(centers,values,dx)
centers=centers(:); values=values(:);
xstep=reshape([(centers-dx/2)';(centers+dx/2)'],[],1);
ystep=repelem(values,2);
end

function q = nash_hydrograph(time_min,peak,time_to_peak,shape)
q=zeros(size(time_min)); positive=time_min>0; p=shape-1;
x=time_min(positive)./time_to_peak;
q(positive)=peak.*x.^p.*exp(p.*(1-x));
end

function value = rmse(a,b)
value=sqrt(mean((double(a(:))-double(b(:))).^2,'omitnan'));
end

function value = nse(obs,sim)
obs=double(obs(:)); sim=double(sim(:));
den=sum((obs-mean(obs,'omitnan')).^2,'omitnan');
if den<=eps; value=double(all(abs(obs-sim)<1e-12));
else; value=1-sum((obs-sim).^2,'omitnan')/den; end
end

function map = bluewhitered()
n=256; half=n/2;
map=[linspace(0.10,1,half)' linspace(0.35,1,half)' ones(half,1); ...
    ones(half,1) linspace(1,0.20,half)' linspace(1,0.10,half)'];
end
