function run_neal_2d_cross_junction_validation(flow_regime)
%RUN_NEAL_2D_CROSS_JUNCTION_VALIDATION Steady three-inlet Neal benchmark.
%   Compares a 10 m explicit local-inertial reference, a 30 m ordinary
%   local-inertial grid, and a 30 m Neal channel-subgrid grid. Equal steady
%   flows enter from the west, north, and south and leave through the east.

if nargin < 1
    flow_regime = "overbank";
end
case_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(fileparts(fileparts(case_dir)));
addpath(fullfile(model_dir, 'HydroPol2D_Functions'));

Cfg = default_config(flow_regime);
if Cfg.flow_regime == "inbank"
    out_dir = fullfile(case_dir, 'Outputs', 'Validation', 'Inbank50');
    fig_dir = fullfile(case_dir, 'Figures', 'Inbank50');
else
    out_dir = fullfile(case_dir, 'Outputs', 'Validation');
    fig_dir = fullfile(case_dir, 'Figures');
end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

G = build_grids(Cfg);
[T, Mass, MaxEta, Snapshots, Final] = run_models(Cfg, G);
Metrics = score_hydrographs(Cfg, T, Mass);
[MapMetrics, MaxDepth] = score_maximum_depth(Cfg, G, MaxEta);
Steady = score_steady_state(Cfg, T);
Confinement = score_confinement(Cfg, G, MaxDepth);
Pass = make_pass_table(Cfg, Metrics, MapMetrics, Mass, Steady, Confinement);

writetable(T, fullfile(out_dir, 'TimeSeries.csv'));
writetable(Mass, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(Metrics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(MapMetrics, fullfile(out_dir, 'MaxDepth_MapMetrics.csv'));
writetable(Steady, fullfile(out_dir, 'Steady_State_Summary.csv'));
writetable(Confinement, fullfile(out_dir, 'Confinement_Summary.csv'));
writetable(Pass, fullfile(out_dir, 'Pass_Fail.csv'));

plot_geometry(Cfg, G, fig_dir);
plot_hydrographs(T, fig_dir);
plot_storage_wet_area(T, fig_dir);
plot_maximum_depth(Cfg, G, MaxDepth, fig_dir);
plot_snapshots(Cfg, G, Snapshots, fig_dir);

save(fullfile(out_dir, 'Results.mat'), 'Cfg', 'G', 'T', 'Mass', ...
    'Metrics', 'MapMetrics', 'Steady', 'Confinement', 'Pass', ...
    'MaxEta', 'MaxDepth', 'Snapshots', ...
    'Final', '-v7.3');

fprintf('\n%s complete. Overall pass = %d.\n', Cfg.case_id, Pass.Pass);
disp(Metrics);
disp(MapMetrics);
disp(Steady);
disp(Confinement);
disp(Mass);
disp(Pass);
end

function Steady = score_steady_state(Cfg,T)
window = T.Time_min >= Cfg.duration_min-20;
models = ["FineExplicit10m";"CoarseOrdinary30m";"NealSubgrid30m"];
labels = {'Fine','Coarse','Neal'};
n = numel(labels);
west = zeros(n,1); north = zeros(n,1); south = zeros(n,1);
main = zeros(n,1); downstream = zeros(n,1); outlet = zeros(n,1);
for i = 1:n
    label = labels{i};
    west(i) = mean(T.(['WestGauge_' label '_m3s'])(window));
    north(i) = mean(T.(['NorthGauge_' label '_m3s'])(window));
    south(i) = mean(T.(['SouthGauge_' label '_m3s'])(window));
    main(i) = mean(T.(['MainGauge_' label '_m3s'])(window));
    downstream(i) = mean(T.(['DownstreamGauge_' label '_m3s'])(window));
    outlet(i) = mean(T.(['Outlet_' label '_m3s'])(window));
end
branch_error = 100.*max(abs([west north south]-Cfg.inflow_m3s),[],2) ...
    ./Cfg.inflow_m3s;
main_error = 100.*(main-3*Cfg.inflow_m3s)./(3*Cfg.inflow_m3s);
outlet_error = 100.*(outlet-3*Cfg.inflow_m3s)./(3*Cfg.inflow_m3s);
Steady = table(repmat(Cfg.case_id,n,1),models,west,north,south,main, ...
    downstream,outlet,branch_error,main_error,outlet_error, ...
    'VariableNames',{'CaseID','Model','WestMean_m3s','NorthMean_m3s', ...
    'SouthMean_m3s','MainMean_m3s','DownstreamMean_m3s','OutletMean_m3s', ...
    'MaxBranchError_pct','MainError_pct','OutletError_pct'});
end

function Cfg = default_config(flow_regime)
Cfg = struct();
Cfg.flow_regime = lower(string(flow_regime));
Cfg.length_m = 990;
Cfg.width_m = 990;
Cfg.fine_dx_m = 10;
Cfg.coarse_dx_m = 30;
Cfg.channel_width_m = 10;
Cfg.channel_depth_m = 1;
Cfg.channel_manning_n = 0.035;
Cfg.floodplain_manning_n = 0.035;
Cfg.channel_slope_mpm = 0.005;
Cfg.valley_side_slope_mpm = 0.005;
Cfg.max_side_relief_m = inf;
Cfg.outlet_elevation_m = 0;
Cfg.duration_min = 180;
Cfg.dt_min = 0.01;
Cfg.record_dt_min = 1;
Cfg.depth_tolerance_m = 1e-6;
Cfg.wet_threshold_m = 0.01;
Cfg.map_thresholds_m = [0.01 0.10];
Cfg.snapshot_times_min = [30 90 150];

% Physical network lines are plotted to the domain edges. Rasterized channel
% endpoints are placed at each grid's own boundary-cell centers.
Cfg.main_vertices = [0 495; 990 495];
Cfg.north_vertices = [495 990; 495 495];
Cfg.south_vertices = [495 0; 495 495];

switch Cfg.flow_regime
    case "overbank"
        Cfg.case_id = "P1-SUBGRID-NEAL-2D-CROSS-001";
        Cfg.total_inflow_m3s = 300;
    case "inbank"
        Cfg.case_id = "P1-SUBGRID-NEAL-2D-CROSS-INBANK-001";
        Cfg.channel_depth_m = 3;
        Cfg.total_inflow_m3s = 50;
    otherwise
        error('Unknown flow regime: %s. Use overbank or inbank.',Cfg.flow_regime);
end
Cfg.inflow_m3s = Cfg.total_inflow_m3s/3;
area_bankfull = Cfg.channel_width_m*Cfg.channel_depth_m;
radius_bankfull = area_bankfull/(Cfg.channel_width_m+2*Cfg.channel_depth_m);
Cfg.bankfull_discharge_m3s = area_bankfull*radius_bankfull^(2/3)* ...
    sqrt(Cfg.channel_slope_mpm)/Cfg.channel_manning_n;
if Cfg.flow_regime == "inbank"
    assert(Cfg.total_inflow_m3s < Cfg.bankfull_discharge_m3s, ...
        'The in-bank forcing must remain below Manning bankfull capacity.');
end
Cfg.q_west = @(t) Cfg.inflow_m3s + 0.*t;
Cfg.q_north = @(t) Cfg.inflow_m3s + 0.*t;
Cfg.q_south = @(t) Cfg.inflow_m3s + 0.*t;

Cfg.west_gauge = [315 495];
Cfg.north_gauge = [495 675];
Cfg.south_gauge = [495 315];
Cfg.main_gauge = [675 495];
Cfg.downstream_gauge = [885 495];
Cfg.confluence = [495 495];
end

function Confinement = score_confinement(Cfg,G,D)
models = ["FineExplicit10m";"CoarseOrdinary30m";"NealSubgrid30m"];
maps = {D.fine,D.coarse,D.neal};
floodplain = ~G.fine.channel_mask;
n = numel(models);
max_depth = zeros(n,1); wet_area = zeros(n,1); confined = false(n,1);
for i = 1:n
    values = maps{i}(floodplain);
    max_depth(i) = max(values,[],'omitnan');
    wet_area(i) = nnz(values>Cfg.wet_threshold_m)*G.fine.area;
    confined(i) = wet_area(i)==0;
end
Confinement = table(repmat(Cfg.case_id,n,1),models, ...
    repmat(Cfg.total_inflow_m3s,n,1),repmat(Cfg.bankfull_discharge_m3s,n,1), ...
    repmat(Cfg.total_inflow_m3s/Cfg.bankfull_discharge_m3s,n,1), ...
    max_depth,wet_area,confined, ...
    'VariableNames',{'CaseID','Model','TotalInflow_m3s', ...
    'ManningBankfullCapacity_m3s','CapacityRatio','MaxFloodplainDepth_m', ...
    'FloodplainWetArea_m2','ConfinedBelow001m'});
end

function G = build_grids(Cfg)
ratio = Cfg.coarse_dx_m / Cfg.fine_dx_m;
assert(ratio == round(ratio), 'Fine-to-coarse ratio must be an integer.');
assert(mod(Cfg.length_m, Cfg.coarse_dx_m) == 0, ...
    'Domain length must be divisible by the coarse resolution.');
assert(mod(Cfg.width_m, Cfg.coarse_dx_m) == 0, ...
    'Domain width must be divisible by the coarse resolution.');

P_fine = network_points(Cfg, Cfg.fine_dx_m);
P_coarse = network_points(Cfg, Cfg.coarse_dx_m);

dx = Cfg.fine_dx_m;
x = (dx/2:dx:Cfg.length_m-dx/2)';
y = (dx/2:dx:Cfg.width_m-dx/2)';
[X,Y] = meshgrid(x,y);
[z_flood,~,~] = valley_surface(X,Y,P_fine,Cfg);
[channel_mask,channel_dout] = rasterize_network(x,y,P_fine);
z_fine = z_flood;
z_fine(channel_mask) = Cfg.outlet_elevation_m + ...
    Cfg.channel_slope_mpm .* channel_dout(channel_mask);
roughness_fine = Cfg.floodplain_manning_n .* ones(size(z_fine));
roughness_fine(channel_mask) = Cfg.channel_manning_n;

G.fine = make_grid("fine", dx, x, y, z_fine, roughness_fine, Cfg);
G.fine.channel_mask = channel_mask;
G.fine.network = P_fine;

dx = Cfg.coarse_dx_m;
xc = (dx/2:dx:Cfg.length_m-dx/2)';
yc = (dx/2:dx:Cfg.width_m-dx/2)';
[Xc,Yc] = meshgrid(xc,yc);
[z_flood_coarse,~,~] = valley_surface(Xc,Yc,P_fine,Cfg);
[channel_mask_coarse,channel_dout_coarse] = rasterize_network(xc,yc,P_coarse);

% Ordinary resampling samples the explicit fine-grid fields at 30 m cell
% centers. Channel cells therefore represent the 10 m channel as a 30 m
% depression, which is the resolution error the subgrid model should reduce.
z_ordinary = interp2(X,Y,z_fine,Xc,Yc,'linear');
n_ordinary = interp2(X,Y,roughness_fine,Xc,Yc,'linear');
G.coarse = make_grid("coarse", dx, xc, yc, z_ordinary, n_ordinary, Cfg);
G.coarse.channel_mask = channel_mask_coarse;
G.coarse.network = P_coarse;

% Neal cells retain the floodplain elevation and embed a rectangular
% channel of the original 10 m width and 1 m bankfull depth.
z_neal = z_flood_coarse;
z_neal(channel_mask_coarse) = Cfg.outlet_elevation_m + ...
    Cfg.channel_slope_mpm .* channel_dout_coarse(channel_mask_coarse) + ...
    Cfg.channel_depth_m;
n_neal = Cfg.floodplain_manning_n .* ones(size(z_neal));
G.neal = make_grid("neal", dx, xc, yc, z_neal, n_neal, Cfg);
G.neal.channel_mask = channel_mask_coarse;
G.neal.network = P_coarse;
G.neal.River_Width(channel_mask_coarse) = Cfg.channel_width_m;
G.neal.River_Depth(channel_mask_coarse) = Cfg.channel_depth_m;
G.neal.n_channel(channel_mask_coarse) = Cfg.channel_manning_n;

models = {'fine','coarse','neal'};
for i = 1:numel(models)
    name = models{i};
    grid = G.(name);
    grid.inlet_rc = [nearest_cell(grid,[grid.x(1) Cfg.confluence(2)]); ...
        nearest_cell(grid,[Cfg.confluence(1) grid.y(end)]); ...
        nearest_cell(grid,[Cfg.confluence(1) grid.y(1)])];
    grid.row_outlet = (1:grid.ny)';
    grid.col_outlet = grid.nx.*ones(grid.ny,1);
    grid.outlet_index = false(size(grid.z));
    grid.outlet_index(:,end) = true;
    assert(any(grid.channel_mask(:,end)), ...
        '%s channel does not reach the eastern boundary.',name);
    G.(name) = grid;
end

assert_connected(G.fine.channel_mask);
assert_connected(G.coarse.channel_mask);
hp2d_neal_validate_geometry(G.neal.River_Width, ...
    G.neal.River_Depth, G.neal.dx);
assert(any(G.neal.channel_mask(:,end)), ...
    'The embedded channel must reach the eastern outlet boundary.');
assert(all(G.neal.outlet_index(:,end)) && nnz(G.neal.outlet_index)==G.neal.ny, ...
    'The complete eastern boundary must use the normal-flow outlet.');
end

function grid = make_grid(mode,dx,x,y,z,roughness,Cfg)
grid = struct();
grid.mode = mode;
grid.dx = dx;
grid.area = dx^2;
grid.x = x(:);
grid.y = y(:);
grid.z = z;
grid.roughness = roughness;
grid.nx = size(z,2);
grid.ny = size(z,1);
grid.n_channel = roughness;
grid.n_flood = Cfg.floodplain_manning_n .* ones(size(z));
grid.River_Width = zeros(size(z));
grid.River_Depth = zeros(size(z));
grid.outlet_index = false(size(z));
grid.outlet_type = 1;
grid.slope_outlet = Cfg.channel_slope_mpm;
grid.row_outlet = [];
grid.col_outlet = [];
grid.inlet_rc = [];
end

function P = network_points(Cfg,spacing)
main_vertices = [spacing/2 Cfg.confluence(2); ...
    Cfg.length_m-spacing/2 Cfg.confluence(2)];
north_vertices = [Cfg.confluence(1) Cfg.width_m-spacing/2; Cfg.confluence];
south_vertices = [Cfg.confluence(1) spacing/2; Cfg.confluence];
[xm,ym,sm] = sample_polyline(main_vertices,spacing);
main_length = sm(end);
dm = main_length - sm;

[~,join_idx] = min(hypot(xm-Cfg.confluence(1),ym-Cfg.confluence(2)));
downstream_length = dm(join_idx);
[xn,yn,sn] = sample_polyline(north_vertices,spacing);
dn = (sn(end) - sn) + downstream_length;
[xs,ys,ss] = sample_polyline(south_vertices,spacing);
ds = (ss(end) - ss) + downstream_length;

xy = [xm ym; xn yn; xs ys];
dout = [dm; dn; ds];
[xy_unique,ia] = unique(xy,'rows','stable');
P.x = xy_unique(:,1);
P.y = xy_unique(:,2);
P.dout = dout(ia);
end

function [x,y,s] = sample_polyline(vertices,spacing)
x = [];
y = [];
s = [];
distance = 0;
for i = 1:size(vertices,1)-1
    delta = vertices(i+1,:) - vertices(i,:);
    length_m = hypot(delta(1),delta(2));
    assert(delta(1) == 0 || delta(2) == 0, ...
        'The Neal benchmark uses orthogonal D4 channel segments.');
    n = length_m / spacing;
    assert(abs(n-round(n)) < 1e-10, ...
        'Every channel segment must align with the grid spacing.');
    f = (0:round(n))' ./ round(n);
    if i > 1; f = f(2:end); end
    x = [x; vertices(i,1) + f.*delta(1)]; %#ok<AGROW>
    y = [y; vertices(i,2) + f.*delta(2)]; %#ok<AGROW>
    s = [s; distance + f.*length_m]; %#ok<AGROW>
    distance = distance + length_m;
end
end

function [z_flood,distance,dout] = valley_surface(X,Y,P,Cfg)
distance2 = inf(size(X));
dout = zeros(size(X));
for i = 1:numel(P.x)
    d2 = (X-P.x(i)).^2 + (Y-P.y(i)).^2;
    update = d2 < distance2;
    distance2(update) = d2(update);
    dout(update) = P.dout(i);
end
distance = sqrt(distance2);
thalweg = Cfg.outlet_elevation_m + Cfg.channel_slope_mpm .* dout;
side_relief = min(Cfg.valley_side_slope_mpm .* distance, ...
    Cfg.max_side_relief_m);
z_flood = thalweg + Cfg.channel_depth_m + side_relief;
end

function [mask,dout] = rasterize_network(x,y,P)
dx = x(2)-x(1);
dy = y(2)-y(1);
mask = false(numel(y),numel(x));
dout = nan(size(mask));
for i = 1:numel(P.x)
    col = round((P.x(i)-x(1))/dx)+1;
    row = round((P.y(i)-y(1))/dy)+1;
    if row >= 1 && row <= numel(y) && col >= 1 && col <= numel(x)
        mask(row,col) = true;
        dout(row,col) = P.dout(i);
    end
end
end

function rc = nearest_cell(grid,xy)
[~,col] = min(abs(grid.x-xy(1)));
[~,row] = min(abs(grid.y-xy(2)));
rc = [row col];
end

function assert_connected(mask)
[row,col] = find(mask,1,'first');
visited = false(size(mask));
queue = zeros(nnz(mask),2);
queue(1,:) = [row col];
visited(row,col) = true;
head = 1;
tail = 1;
while head <= tail
    rc = queue(head,:);
    head = head + 1;
    neighbors = [rc(1)-1 rc(2); rc(1)+1 rc(2); ...
        rc(1) rc(2)-1; rc(1) rc(2)+1];
    for j = 1:4
        r = neighbors(j,1); c = neighbors(j,2);
        if r >= 1 && r <= size(mask,1) && c >= 1 && c <= size(mask,2) && ...
                mask(r,c) && ~visited(r,c)
            tail = tail + 1;
            queue(tail,:) = [r c];
            visited(r,c) = true;
        end
    end
end
assert(nnz(visited) == nnz(mask), 'The channel network is not D4 connected.');
end

function [T,Mass,MaxEta,Snapshots,Final] = run_models(Cfg,G)
dt_s = Cfg.dt_min*60;
nsteps = round(Cfg.duration_min/Cfg.dt_min);
record_every = round(Cfg.record_dt_min/Cfg.dt_min);
nrecords = floor(nsteps/record_every)+1;

models = {'fine','coarse','neal'};
for i = 1:3
    S.(models{i}) = initialize_state(G.(models{i}));
end
initial_storage = zeros(3,1);
cum_out = zeros(3,1);
for i = 1:3
    initial_storage(i) = compute_storage(S.(models{i}),G.(models{i}));
end
cum_input = 0;

R = initialize_records(nrecords);
MaxEta.fine = -inf(size(G.fine.z));
MaxEta.coarse = -inf(size(G.coarse.z));
MaxEta.neal = -inf(size(G.neal.z));
Snapshots.time_min = Cfg.snapshot_times_min(:)';
for i = 1:3
    Snapshots.(models{i}) = nan(G.fine.ny,G.fine.nx, ...
        numel(Cfg.snapshot_times_min));
end

[R,MaxEta,Snapshots] = store_record(R,MaxEta,Snapshots,1,0,S,G,Cfg);
rec = 1;

for it = 1:nsteps
    time_min = (it-1)*Cfg.dt_min;
    q_west = Cfg.q_west(time_min);
    q_north = Cfg.q_north(time_min);
    q_south = Cfg.q_south(time_min);
    for i = 1:3
        name = models{i};
        S.(name) = apply_inflows(S.(name),G.(name), ...
            [q_west q_north q_south],dt_s);
        S.(name) = step_model(S.(name),G.(name),Cfg.dt_min);
        qout = outlet_discharge(S.(name),G.(name));
        cum_out(i) = cum_out(i) + qout*dt_s;
    end
    cum_input = cum_input + (q_west+q_north+q_south)*dt_s;

    if mod(it,record_every) == 0 || it == nsteps
        rec = rec+1;
        [R,MaxEta,Snapshots] = store_record(R,MaxEta,Snapshots,rec, ...
            it*Cfg.dt_min,S,G,Cfg);
    end
end

T = struct2table(R);
final_storage = zeros(3,1);
for i = 1:3
    final_storage(i) = compute_storage(S.(models{i}),G.(models{i}));
end
residual = initial_storage + cum_input - cum_out - final_storage;
Mass = table(repmat(Cfg.case_id,3,1), ...
    ["FineExplicit10m";"CoarseOrdinary30m";"NealSubgrid30m"], ...
    repmat(cum_input,3,1),initial_storage,cum_out,final_storage,residual, ...
    100.*residual./max(cum_input,1), ...
    'VariableNames',{'CaseID','Model','InputVolume_m3','InitialStorage_m3', ...
    'OutletVolume_m3','FinalStorage_m3','Residual_m3','MassResidual_pct'});
Final = S;
end

function R = initialize_records(n)
names = {'Time_min','WestInflow_m3s','NorthInflow_m3s', ...
    'SouthInflow_m3s','TotalInflow_m3s', ...
    'Outlet_Fine_m3s','Outlet_Coarse_m3s','Outlet_Neal_m3s', ...
    'WestGauge_Fine_m3s','WestGauge_Coarse_m3s','WestGauge_Neal_m3s', ...
    'NorthGauge_Fine_m3s','NorthGauge_Coarse_m3s','NorthGauge_Neal_m3s', ...
    'SouthGauge_Fine_m3s','SouthGauge_Coarse_m3s','SouthGauge_Neal_m3s', ...
    'MainGauge_Fine_m3s','MainGauge_Coarse_m3s','MainGauge_Neal_m3s', ...
    'DownstreamGauge_Fine_m3s','DownstreamGauge_Coarse_m3s','DownstreamGauge_Neal_m3s', ...
    'Storage_Fine_m3','Storage_Coarse_m3','Storage_Neal_m3', ...
    'WetArea_Fine_m2','WetArea_Coarse_m2','WetArea_Neal_m2', ...
    'ConfluenceDepth_Fine_m','ConfluenceDepth_Coarse_m','ConfluenceDepth_Neal_m'};
for i = 1:numel(names)
    R.(names{i}) = zeros(n,1);
end
end

function state = initialize_state(grid)
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
        grid.River_Width,grid.River_Depth,grid.dx);
else
    state.C_a = grid.area.*ones(grid.ny,grid.nx);
end
end

function state = apply_inflows(state,grid,q_m3s,dt_s)
for i = 1:numel(q_m3s)
    if q_m3s(i) <= 0; continue; end
    row = grid.inlet_rc(i,1);
    col = grid.inlet_rc(i,2);
    volume = q_m3s(i)*dt_s;
    if grid.mode == "neal"
        h = state.d_tot(row,col)/1000;
        h = hp2d_neal_apply_volume_change(h,volume, ...
            grid.River_Width(row,col),grid.River_Depth(row,col),grid.dx);
        state.d_tot(row,col) = 1000*h;
    else
        state.d_tot(row,col) = state.d_tot(row,col) + 1000*volume/grid.area;
    end
end
end

function state = step_model(state,grid,dt_min)
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

function [R,MaxEta,Snapshots] = store_record(R,MaxEta,Snapshots,i,time,S,G,Cfg)
R.Time_min(i) = time;
R.WestInflow_m3s(i) = Cfg.q_west(time);
R.NorthInflow_m3s(i) = Cfg.q_north(time);
R.SouthInflow_m3s(i) = Cfg.q_south(time);
R.TotalInflow_m3s(i) = R.WestInflow_m3s(i)+R.NorthInflow_m3s(i)+ ...
    R.SouthInflow_m3s(i);

models = {'fine','coarse','neal'};
labels = {'Fine','Coarse','Neal'};
for k = 1:3
    name = models{k}; label = labels{k};
    grid = G.(name); state = S.(name);
    R.(['Outlet_' label '_m3s'])(i) = outlet_discharge(state,grid);
    R.(['WestGauge_' label '_m3s'])(i) = section_discharge_x( ...
        state,grid,Cfg.west_gauge(1));
    R.(['NorthGauge_' label '_m3s'])(i) = section_discharge_y( ...
        state,grid,Cfg.north_gauge(2),1);
    R.(['SouthGauge_' label '_m3s'])(i) = section_discharge_y( ...
        state,grid,Cfg.south_gauge(2),-1);
    R.(['MainGauge_' label '_m3s'])(i) = section_discharge_x( ...
        state,grid,Cfg.main_gauge(1));
    R.(['DownstreamGauge_' label '_m3s'])(i) = section_discharge_x( ...
        state,grid,Cfg.downstream_gauge(1));
    R.(['Storage_' label '_m3'])(i) = compute_storage(state,grid);
    depth_fine = project_depth_to_fine(state,grid,G.fine,Cfg);
    R.(['WetArea_' label '_m2'])(i) = nnz(depth_fine>Cfg.wet_threshold_m)*G.fine.area;
    R.(['ConfluenceDepth_' label '_m'])(i) = channel_depth_at( ...
        state,grid,Cfg.confluence);

    eta = model_wse(state,grid);
    wet = state.d_tot > Cfg.depth_tolerance_m*1000;
    MaxEta.(name)(wet) = max(MaxEta.(name)(wet),eta(wet));

    snap_idx = find(abs(Snapshots.time_min-time)<1e-9,1);
    if ~isempty(snap_idx)
        Snapshots.(name)(:,:,snap_idx) = depth_fine;
    end
end
end

function q = outlet_discharge(state,grid)
q = sum(state.outlet_flow(:),'omitnan')*grid.area/1000/3600;
end

function q = section_discharge_x(state,grid,x_location)
[~,col] = min(abs(grid.x-x_location));
col = min(col,grid.nx-1);
v = state.outflow(:,col,1);
q = sum(v(:),'omitnan')*grid.area/1000/3600;
end

function q = section_discharge_y(state,grid,y_location,direction_sign)
[~,row] = min(abs(grid.y-y_location));
row = min(row,grid.ny-1);
v = state.outflow(row,:,2);
q = direction_sign*sum(v(:),'omitnan')*grid.area/1000/3600;
end

function d = channel_depth_at(state,grid,xy)
rc = nearest_cell(grid,xy);
d = max(state.d_tot(rc(1),rc(2))/1000,0);
end

function storage = compute_storage(state,grid)
h = max(state.d_tot/1000,0);
if grid.mode == "neal"
    V = hp2d_neal_cell_volume(h,grid.River_Width,grid.River_Depth,grid.dx);
else
    V = h.*grid.area;
end
storage = sum(V(:),'omitnan');
end

function eta = model_wse(state,grid)
h = max(state.d_tot/1000,0);
if grid.mode == "neal"
    eta = grid.z-grid.River_Depth+h;
else
    eta = grid.z+h;
end
end

function depth = project_depth_to_fine(state,grid,fine,Cfg)
if grid.mode == "fine"
    depth = max(state.d_tot/1000,0);
    return;
end
ratio = grid.dx/fine.dx;
assert(ratio == round(ratio),'Projection ratio must be an integer.');
eta = model_wse(state,grid);
wet = state.d_tot > Cfg.depth_tolerance_m*1000;
eta_fine = repelem(eta,ratio,ratio);
wet_fine = repelem(wet,ratio,ratio);
depth = max(eta_fine-fine.z,0);
depth(~wet_fine) = 0;
end

function Metrics = score_hydrographs(Cfg,T,Mass)
models = ["FineExplicit10m";"CoarseOrdinary30m";"NealSubgrid30m"];
out = {T.Outlet_Fine_m3s,T.Outlet_Coarse_m3s,T.Outlet_Neal_m3s};
west = {T.WestGauge_Fine_m3s,T.WestGauge_Coarse_m3s,T.WestGauge_Neal_m3s};
north = {T.NorthGauge_Fine_m3s,T.NorthGauge_Coarse_m3s,T.NorthGauge_Neal_m3s};
south = {T.SouthGauge_Fine_m3s,T.SouthGauge_Coarse_m3s,T.SouthGauge_Neal_m3s};
main = {T.MainGauge_Fine_m3s,T.MainGauge_Coarse_m3s,T.MainGauge_Neal_m3s};
down = {T.DownstreamGauge_Fine_m3s,T.DownstreamGauge_Coarse_m3s,T.DownstreamGauge_Neal_m3s};
storage = {T.Storage_Fine_m3,T.Storage_Coarse_m3,T.Storage_Neal_m3};
confluence = {T.ConfluenceDepth_Fine_m,T.ConfluenceDepth_Coarse_m,T.ConfluenceDepth_Neal_m};

n = 3;
out_rmse = zeros(n,1); out_nse = ones(n,1);
west_rmse = zeros(n,1); west_nse = ones(n,1);
north_rmse = zeros(n,1); north_nse = ones(n,1);
south_rmse = zeros(n,1); south_nse = ones(n,1);
main_rmse = zeros(n,1); main_nse = ones(n,1);
down_rmse = zeros(n,1); down_nse = ones(n,1);
plateau_error = zeros(n,1); t95_error = zeros(n,1);
volume_error = zeros(n,1); storage_rmse = zeros(n,1);
peak_storage_error = zeros(n,1); wet_area_rmse = zeros(n,1);
peak_wet_area_error = zeros(n,1); confluence_rmse = zeros(n,1);
ref_volume = trapz(T.Time_min*60,out{1});
steady_window = T.Time_min >= Cfg.duration_min-20;
ref_plateau = mean(out{1}(steady_window));
ref_t95 = first_crossing_time(T.Time_min,out{1},0.95*ref_plateau);
ref_peak_storage = max(storage{1});
wet_area = {T.WetArea_Fine_m2,T.WetArea_Coarse_m2,T.WetArea_Neal_m2};
ref_peak_wet_area = max(wet_area{1});
for i = 1:n
    out_rmse(i) = rmse(out{1},out{i});
    out_nse(i) = nse(out{1},out{i});
    west_rmse(i) = rmse(west{1},west{i});
    west_nse(i) = nse(west{1},west{i});
    north_rmse(i) = rmse(north{1},north{i});
    north_nse(i) = nse(north{1},north{i});
    south_rmse(i) = rmse(south{1},south{i});
    south_nse(i) = nse(south{1},south{i});
    main_rmse(i) = rmse(main{1},main{i});
    main_nse(i) = nse(main{1},main{i});
    down_rmse(i) = rmse(down{1},down{i});
    down_nse(i) = nse(down{1},down{i});
    plateau = mean(out{i}(steady_window));
    plateau_error(i) = 100*(plateau-ref_plateau)/max(ref_plateau,eps);
    t95_error(i) = first_crossing_time( ...
        T.Time_min,out{i},0.95*ref_plateau)-ref_t95;
    volume_error(i) = 100*(trapz(T.Time_min*60,out{i})-ref_volume)/max(ref_volume,eps);
    storage_rmse(i) = rmse(storage{1},storage{i});
    peak_storage_error(i) = 100*(max(storage{i})-ref_peak_storage)/max(ref_peak_storage,eps);
    wet_area_rmse(i) = rmse(wet_area{1},wet_area{i});
    peak_wet_area_error(i) = 100*(max(wet_area{i})-ref_peak_wet_area)/max(ref_peak_wet_area,eps);
    confluence_rmse(i) = rmse(confluence{1},confluence{i});
end

Metrics = table(repmat(Cfg.case_id,n,1),models,out_rmse,out_nse, ...
    west_rmse,west_nse,north_rmse,north_nse,south_rmse,south_nse, ...
    main_rmse,main_nse,down_rmse,down_nse,plateau_error,t95_error, ...
    volume_error,storage_rmse,peak_storage_error,wet_area_rmse, ...
    peak_wet_area_error,confluence_rmse,Mass.MassResidual_pct, ...
    'VariableNames',{'CaseID','Model','Outlet_RMSE_m3s','Outlet_NSE', ...
    'WestGauge_RMSE_m3s','WestGauge_NSE','NorthGauge_RMSE_m3s', ...
    'NorthGauge_NSE','SouthGauge_RMSE_m3s','SouthGauge_NSE', ...
    'MainGauge_RMSE_m3s','MainGauge_NSE','DownstreamGauge_RMSE_m3s', ...
    'DownstreamGauge_NSE','OutletPlateauError_pct','TimeTo95Error_min', ...
    'OutletVolumeError_pct','Storage_RMSE_m3','PeakStorageError_pct', ...
    'WetArea_RMSE_m2','PeakWetAreaError_pct','ConfluenceDepth_RMSE_m', ...
    'MassResidual_pct'});
end

function time = first_crossing_time(t,q,threshold)
idx = find(q>=threshold,1,'first');
if isempty(idx)
    time = t(end);
else
    time = t(idx);
end
end

function [MapMetrics,MaxDepth] = score_maximum_depth(Cfg,G,MaxEta)
MaxDepth.fine = max(MaxEta.fine-G.fine.z,0);
MaxDepth.fine(~isfinite(MaxDepth.fine)) = 0;
MaxDepth.coarse = project_max_eta(MaxEta.coarse,G.coarse,G.fine);
MaxDepth.neal = project_max_eta(MaxEta.neal,G.neal,G.fine);

models = ["FineExplicit10m";"CoarseOrdinary30m";"NealSubgrid30m"];
maps = {MaxDepth.fine,MaxDepth.coarse,MaxDepth.neal};
n = 3;
depth_rmse = zeros(n,1); mae = zeros(n,1); bias = zeros(n,1);
max_error = zeros(n,1); wet_error_001 = zeros(n,1); wet_error_010 = zeros(n,1);
csi_001 = ones(n,1); csi_010 = ones(n,1);
for i = 1:n
    delta = maps{i}-maps{1};
    depth_rmse(i) = sqrt(mean(delta(:).^2));
    mae(i) = mean(abs(delta(:)));
    bias(i) = mean(delta(:));
    max_error(i) = max(abs(delta(:)));
    [wet_error_001(i),csi_001(i)] = wet_scores( ...
        maps{1},maps{i},Cfg.map_thresholds_m(1));
    [wet_error_010(i),csi_010(i)] = wet_scores( ...
        maps{1},maps{i},Cfg.map_thresholds_m(2));
end
MapMetrics = table(repmat(Cfg.case_id,n,1),models,depth_rmse,mae,bias, ...
    max_error,wet_error_001,csi_001,wet_error_010,csi_010, ...
    'VariableNames',{'CaseID','Model','MaxDepth_RMSE_m','MaxDepth_MAE_m', ...
    'MaxDepth_Bias_m','MaxDepth_MaxError_m','WetAreaError_001_pct', ...
    'CSI_001','WetAreaError_010_pct','CSI_010'});
end

function depth = project_max_eta(max_eta,grid,fine)
ratio = grid.dx/fine.dx;
wet = isfinite(max_eta);
eta_fine = repelem(max_eta,ratio,ratio);
wet_fine = repelem(wet,ratio,ratio);
depth = max(eta_fine-fine.z,0);
depth(~wet_fine) = 0;
depth(~isfinite(depth)) = 0;
end

function [area_error,csi] = wet_scores(reference,simulation,threshold)
ref_wet = reference>threshold;
sim_wet = simulation>threshold;
area_error = 100*(nnz(sim_wet)-nnz(ref_wet))/max(nnz(ref_wet),1);
csi = nnz(ref_wet & sim_wet)/max(nnz(ref_wet | sim_wet),1);
end

function Pass = make_pass_table(Cfg,Metrics,MapMetrics,Mass,Steady,Confinement)
coarse = Metrics(Metrics.Model=="CoarseOrdinary30m",:);
neal = Metrics(Metrics.Model=="NealSubgrid30m",:);
coarse_map = MapMetrics(MapMetrics.Model=="CoarseOrdinary30m",:);
neal_map = MapMetrics(MapMetrics.Model=="NealSubgrid30m",:);
neal_steady = Steady(Steady.Model=="NealSubgrid30m",:);

outlet_better = neal.Outlet_RMSE_m3s < coarse.Outlet_RMSE_m3s && ...
    neal.Outlet_NSE > coarse.Outlet_NSE;
main_better = neal.MainGauge_RMSE_m3s < coarse.MainGauge_RMSE_m3s && ...
    neal.MainGauge_NSE > coarse.MainGauge_NSE;
branch_better = all([neal.WestGauge_RMSE_m3s neal.NorthGauge_RMSE_m3s ...
    neal.SouthGauge_RMSE_m3s] < [coarse.WestGauge_RMSE_m3s ...
    coarse.NorthGauge_RMSE_m3s coarse.SouthGauge_RMSE_m3s]);
map_better = neal_map.MaxDepth_RMSE_m < coarse_map.MaxDepth_RMSE_m;
if Cfg.flow_regime == "inbank"
    extent_better = neal_map.CSI_010 >= coarse_map.CSI_010;
else
    extent_better = neal_map.CSI_010 > coarse_map.CSI_010;
end
storage_better = neal.Storage_RMSE_m3 < coarse.Storage_RMSE_m3;
absolute_hydrograph = neal.Outlet_NSE > 0.90;
wet_area_magnitude = abs(neal_map.WetAreaError_010_pct) < 10;
absolute_map = neal_map.MaxDepth_RMSE_m < 0.15 && ...
    neal_map.CSI_010 > 0.80 && wet_area_magnitude;
mass_pass = all(abs(Mass.MassResidual_pct)<0.1);
steady_pass = neal_steady.MaxBranchError_pct < 2 && ...
    abs(neal_steady.OutletError_pct) < 2;
if Cfg.flow_regime == "inbank"
    neal_confinement = Confinement(Confinement.Model=="NealSubgrid30m",:);
    confinement_pass = Cfg.total_inflow_m3s < Cfg.bankfull_discharge_m3s && ...
        neal_confinement.ConfinedBelow001m;
else
    confinement_pass = true;
end
pass = outlet_better && branch_better && main_better && map_better && extent_better && ...
    storage_better && absolute_hydrograph && absolute_map && mass_pass && ...
    steady_pass && confinement_pass;

Pass = table(Cfg.case_id,outlet_better,branch_better,main_better, ...
    map_better,extent_better,storage_better,absolute_hydrograph, ...
    wet_area_magnitude,absolute_map,mass_pass,steady_pass,confinement_pass,pass, ...
    'VariableNames',{'CaseID','OutletHydrographImproved','AllBranchGaugesImproved', ...
    'MainGaugeHydrographImproved','MaximumDepthImproved', ...
    'WetExtentCriterionPassed','StorageRMSEImproved','OutletNSEAbove090', ...
    'WetAreaErrorWithin10pct','MapThresholdsPassed', ...
    'MassResidualUnder01pct','SteadyFlowErrorsUnder2pct', ...
    'InbankConfinementPassed','Pass'});
end

function plot_geometry(Cfg,G,fig_dir)
fig = figure('Color','w','Position',[60 60 1500 500]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
grids = {G.fine,G.coarse,G.neal};
titles = {'(a) Fine explicit terrain','(b) Ordinary 30 m terrain', ...
    '(c) Neal (2012) 30 m terrain'};
for i = 1:3
    ax = nexttile;
    imagesc(ax,grids{i}.x,grids{i}.y,grids{i}.z);
    axis(ax,'xy','equal','tight');
    hold(ax,'on');
    plot_network(ax,Cfg,'k-',1.4);
    if i == 1
        plot(ax,Cfg.main_vertices(1,1),Cfg.main_vertices(1,2),'wo', ...
            'MarkerFaceColor',[0.00 0.45 0.74],'MarkerSize',6);
        plot(ax,Cfg.north_vertices(1,1),Cfg.north_vertices(1,2),'wo', ...
            'MarkerFaceColor',[0.85 0.33 0.10],'MarkerSize',6);
        plot(ax,Cfg.south_vertices(1,1),Cfg.south_vertices(1,2),'wo', ...
            'MarkerFaceColor',[0.49 0.18 0.56],'MarkerSize',6);
        plot(ax,Cfg.main_vertices(end,1),Cfg.main_vertices(end,2),'ks', ...
            'MarkerFaceColor','w','MarkerSize',6);
    end
    title(ax,titles{i},'FontWeight','normal');
    xlabel(ax,'Easting (m)');
    if i==1; ylabel(ax,'Northing (m)'); else; yticklabels(ax,[]); end
    colormap(ax,turbo(256));
    cb = colorbar(ax); cb.Label.String = 'Elevation (m)';
    set(ax,'FontName','Arial','FontSize',10,'Box','on','Layer','top');
end
export_figure(fig,fig_dir,'Domain_Geometry');
end

function plot_network(ax,Cfg,style,width)
plot(ax,Cfg.north_vertices(:,1),Cfg.north_vertices(:,2),style,'LineWidth',width);
plot(ax,Cfg.south_vertices(:,1),Cfg.south_vertices(:,2),style,'LineWidth',width);
plot(ax,Cfg.main_vertices(:,1),Cfg.main_vertices(:,2),style,'LineWidth',width);
end

function plot_hydrographs(T,fig_dir)
c = model_colors();
fig = figure('Color','w','Position',[50 40 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

ax = nexttile;
plot_series(ax,T.Time_min,T.WestGauge_Fine_m3s,T.WestGauge_Coarse_m3s, ...
    T.WestGauge_Neal_m3s,c);
plot(ax,T.Time_min,T.WestInflow_m3s,'--','Color',[0.2 0.2 0.2],'LineWidth',1.1);
title(ax,'(a) West inlet branch','FontWeight','normal');
ylabel(ax,'Discharge (m^3 s^{-1})');

ax = nexttile;
plot_series(ax,T.Time_min,T.NorthGauge_Fine_m3s,T.NorthGauge_Coarse_m3s, ...
    T.NorthGauge_Neal_m3s,c);
plot(ax,T.Time_min,T.NorthInflow_m3s,'--','Color',[0.2 0.2 0.2],'LineWidth',1.1);
title(ax,'(b) North inlet branch','FontWeight','normal');

ax = nexttile;
plot_series(ax,T.Time_min,T.SouthGauge_Fine_m3s,T.SouthGauge_Coarse_m3s, ...
    T.SouthGauge_Neal_m3s,c);
plot(ax,T.Time_min,T.SouthInflow_m3s,'--','Color',[0.2 0.2 0.2],'LineWidth',1.1);
title(ax,'(c) South inlet branch','FontWeight','normal');

ax = nexttile;
plot_series(ax,T.Time_min,T.MainGauge_Fine_m3s,T.MainGauge_Coarse_m3s, ...
    T.MainGauge_Neal_m3s,c);
title(ax,'(d) Main channel below junction','FontWeight','normal');
xlabel(ax,'Time (min)'); ylabel(ax,'Discharge (m^3 s^{-1})');

ax = nexttile;
plot_series(ax,T.Time_min,T.DownstreamGauge_Fine_m3s, ...
    T.DownstreamGauge_Coarse_m3s,T.DownstreamGauge_Neal_m3s,c);
title(ax,'(e) Downstream internal gauge','FontWeight','normal');
xlabel(ax,'Time (min)');

ax = nexttile;
plot_series(ax,T.Time_min,T.Outlet_Fine_m3s,T.Outlet_Coarse_m3s, ...
    T.Outlet_Neal_m3s,c);
title(ax,'(f) East outlet','FontWeight','normal');
xlabel(ax,'Time (min)');
legend(ax,{'Fine explicit 10 m','Ordinary 30 m','Neal (2012) 30 m'}, ...
    'Location','best','Box','off');
export_figure(fig,fig_dir,'Hydrograph_Comparison');
end

function plot_series(ax,t,ref,coarse,neal,c)
hold(ax,'on');
plot(ax,t,ref,'Color',c.ref,'LineWidth',2.0);
plot(ax,t,coarse,'Color',c.coarse,'LineWidth',1.5);
plot(ax,t,neal,'Color',c.neal,'LineWidth',1.8);
grid(ax,'on'); box(ax,'on'); xlim(ax,[t(1) t(end)]);
set(ax,'FontName','Arial','FontSize',10,'Layer','top');
end

function plot_storage_wet_area(T,fig_dir)
c = model_colors();
fig = figure('Color','w','Position',[100 100 1200 470]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile;
plot_series(ax,T.Time_min,T.Storage_Fine_m3,T.Storage_Coarse_m3, ...
    T.Storage_Neal_m3,c);
xlabel(ax,'Time (min)'); ylabel(ax,'Surface-water storage (m^3)');
title(ax,'(a) Domain storage','FontWeight','normal');
ax = nexttile;
plot_series(ax,T.Time_min,T.WetArea_Fine_m2/1e6,T.WetArea_Coarse_m2/1e6, ...
    T.WetArea_Neal_m2/1e6,c);
xlabel(ax,'Time (min)'); ylabel(ax,'Wet area (km^2)');
title(ax,'(b) Area with depth > 0.01 m','FontWeight','normal');
legend(ax,{'Fine explicit 10 m','Ordinary 30 m','Neal (2012) 30 m'}, ...
    'Location','best','Box','off');
export_figure(fig,fig_dir,'Storage_WetArea');
end

function plot_maximum_depth(Cfg,G,D,fig_dir)
fig = figure('Color','w','Position',[40 40 1450 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
maps = {D.fine,D.coarse,D.neal};
titles = {'(a) Fine explicit 10 m','(b) Ordinary 30 m projected to 10 m', ...
    '(c) Neal (2012) 30 m projected to 10 m'};
dmax = max(D.fine(:));
for i = 1:3
    ax = nexttile;
    imagesc(ax,G.fine.x,G.fine.y,maps{i}); axis(ax,'xy','equal','tight');
    clim(ax,[0 max(dmax,0.01)]); colormap(ax,parula(256));
    title(ax,titles{i},'FontWeight','normal');
    xlabel(ax,'Easting (m)'); if i==1; ylabel(ax,'Northing (m)'); end
    cb = colorbar(ax); cb.Label.String = 'Maximum depth (m)';
    set(ax,'FontName','Arial','FontSize',9,'Box','on');
end

errors = {D.coarse-D.fine,D.neal-D.fine};
error_titles = {'(d) Ordinary minus reference','(e) Neal (2012) minus reference'};
emax = max(abs([errors{1}(:);errors{2}(:)]));
for i = 1:2
    ax = nexttile;
    imagesc(ax,G.fine.x,G.fine.y,errors{i}); axis(ax,'xy','equal','tight');
    clim(ax,[-max(emax,0.01) max(emax,0.01)]); colormap(ax,bluewhitered());
    title(ax,error_titles{i},'FontWeight','normal');
    xlabel(ax,'Easting (m)'); if i==1; ylabel(ax,'Northing (m)'); end
    cb = colorbar(ax); cb.Label.String = 'Depth error (m)';
    set(ax,'FontName','Arial','FontSize',9,'Box','on');
end

ax = nexttile;
ref_wet = D.fine>Cfg.map_thresholds_m(2);
neal_wet = D.neal>Cfg.map_thresholds_m(2);
classification = zeros(size(ref_wet));
classification(ref_wet & ~neal_wet) = -1;
classification(~ref_wet & neal_wet) = 1;
imagesc(ax,G.fine.x,G.fine.y,classification); axis(ax,'xy','equal','tight');
clim(ax,[-1 1]); colormap(ax,[0.80 0.20 0.20; 0.88 0.88 0.88; 0.15 0.45 0.75]);
title(ax,'(f) Neal (2012) wet-extent error at 0.10 m','FontWeight','normal');
xlabel(ax,'Easting (m)');
colorbar(ax,'Ticks',[-1 0 1], ...
    'TickLabels',{'Missed','Agreement','Extra'});
set(ax,'FontName','Arial','FontSize',9,'Box','on');
export_figure(fig,fig_dir,'Maximum_Depth_Comparison');
end

function plot_snapshots(~,G,S,fig_dir)
fig = figure('Color','w','Position',[30 30 1350 1000]);
tiledlayout(3,numel(S.time_min),'TileSpacing','compact','Padding','compact');
models = {'fine','coarse','neal'};
row_labels = {'Fine explicit 10 m','Ordinary 30 m','Neal (2012) 30 m'};
dmax = 0;
for i = 1:3; dmax = max(dmax,max(S.(models{i})(:),[],'omitnan')); end
for i = 1:3
    for j = 1:numel(S.time_min)
        ax = nexttile;
        imagesc(ax,G.fine.x,G.fine.y,S.(models{i})(:,:,j));
        axis(ax,'xy','equal','tight'); clim(ax,[0 max(dmax,0.01)]);
        colormap(ax,parula(256));
        if i==1
            title(ax,sprintf('(%c) %d min',char('a'+j-1),S.time_min(j)), ...
                'FontWeight','normal');
        end
        if j==1; ylabel(ax,row_labels{i}); else; yticklabels(ax,[]); end
        if i==3; xlabel(ax,'Easting (m)'); else; xticklabels(ax,[]); end
        set(ax,'FontName','Arial','FontSize',9,'Box','on');
    end
end
cb = colorbar(nexttile(3),'Location','eastoutside');
cb.Label.String = 'Water depth (m)';
export_figure(fig,fig_dir,'Depth_Snapshots');
end

function export_figure(fig,fig_dir,name)
exportgraphics(fig,fullfile(fig_dir,[name '.png']),'Resolution',300);
exportgraphics(fig,fullfile(fig_dir,[name '.pdf']),'ContentType','vector');
close(fig);
end

function value = rmse(reference,simulation)
delta = simulation-reference;
value = sqrt(mean(delta.^2,'omitnan'));
end

function value = nse(reference,simulation)
den = sum((reference-mean(reference,'omitnan')).^2,'omitnan');
if den <= eps
    value = double(all(abs(reference-simulation)<1e-12));
else
    value = 1-sum((simulation-reference).^2,'omitnan')/den;
end
end

function c = model_colors()
c.ref = [0.05 0.05 0.05];
c.coarse = [0.82 0.26 0.18];
c.neal = [0.05 0.48 0.28];
end

function map = bluewhitered()
n = 256;
lower = [linspace(0.12,1,n/2)' linspace(0.35,1,n/2)' ones(n/2,1)];
upper = [ones(n/2,1) linspace(1,0.20,n/2)' linspace(1,0.15,n/2)'];
map = [lower;upper];
end
