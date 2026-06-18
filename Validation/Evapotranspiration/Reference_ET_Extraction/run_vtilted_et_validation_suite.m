clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

registry_path = fullfile(case_dir, 'ET_Case_Registry.csv');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
cell_dir = fullfile(out_dir, 'CellResults');
fig_dir = fullfile(out_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(cell_dir, 'dir'); mkdir(cell_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

Cases = readtable(registry_path, 'TextType', 'string');
Diagnostics = table();
PassFail = table();

for i = 1:height(Cases)
    Cfg = Cases(i,:);

    if Cfg.regime == "potential_et_formula"
        [DiagRow, CellRows, passed] = run_potential_et_case(Cfg);
    else
        [DiagRow, CellRows, passed] = run_extraction_case(Cfg, functions_dir);
    end

    Diagnostics = [Diagnostics; DiagRow]; %#ok<AGROW>
    writetable(CellRows, fullfile(cell_dir, Cfg.case_id + "_cells.csv"));

    status = "fail";
    if passed
        status = "pass";
    end

    PassFail = [PassFail; table( ...
        Cfg.case_id, Cfg.case_name, Cfg.regime, status, passed, ...
        'VariableNames', {'case_id','case_name','regime','status','report_ready'})]; %#ok<AGROW>
end

writetable(Diagnostics, fullfile(out_dir, 'VTilted_ET_Diagnostics.csv'));
writetable(PassFail, fullfile(out_dir, 'VTilted_ET_Pass_Fail.csv'));
make_et_figures(Diagnostics, cell_dir, fig_dir);

disp(Diagnostics);
disp(PassFail);

function [DiagRow, CellRows, passed] = run_potential_et_case(Cfg)

[DEM, zone_id, zone_name, idx_nan, albedo] = make_vtilted_domain();

Temp = [24 25 24; 23 24 23; 22 23 22];
Temp_max = Temp + 7;
Temp_min = Temp - 8;
U_2 = 2.0 + 0.2 * double(zone_id == 2);
U_R = 60 * ones(size(DEM));
G = zeros(size(DEM));
lat = 39.0;
Krs = 0.16;
day_year = 180;

[ETP_modelled, Ep_modelled] = Evapotranspiration( ...
    DEM, Temp, Temp_max, Temp_min, day_year, lat, U_2, U_R, Krs, albedo, G);
[ETP_ref, Ep_ref] = reference_evapotranspiration( ...
    DEM, Temp, Temp_max, Temp_min, day_year, lat, U_2, Krs, albedo, G);

etp_error = ETP_modelled - ETP_ref;
ep_error = Ep_modelled - Ep_ref;
max_etp_error = max_abs_omitnan(etp_error);
max_ep_error = max_abs_omitnan(ep_error);
nan_mismatch = count_nan_mismatch(ETP_modelled, ETP_ref) + ...
    count_nan_mismatch(Ep_modelled, Ep_ref);

passed = max_etp_error < 1e-10 && max_ep_error < 1e-10 && nan_mismatch == 0;

CellRows = make_cell_rows(Cfg, zone_id, zone_name, idx_nan, ...
    ETP_ref, ETP_modelled, Ep_ref, Ep_modelled, ...
    nan(size(DEM)), nan(size(DEM)), nan(size(DEM)), nan(size(DEM)), ...
    nan(size(DEM)), nan(size(DEM)), nan(size(DEM)), nan(size(DEM)));

DiagRow = table( ...
    Cfg.case_id, Cfg.case_name, Cfg.regime, ...
    mean_omitnan(ETP_ref), mean_omitnan(ETP_modelled), ...
    nan, nan, nan, nan, nan, nan, nan, nan, ...
    max_etp_error, max_ep_error, nan, nan, nan, 0, 0, nan_mismatch, passed, ...
    'VariableNames', diagnostic_names());
end

function [DiagRow, CellRows, passed] = run_extraction_case(Cfg, functions_dir)

[elevation, zone_id, zone_name, idx_nan, ~] = make_vtilted_domain();
dt_days = Cfg.dt_h / 24;
time_step = Cfg.dt_h * 60; %#ok<NASGU>
k = 1; %#ok<NASGU>
C_a = 400 * ones(size(elevation));

flags = struct();
flags.flag_ETP = 1;
flags.flag_input_ETP_map = double(Cfg.mode == "map"); %#ok<STRNU>

LULC_Properties = struct();
LULC_Properties.idx_imp = false(size(elevation));

depths = struct();
depths.d_t = Cfg.ponded_depth_mm * ones(size(elevation));

if Cfg.regime == "internal_etp_mode"
    depths.d_t(:) = 0;
    depths.d_t(zone_id == 2) = Cfg.ponded_depth_mm;
end

if Cfg.regime == "mask_handling"
    idx_nan(end,end) = true;
    LULC_Properties.idx_imp(zone_id == 2) = true;
end

depths.d_t(idx_nan) = nan;

input_evaporation = Cfg.evap_map_mm_d * ones(size(elevation));
input_transpiration = Cfg.transp_map_mm_d * ones(size(elevation));
input_evaporation(idx_nan) = nan;
input_transpiration(idx_nan) = nan;

Hydro_States = struct();
Hydro_States.ETP = Cfg.internal_etp_mm_d * ones(size(elevation));
Hydro_States.Ep = Cfg.open_water_ep_mm_d * ones(size(elevation));
Hydro_States.ETP(idx_nan) = nan;
Hydro_States.Ep(idx_nan) = nan;

BC_States = struct();
errors = zeros(3,1);

Soil_Properties = make_layered_soil(Cfg, idx_nan, zone_id);
min_soil_moisture = Cfg.min_soil_storage_mm * ones(size(elevation));

if Cfg.regime == "mask_handling"
    Soil_Properties.Layers.near_surface_storage_mm(zone_id == 1) = 1;
    Soil_Properties.Layers.root_zone_storage_mm(zone_id == 1) = 0;
    min_soil_moisture(zone_id == 1) = 2;
    Soil_Properties.Layers.near_surface_storage_mm(zone_id == 2) = 10;
    Soil_Properties.Layers.root_zone_storage_mm(zone_id == 2) = 10;
    Soil_Properties.Layers.near_surface_storage_mm(zone_id == 3) = 10;
    Soil_Properties.Layers.root_zone_storage_mm(zone_id == 3) = 10;
    min_soil_moisture(zone_id ~= 1) = 0;
    Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end

min_soil_moisture(idx_nan) = nan;

Initial = capture_state(Soil_Properties, depths);
Reference = reference_extraction(Cfg, Initial, depths.d_t, Hydro_States, ...
    input_evaporation, input_transpiration, LULC_Properties.idx_imp, ...
    min_soil_moisture, idx_nan, dt_days);

run(fullfile(functions_dir, 'Evaporation_Evapotranspiration_Module.m'));

Model = capture_state(Soil_Properties, depths);
soil_et_model_mm = Hydro_States.ETR .* dt_days;
surface_evap_model_mm = BC_States.delta_E;

soil_et_error = soil_et_model_mm - Reference.soil_et_mm;
surface_evap_error = surface_evap_model_mm - Reference.surface_evap_mm;
soil_storage_error = Model.soil_storage_mm - Reference.soil_storage_final_mm;
surface_storage_error = Model.surface_storage_mm - Reference.surface_storage_final_mm;

initial_volume_m3 = storage_volume_m3(Initial.soil_storage_mm + Initial.surface_storage_mm, C_a, idx_nan);
final_volume_m3 = storage_volume_m3(Model.soil_storage_mm + Model.surface_storage_mm, C_a, idx_nan);
removed_volume_m3 = storage_volume_m3(soil_et_model_mm + surface_evap_model_mm, C_a, idx_nan);
mass_residual_m3 = initial_volume_m3 - final_volume_m3 - removed_volume_m3;

module_error_m3 = errors(3);
max_soil_et_error = max_abs_omitnan(soil_et_error);
max_surface_evap_error = max_abs_omitnan(surface_evap_error);
max_soil_storage_error = max_abs_omitnan(soil_storage_error);
max_surface_storage_error = max_abs_omitnan(surface_storage_error);
max_storage_error = max(max_soil_storage_error, max_surface_storage_error);
nan_mismatch = count_nan_mismatch(soil_et_model_mm, Reference.soil_et_mm) + ...
    count_nan_mismatch(surface_evap_model_mm, Reference.surface_evap_mm);

passed = max_soil_et_error < 1e-10 && ...
    max_surface_evap_error < 1e-10 && ...
    max_storage_error < 1e-10 && ...
    abs(mass_residual_m3) < 1e-6 && ...
    abs(module_error_m3) < 1e-6 && ...
    nan_mismatch == 0;

CellRows = make_cell_rows(Cfg, zone_id, zone_name, idx_nan, ...
    nan(size(elevation)), nan(size(elevation)), nan(size(elevation)), nan(size(elevation)), ...
    Reference.soil_et_mm, soil_et_model_mm, ...
    Reference.surface_evap_mm, surface_evap_model_mm, ...
    Initial.soil_storage_mm, Model.soil_storage_mm, ...
    Initial.surface_storage_mm, Model.surface_storage_mm);

DiagRow = table( ...
    Cfg.case_id, Cfg.case_name, Cfg.regime, ...
    nan, nan, ...
    sum_omitnan(Reference.soil_et_mm), sum_omitnan(soil_et_model_mm), ...
    sum_omitnan(Reference.surface_evap_mm), sum_omitnan(surface_evap_model_mm), ...
    mean_omitnan(Initial.soil_storage_mm), mean_omitnan(Model.soil_storage_mm), ...
    mean_omitnan(Initial.surface_storage_mm), mean_omitnan(Model.surface_storage_mm), ...
    nan, nan, max_soil_et_error, max_surface_evap_error, ...
    max_storage_error, mass_residual_m3, module_error_m3, nan_mismatch, passed, ...
    'VariableNames', diagnostic_names());
end

function names = diagnostic_names()
names = {'case_id','case_name','regime', ...
    'mean_reference_etp_mm_d','mean_model_etp_mm_d', ...
    'soil_et_expected_sum_mm','soil_et_model_sum_mm', ...
    'surface_evap_expected_sum_mm','surface_evap_model_sum_mm', ...
    'initial_soil_storage_mean_mm','final_soil_storage_mean_mm', ...
    'initial_surface_storage_mean_mm','final_surface_storage_mean_mm', ...
    'max_etp_error_mm_d','max_ep_error_mm_d', ...
    'max_soil_et_error_mm','max_surface_evap_error_mm', ...
    'max_storage_error_mm','mass_residual_m3','module_error_m3', ...
    'nan_mismatch_count','passed'};
end

function [DEM, zone_id, zone_name, idx_nan, albedo] = make_vtilted_domain()
DEM = [112 100 112; 111 100 111; 110 100 110];
zone_id = [1 2 3; 1 2 3; 1 2 3];
zone_name = strings(size(zone_id));
zone_name(zone_id == 1) = "left_hillslope";
zone_name(zone_id == 2) = "channel_strip";
zone_name(zone_id == 3) = "right_hillslope";
idx_nan = false(size(DEM));
albedo = 0.23 * ones(size(DEM));
albedo(zone_id == 2) = 0.08;
end

function Soil_Properties = make_layered_soil(Cfg, idx_nan, zone_id)
template = double(zone_id);
Soil_Properties = struct();
layers = struct();
layers.root_depth_m = Cfg.root_depth_m * ones(size(template));
layers.near_surface_thickness_m = Cfg.near_thickness_m * ones(size(template));
layers.root_zone_thickness_m = Cfg.root_thickness_m * ones(size(template));
layers.transmission_thickness_m = Cfg.trans_thickness_m * ones(size(template));
layers.near_surface_capacity_mm = 80 * ones(size(template));
layers.root_zone_capacity_mm = 200 * ones(size(template));
layers.transmission_capacity_mm = 400 * ones(size(template));
layers.near_surface_storage_mm = Cfg.near_storage_mm * ones(size(template));
layers.root_zone_storage_mm = Cfg.root_storage_mm * ones(size(template));
layers.transmission_storage_mm = Cfg.trans_storage_mm * ones(size(template));
layers.near_surface_capacity_mm(idx_nan) = nan;
layers.root_zone_capacity_mm(idx_nan) = nan;
layers.transmission_capacity_mm(idx_nan) = nan;
layers.near_surface_storage_mm(idx_nan) = nan;
layers.root_zone_storage_mm(idx_nan) = nan;
layers.transmission_storage_mm(idx_nan) = nan;
Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end

function State = capture_state(Soil_Properties, depths)
State = struct();
State.near_mm = Soil_Properties.Layers.near_surface_storage_mm;
State.root_mm = Soil_Properties.Layers.root_zone_storage_mm;
State.trans_mm = Soil_Properties.Layers.transmission_storage_mm;
State.soil_storage_mm = Soil_Properties.I_t;
State.surface_storage_mm = depths.d_t;
end

function Reference = reference_extraction(Cfg, Initial, d_t, Hydro_States, ...
    input_evaporation, input_transpiration, idx_imp, min_soil_moisture, idx_nan, dt_days)

surface_evap_mm = zeros(size(d_t));
near_final = Initial.near_mm;
root_final = Initial.root_mm;
trans_final = Initial.trans_mm;
surface_final = d_t;

idx_open = d_t > 0 & ~idx_nan;
surface_evap_mm(idx_open) = min(d_t(idx_open), Hydro_States.Ep(idx_open) .* dt_days);
surface_final(idx_open) = d_t(idx_open) - surface_evap_mm(idx_open);

if Cfg.mode == "map"
    demand_mm = (input_evaporation + input_transpiration) .* dt_days;
else
    demand_mm = Hydro_States.ETP .* dt_days;
end

soil_storage = Initial.near_mm + Initial.root_mm + Initial.trans_mm;
idx_soil = d_t == 0 & soil_storage > min_soil_moisture & ~idx_imp & ~idx_nan;

near_access_fraction = ones(size(d_t));
idx_partial_near = Cfg.near_thickness_m > 0 & ...
    Cfg.root_depth_m > 0 & Cfg.root_depth_m < Cfg.near_thickness_m;
if idx_partial_near
    near_access_fraction(:) = Cfg.root_depth_m / Cfg.near_thickness_m;
end
if Cfg.root_depth_m <= 0
    near_access_fraction(:) = 0;
end

remaining = max(demand_mm, 0);
remaining(~idx_soil) = 0;

near_available = max(near_final .* near_access_fraction, 0);
take_near = min(remaining, near_available);
near_final = near_final - take_near;
remaining = remaining - take_near;

root_available = max(root_final, 0);
take_root = min(remaining, root_available);
root_final = root_final - take_root;
remaining = remaining - take_root;

soil_et_mm = max(demand_mm, 0) - max(remaining, 0);
soil_et_mm(~idx_soil) = 0;

soil_et_mm(idx_nan) = nan;
surface_evap_mm(idx_nan) = nan;
near_final(idx_nan) = nan;
root_final(idx_nan) = nan;
trans_final(idx_nan) = nan;
surface_final(idx_nan) = nan;

Reference = struct();
Reference.soil_et_mm = soil_et_mm;
Reference.surface_evap_mm = surface_evap_mm;
Reference.near_final_mm = near_final;
Reference.root_final_mm = root_final;
Reference.trans_final_mm = trans_final;
Reference.soil_storage_final_mm = near_final + root_final + trans_final;
Reference.surface_storage_final_mm = surface_final;
end

function [ETP, Ep] = reference_evapotranspiration(DEM, Temp, Temp_max, Temp_min, ...
    Day_year, lat, U_2_input, Krs, albedo, G_input)

theta_S_B = 4.903e-9;
idx = isinf(DEM) | isnan(DEM) | DEM < -200 | DEM > 99999;

T = Temp;
Tmax = Temp_max;
Tmin = Temp_min;
U_2 = U_2_input;
G = G_input;
T(idx) = nan;
Tmax(idx) = nan;
Tmin(idx) = nan;
U_2(idx) = nan;
G(idx) = nan;
albedo(idx) = nan;

J = Day_year * ones(size(DEM));
J(idx) = nan;

phi = lat * pi / 180;
dec_sol = 0.409 * sin((2 * pi / 365) .* J - 1.39);
X = 1 - tan(phi).^2 .* tan(dec_sol).^2;
X(X <= 0) = 0.00001;
ws = (pi / 2) - atan((-tan(phi) .* tan(dec_sol)) ./ sqrt(X));
dr = 1 + 0.033 * cos((2 * pi / 365) .* J);
Ra = (118.08 / pi) .* dr .* ...
    (ws .* sin(phi) .* sin(dec_sol) + cos(phi) .* cos(dec_sol) .* sin(ws));
Rs = Krs .* Ra .* sqrt(Tmax - Tmin);
Rso = (0.75 + 2e-5 .* DEM) .* Ra;
Rns = (1 - albedo) .* Rs;
e_s = 0.6108 * exp((17.27 .* T) ./ (T + 237.3));
e_a = 0.61 * exp((17.27 .* Tmin) ./ (Tmin + 237.3));
Rnl = theta_S_B .* (((Tmax + 273.16).^4 + (Tmin + 273.16).^4) ./ 2) .* ...
    (0.34 - 0.14 .* sqrt(e_a)) .* (1.35 .* (Rs ./ Rso) - 0.35);
Rn = Rns - Rnl;
delta = (4098 .* (0.6108 .* exp((17.27 .* T) ./ (T + 237.3)))) ./ ...
    ((T + 237.3).^2);
Patm = 101.3 .* (((293 - 0.0065 .* DEM) ./ 293).^5.26);
gamma = 0.665e-3 .* Patm;

ETP = (0.408 .* delta .* (Rn - G) + ...
    (gamma .* 900 .* U_2 .* (e_s - e_a)) ./ (T + 273)) ./ ...
    (delta + gamma .* (1 + 0.34 .* U_2));
Ep = (0.408 .* delta .* (Rn - G) + ...
    (gamma .* 900 .* U_2 .* (e_s - e_a)) ./ (T + 273)) ./ ...
    (delta + gamma);

ETP(idx) = nan;
Ep(idx) = nan;
end

function CellRows = make_cell_rows(Cfg, zone_id, zone_name, idx_nan, ...
    etp_ref, etp_model, ep_ref, ep_model, soil_et_ref, soil_et_model, ...
    surface_evap_ref, surface_evap_model, soil_initial, soil_final, ...
    surface_initial, surface_final)

[row_idx, col_idx] = ndgrid(1:size(zone_id,1), 1:size(zone_id,2));
CellRows = table( ...
    repmat(Cfg.case_id, numel(zone_id), 1), ...
    repmat(Cfg.case_name, numel(zone_id), 1), ...
    repmat(Cfg.regime, numel(zone_id), 1), ...
    row_idx(:), col_idx(:), zone_id(:), zone_name(:), idx_nan(:), ...
    etp_ref(:), etp_model(:), ep_ref(:), ep_model(:), ...
    soil_et_ref(:), soil_et_model(:), ...
    surface_evap_ref(:), surface_evap_model(:), ...
    soil_initial(:), soil_final(:), surface_initial(:), surface_final(:), ...
    'VariableNames', {'case_id','case_name','regime','row','col','zone_id', ...
    'zone_name','is_invalid','etp_reference_mm_d','etp_model_mm_d', ...
    'ep_reference_mm_d','ep_model_mm_d','soil_et_reference_mm', ...
    'soil_et_model_mm','surface_evap_reference_mm', ...
    'surface_evap_model_mm','soil_storage_initial_mm', ...
    'soil_storage_final_mm','surface_storage_initial_mm', ...
    'surface_storage_final_mm'});
end

function make_et_figures(Diagnostics, cell_dir, fig_dir)

close all force

f = figure('Visible','off', 'Color','w', 'Units','pixels', ...
    'Position', [100 100 1650 950]);
tiledlayout(f, 2, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile
idx_formula = Diagnostics.regime == "potential_et_formula";
bar([Diagnostics.mean_reference_etp_mm_d(idx_formula), ...
    Diagnostics.mean_model_etp_mm_d(idx_formula)])
set(gca, 'XTickLabel', {'Reference','HydroPol2D'})
ylabel('Mean ETP (mm/day)')
title('Potential ET formula')
grid on

nexttile
idx = Diagnostics.regime ~= "potential_et_formula";
bar(categorical(Diagnostics.case_id(idx)), ...
    [Diagnostics.soil_et_expected_sum_mm(idx), ...
    Diagnostics.soil_et_model_sum_mm(idx), ...
    Diagnostics.surface_evap_expected_sum_mm(idx), ...
    Diagnostics.surface_evap_model_sum_mm(idx)])
ylabel('Domain sum (mm over cells)')
title('Actual ET and open-water evaporation')
legend({'Soil ref','Soil model','Surface ref','Surface model'}, 'Location','best')
set(gca, 'XTickLabelRotation', 25)
grid on

nexttile
bar(categorical(Diagnostics.case_id(idx)), ...
    [Diagnostics.initial_soil_storage_mean_mm(idx), ...
    Diagnostics.final_soil_storage_mean_mm(idx), ...
    Diagnostics.initial_surface_storage_mean_mm(idx), ...
    Diagnostics.final_surface_storage_mean_mm(idx)])
ylabel('Mean storage (mm)')
title('Storage before and after ET')
legend({'Soil initial','Soil final','Surface initial','Surface final'}, 'Location','best')
set(gca, 'XTickLabelRotation', 25)
grid on

nexttile
error_floor = 1e-16;
error_values = [Diagnostics.max_etp_error_mm_d, Diagnostics.max_soil_et_error_mm, ...
    Diagnostics.max_surface_evap_error_mm, Diagnostics.max_storage_error_mm];
error_values(~isfinite(error_values)) = error_floor;
error_values = max(error_values, error_floor);
bar(categorical(Diagnostics.case_id), error_values)
set(gca, 'YScale', 'log')
ylabel('Absolute error')
title('Validation errors')
legend({'ETP mm/day','Soil ET mm','Surface E mm','Storage mm'}, 'Location','best')
ylim([error_floor 1e-8])
set(gca, 'XTickLabelRotation', 25)
grid on

exportgraphics(f, fullfile(fig_dir, 'VTilted_ET_Overview.png'), 'Resolution', 220);
exportgraphics(f, fullfile(fig_dir, 'VTilted_ET_Overview.pdf'), 'ContentType','vector');
close(f)

files = dir(fullfile(cell_dir, 'P1-ET-VT-*_cells.csv'));
for i = 1:numel(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'TextType', 'string');
    make_case_map_figure(T, fig_dir);
end
end

function make_case_map_figure(T, fig_dir)
case_id = T.case_id(1);
regime = T.regime(1);

soil_et = reshape(T.soil_et_model_mm, 3, 3);
surface_e = reshape(T.surface_evap_model_mm, 3, 3);
soil_final = reshape(T.soil_storage_final_mm, 3, 3);
etp = reshape(T.etp_model_mm_d, 3, 3);

f = figure('Visible','off', 'Color','w', 'Units','pixels', ...
    'Position', [100 100 1000 760]);
tiledlayout(f, 2, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile
imagesc(etp); colorbar; axis image
title('ETP (mm/day)')

nexttile
imagesc(soil_et); colorbar; axis image
title('Soil ET (mm)')

nexttile
imagesc(surface_e); colorbar; axis image
title('Surface evaporation (mm)')

nexttile
imagesc(soil_final); colorbar; axis image
title('Final soil storage (mm)')

sgtitle(sprintf('%s %s', case_id, regime), 'Interpreter','none')
safe_id = matlab.lang.makeValidName(char(case_id));
exportgraphics(f, fullfile(fig_dir, [safe_id '_et_maps.png']), 'Resolution', 220);
exportgraphics(f, fullfile(fig_dir, [safe_id '_et_maps.pdf']), 'ContentType','vector');
close(f)
end

function value = storage_volume_m3(depth_mm, C_a, idx_nan)
depth_mm(idx_nan) = nan;
value = sum_omitnan(C_a .* depth_mm ./ 1000);
end

function value = max_abs_omitnan(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    value = nan;
else
    value = max(abs(x));
end
end

function value = mean_omitnan(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    value = nan;
else
    value = mean(x);
end
end

function value = sum_omitnan(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    value = 0;
else
    value = sum(x);
end
end

function count = count_nan_mismatch(a, b)
count = nnz(isnan(a) ~= isnan(b));
end
