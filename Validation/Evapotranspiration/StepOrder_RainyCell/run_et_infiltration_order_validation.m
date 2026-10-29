% VAL-ET-ORDER-001: ET is taken after infiltration in the legacy step.
%
% One Hydrological_Model step on a 3 x 3 grid with 2 mm of rain, ET on and
% the legacy soil bucket. The order is snow, canopy, infiltration, open-water
% evaporation, soil ET, as in HydroPol2D-Python (core/engine.py). Expected:
%   - a cell whose soil takes all the rain ends the step dry and loses soil
%     ET (ETP x Kc), with no open-water evaporation;
%   - a cell whose soil is full keeps the rain ponded and loses open-water
%     evaporation (Ep), with no soil ET;
%   - a cell with no rain loses soil ET;
%   - the ledger rain = d(soil) + d(surface) + evaporation + soil ET closes.
% Before the fix, ET ran first: the rained-on cell counted as open water, so
% it lost Ep from the rain and no soil ET.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

addpath(functions_dir, '-begin');

[Cells, Ledger] = run_one_step(functions_dir);

tol = 1e-12;
infiltrating = Cells.zone == "infiltrating";
ponded = Cells.zone == "ponded";
no_rain = Cells.zone == "no_rain";

checks = table();
checks = [checks; check_row("Infiltrating cell ends dry", ...
    max(Cells.surface_after_mm(infiltrating)), tol)];
checks = [checks; check_row("Infiltrating cell: no open-water evaporation", ...
    max(Cells.open_water_evaporation_mm(infiltrating)), tol)];
checks = [checks; check_row("Infiltrating cell: soil ET = ETP x Kc", ...
    max(abs(Cells.soil_et_mm_day(infiltrating) - Cells.expected_soil_et_mm_day(infiltrating))), tol)];
checks = [checks; check_row("Ponded cell: open-water evaporation = Ep", ...
    max(abs(Cells.open_water_evaporation_mm(ponded) - Cells.expected_open_water_mm(ponded))), tol)];
checks = [checks; check_row("Ponded cell: no soil ET", ...
    max(abs(Cells.soil_et_mm_day(ponded))), tol)];
checks = [checks; check_row("No-rain cell: soil ET = ETP x Kc", ...
    max(abs(Cells.soil_et_mm_day(no_rain) - Cells.expected_soil_et_mm_day(no_rain))), tol)];
checks = [checks; check_row("Ledger residual (m3)", abs(Ledger.residual_m3), 1e-10)];

passed = all(checks.passed);
Diagnostics = table("VAL-ET-ORDER-001", "ET after infiltration in the legacy step", ...
    max(checks.value(1:end-1)), Ledger.residual_m3, passed, ...
    'VariableNames', {'case_id','case_name','max_abs_error','mass_residual_m3','passed'});
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-ET-ORDER-001", "ET after infiltration in the legacy step", ...
    status, true, 'VariableNames', {'case_id','case_name','status','report_ready'});

writetable(Cells, fullfile(out_dir, 'ET_Order_Cells.csv'));
writetable(checks, fullfile(out_dir, 'ET_Order_Checks.csv'));
writetable(Ledger, fullfile(out_dir, 'ET_Order_Ledger.csv'));
writetable(Diagnostics, fullfile(out_dir, 'ET_Order_Diagnostics.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));

disp(Cells);
disp(checks);
disp(Diagnostics);
if ~passed
    error('VAL-ET-ORDER-001 failed.');
end
fprintf('VAL-ET-ORDER-001 passed.\n');

function [Cells, Ledger] = run_one_step(functions_dir)
% Row 1: soil takes all the rain. Row 2: soil already full, the rain ponds.
% Row 3: no rain.
ny = 3; nx = 3;
zone = repmat(["infiltrating"; "ponded"; "no_rain"], 1, nx);
rain_mm = 2 * double(zone ~= "no_rain");

resolution = 10;
time_step = 15; %#ok<NASGU> % min
k = 1; %#ok<NASGU>
dt_days = 15 / 60 / 24;
catch_index = 0; %#ok<NASGU>
errors = zeros(6,1); %#ok<NASGU>

flags = struct('flag_subgrid', 0, 'flag_overbanks', 0, 'flag_snow_modeling', 0, ...
    'flag_ETP', 1, 'flag_abstraction', 0, 'flag_infiltration', 1, ...
    'flag_groundwater_modeling', 0, 'flag_layered_subsurface', 0, ...
    'flag_input_ETP_map', 0); %#ok<NASGU>

elevation = 100 * ones(ny, nx);
idx_nan = false(ny, nx);
zero_matrix = zeros(ny, nx); %#ok<NASGU>
DEM_raster = struct('Z', elevation, 'cellsize', resolution); %#ok<NASGU>
Elevation_Properties = struct('elevation_cell', elevation); %#ok<NASGU>
Wshed_Properties = struct('Resolution', resolution, 'cell_area', resolution^2); %#ok<NASGU>
LAI_raster = struct('Z', zeros(ny, nx)); %#ok<NASGU>
C_a = resolution^2 * ones(ny, nx);

soil_depth_m = 1.0;
theta_sat = 0.45;
theta_r = 0.05;
capacity_mm = 1000 * soil_depth_m * (theta_sat - theta_r);
I_0 = 100 * ones(ny, nx);
I_0(zone == "ponded") = capacity_mm;
Soil_Properties = struct('I_t', I_0, 'Soil_Depth', soil_depth_m * ones(ny, nx), ...
    'theta_sat', theta_sat * ones(ny, nx), 'theta_r', theta_r * ones(ny, nx), ...
    'theta_i', 0.2 * ones(ny, nx), 'ksat', 50 * ones(ny, nx), ...
    'alpha_vg', 3 * ones(ny, nx), 'n_vg', 1.5 * ones(ny, nx));
min_soil_moisture = 10 * ones(ny, nx); %#ok<NASGU>

etp_mm_day = 4;
ep_mm_day = 5;
kc = 0.8;
LULC_Properties = struct('frac_perv', ones(ny, nx), 'idx_imp', false(ny, nx), ...
    'Kc', kc * ones(ny, nx)); %#ok<NASGU>
Hydro_States = struct('ETP', etp_mm_day * ones(ny, nx), 'Ep', ep_mm_day * ones(ny, nx), ...
    'S', zeros(ny, nx)); %#ok<NASGU>

depths = struct('d_t', zeros(ny, nx), 'd_p', zeros(ny, nx));
% Water table at the soil base: the whole soil column is unsaturated.
BC_States = struct('delta_p_agg', rain_mm, 'inflow', zeros(ny, nx), ...
    'h_t', elevation - soil_depth_m);

d_0 = depths.d_t;
run(fullfile(functions_dir, 'Hydrological_Model.m'));

open_water_mm = BC_States.delta_E;
soil_et_mm_day = Hydro_States.ETR;
expected_soil_et = etp_mm_day * kc * double(zone ~= "ponded");
expected_open_water = ep_mm_day * dt_days * double(zone == "ponded");

Cells = table(zone(:), rain_mm(:), I_0(:), Soil_Properties.I_t(:), d_0(:), depths.d_t(:), ...
    open_water_mm(:), soil_et_mm_day(:), expected_soil_et(:), expected_open_water(:), ...
    'VariableNames', {'zone','rain_mm','soil_before_mm','soil_after_mm', ...
    'surface_before_mm','surface_after_mm','open_water_evaporation_mm', ...
    'soil_et_mm_day','expected_soil_et_mm_day','expected_open_water_mm'});

area = resolution^2;
rain_m3 = sum(rain_mm(:)) / 1000 * area;
d_soil_m3 = sum(Soil_Properties.I_t(:) - I_0(:)) / 1000 * area;
d_surface_m3 = sum(depths.d_t(:) - d_0(:)) / 1000 * area;
evaporation_m3 = sum(open_water_mm(:)) / 1000 * area;
soil_et_m3 = sum(soil_et_mm_day(:)) * dt_days / 1000 * area;
residual_m3 = rain_m3 - d_soil_m3 - d_surface_m3 - evaporation_m3 - soil_et_m3;
Ledger = table(rain_m3, d_soil_m3, d_surface_m3, evaporation_m3, soil_et_m3, residual_m3);
end

function Row = check_row(name, value, tolerance)
Row = table(string(name), value, tolerance, value <= tolerance, ...
    'VariableNames', {'check','value','tolerance','passed'});
end
