function report = run_deep_aquifer_validation()
%RUN_DEEP_AQUIFER_VALIDATION Focused checks for distributed deep baseflow.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'HydroPol2D_Functions'));

subbasin = [1 1; 2 2];
weight = [1 3; 2 1];
zero_recharge = zeros(2);
cell_area_m2 = 25;
fraction = 0.2;
tau_days = 10/86400;

[initial,~,~] = hp2d_deep_aquifer_step(struct(),zero_recharge,subbasin, ...
    weight,fraction,tau_days,cell_area_m2,1);
initial.deep_aquifer_storage_m3 = [100;50];
initial.cumulative_deep_aquifer_recharge_m3(:) = 0;
initial.cumulative_deep_aquifer_discharge_m3(:) = 0;

[one_step,~,surface_one_mm] = hp2d_deep_aquifer_step(initial,zero_recharge, ...
    subbasin,weight,fraction,tau_days,cell_area_m2,10);
many_steps = initial;
surface_many_m3 = 0;
for k = 1:10
    [many_steps,~,surface_mm] = hp2d_deep_aquifer_step(many_steps,zero_recharge, ...
        subbasin,weight,fraction,tau_days,cell_area_m2,1);
    surface_many_m3 = surface_many_m3 + sum(surface_mm,'all')/1000*cell_area_m2;
end
one_discharge_m3 = sum(surface_one_mm,'all')/1000*cell_area_m2;
assert(max(abs(one_step.deep_aquifer_storage_m3- ...
    many_steps.deep_aquifer_storage_m3)) < 1e-10);
assert(abs(one_discharge_m3-surface_many_m3) < 1e-10);

recharge_rate_m_s = [0.01 0.02;0.03 0.04];
[forced,shallow_rate,surface_mm] = hp2d_deep_aquifer_step(initial, ...
    recharge_rate_m_s,subbasin,weight,fraction,tau_days,cell_area_m2,10);
input_m3 = sum(recharge_rate_m_s,'all')*10*cell_area_m2*fraction;
storage_change_m3 = sum(forced.deep_aquifer_storage_m3-initial.deep_aquifer_storage_m3);
discharge_m3 = sum(surface_mm,'all')/1000*cell_area_m2;
assert(abs(input_m3-storage_change_m3-discharge_m3) < 1e-10);
assert(max(abs(shallow_rate-(1-fraction)*recharge_rate_m_s),[],'all') < 1e-12);

fixture_path = fullfile(fileparts(mfilename('fullpath')), 'deep_aquifer_step.json');
fixture = jsondecode(fileread(fixture_path));
fixture_state = struct( ...
    'deep_aquifer_storage_m3', double(fixture.storage_m3(:)), ...
    'cumulative_deep_aquifer_recharge_m3', zeros(2,1), ...
    'cumulative_deep_aquifer_discharge_m3', zeros(2,1));
[fixture_state,fixture_shallow,fixture_surface_mm] = hp2d_deep_aquifer_step( ...
    fixture_state,double(fixture.recharge_rate_m_s),double(fixture.subbasin_id), ...
    double(fixture.channel_weight),double(fixture.recharge_fraction), ...
    double(fixture.response_time_days),double(fixture.cell_area_m2),double(fixture.dt_s));
fixture_surface_m3 = fixture_surface_mm ./ 1000 .* double(fixture.cell_area_m2);
fixture_expected = fixture.expected;
assert(local_relative_error(fixture_state.deep_aquifer_storage_m3, ...
    double(fixture_expected.storage_m3(:))) < 1e-8);
assert(local_relative_error(fixture_state.deep_aquifer_discharge_m3_s .* fixture.dt_s, ...
    double(fixture_expected.discharge_m3(:))) < 1e-8);
assert(local_relative_error(fixture_surface_m3, ...
    double(fixture_expected.surface_addition_m3)) < 1e-8);
assert(local_relative_error(fixture_shallow, ...
    double(fixture_expected.shallow_recharge_rate_m_s)) < 1e-8);

report = struct( ...
    'timestep_invariant',true, ...
    'mass_conservative',true, ...
    'distributed_return_conservative',true, ...
    'recharge_partition_correct',true, ...
    'matlab_python_fixture_relative_error_below_1e_8',true);
fprintf('Deep-aquifer validation passed: timestep invariance, partition, and mass balance.\n');
end

function error_value = local_relative_error(actual, expected)
actual = double(actual);
expected = double(expected);
scale = max(abs(expected), 1e-15);
error_value = max(abs(actual-expected) ./ scale, [], 'all');
end
