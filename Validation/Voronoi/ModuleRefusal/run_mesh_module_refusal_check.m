% VAL-MESH-REFUSE-001: a mesh run stops when a module it would skip is on.
%
% A MATLAB mesh run (flag_voronoi = 1, or a Case Builder quadtree package)
% reads its physics from the prepared case, not from the grid inputs. Before
% this check, snow, reservoirs, water quality, human risk and the deep
% aquifer switched on in the grid inputs were silently left out of the mesh
% run. HydroPol2D-Python stops such a run (MeshFeatureNotSupportedError).
% The check covers:
%   A. hp2d_refuse_mesh_modules refuses each module switch on its own, and
%      lets a run with every switch off through.
%   B. A module the prepared case carries (config.reservoirs) is not refused;
%      an empty one is.
%   C. HydroPol2D_Voronoi_Selection passes the bypass module switches on.
%   D. HydroPol2D_Run_Prepared_Quadtree refuses a Builder deep-aquifer switch
%      before it reads any raster.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

refusal_id = 'HydroPol2D:MeshModuleNotSupported';
module_flags = {'flag_snow_modeling','flag_reservoir','flag_waterquality', ...
    'flag_human_instability','flag_deep_aquifer','flag_dam_break', ...
    'flag_detention_ponds','flag_urban_drainage','flag_green_infrastructure'};

% A. Each switch on its own is refused; all off runs.
all_off = cell2struct(num2cell(zeros(1, numel(module_flags))), module_flags, 2);
refused_alone = false(numel(module_flags), 1);
for k = 1:numel(module_flags)
    one_on = all_off;
    one_on.(module_flags{k}) = 1;
    refused_alone(k) = strcmp(error_id(@() hp2d_refuse_mesh_modules(one_on, struct(), 'check')), refusal_id);
end
all_off_runs = isempty(error_id(@() hp2d_refuse_mesh_modules(all_off, struct(), 'check')));
human_risk_mode3 = all_off; human_risk_mode3.flag_human_instability = 3;
mode3_refused = strcmp(error_id(@() hp2d_refuse_mesh_modules(human_risk_mode3, struct(), 'check')), refusal_id);

% B. A module the prepared case carries is not refused.
reservoir_on = all_off; reservoir_on.flag_reservoir = 1;
carried_runs = isempty(error_id(@() hp2d_refuse_mesh_modules(reservoir_on, ...
    struct('reservoirs', struct('source_cell', 1)), 'check')));
empty_refused = strcmp(error_id(@() hp2d_refuse_mesh_modules(reservoir_on, ...
    struct('reservoirs', []), 'check')), refusal_id);

% C. The selection passes the bypass module switches on.
bypass_dir = tempname; mkdir(bypass_dir);
bypass_script = fullfile(bypass_dir, 'mesh_refusal_bypass.m');
fid = fopen(bypass_script, 'w');
fprintf(fid, ['InputData_Bypass = struct();\n' ...
    'InputData_Bypass.flags.flag_voronoi = 1;\n' ...
    'InputData_Bypass.flags.flag_snow_modeling = 1;\n' ...
    'InputData_Bypass.flags.flag_reservoir = 0;\n' ...
    'InputData_Bypass.Voronoi.case_file = ''prepared_case.mat'';\n']);
fclose(fid);
selection = HydroPol2D_Voronoi_Selection('bypass', '', bypass_script, struct());
selection_passes = isfield(selection, 'module_flags') && ...
    isfield(selection.module_flags, 'flag_snow_modeling') && ...
    selection.module_flags.flag_snow_modeling == 1 && ...
    isfield(selection.module_flags, 'flag_reservoir') && ...
    selection.module_flags.flag_reservoir == 0;
selection_refuses = selection_passes && strcmp(error_id(@() hp2d_refuse_mesh_modules( ...
    selection.module_flags, struct(), 'check')), refusal_id);

% D. The quadtree runner refuses a Builder deep-aquifer switch first.
quadtree_refuses = quadtree_error(true);
quadtree_off_id = quadtree_error(false);
quadtree_off_passes_refusal = ~strcmp(quadtree_off_id, refusal_id);
quadtree_refuses = strcmp(quadtree_refuses, refusal_id);

Checks = table(string(module_flags'), refused_alone, 'VariableNames', {'flag','refused_alone'});
passed = all(refused_alone) && all_off_runs && mode3_refused && carried_runs && ...
    empty_refused && selection_passes && selection_refuses && quadtree_refuses && ...
    quadtree_off_passes_refusal;
status = "fail";
if passed; status = "pass"; end

PassFail = table("VAL-MESH-REFUSE-001", "Mesh runs refuse modules they would skip", ...
    sum(refused_alone), numel(module_flags), all_off_runs, carried_runs, empty_refused, ...
    selection_refuses, quadtree_refuses, status, true, ...
    'VariableNames', {'case_id','case_name','refused_switches','module_switches', ...
    'all_off_runs','carried_module_runs','empty_module_refused', ...
    'selection_refuses','quadtree_refuses','status','report_ready'});
writetable(Checks, fullfile(out_dir, 'Mesh_Module_Refusals.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(Checks);
disp(PassFail);
rmdir(bypass_dir, 's');
if ~passed
    error('VAL-MESH-REFUSE-001 failed.');
end
fprintf('VAL-MESH-REFUSE-001 passed.\n');

function id = error_id(action)
id = '';
try
    action();
catch caught
    id = caught.identifier;
end
end

function id = quadtree_error(deep_aquifer)
case_root = tempname; mkdir(case_root);
config = struct('flags', struct('deep_aquifer', deep_aquifer, 'groundwater', true));
fid = fopen(fullfile(case_root, 'hydropol_config.json'), 'w');
fprintf(fid, '%s', jsonencode(config));
fclose(fid);
id = error_id(@() HydroPol2D_Run_Prepared_Quadtree(case_root));
if isempty(id), id = 'ran'; end
rmdir(case_root, 's');
end
