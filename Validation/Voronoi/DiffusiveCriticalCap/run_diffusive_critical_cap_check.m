% VAL-MESH-DIFF-CRIT-001: mesh diffusive flux capped at critical flow.
%
% HydroPol2D-Python caps the mesh diffusive-wave face discharge at the
% critical discharge h sqrt(g h) when critical_flow is on (voronoi.py).
% MATLAB had no cap (audit L24), so on steep, smooth ground the diffusive
% wave moved more water than a critical section can pass.
% Two 10 m cells, 0.2 m of water on a 1 m drop, Manning n 0.01: Manning's
% discharge per width is far above critical. With critical_flow the face
% discharge must equal 0.2 sqrt(9.81 x 0.2); without it, Manning's value.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

mesh = struct('n_cells', 2, 'n_edges', 1, 'edge_owner', 1, 'edge_neighbor', 2, ...
    'edge_length', 10, 'edge_distance', 10, 'surface_area', [100; 100], 'surface_bed', [1; 0]);
volume = [0.2 * 100; 0];
n = 0.01; g = 9.81;
capped = Voronoi_Diffusive_Face_Flux(mesh, volume, n, critical_flow=true, gravity=g);
uncapped = Voronoi_Diffusive_Face_Flux(mesh, volume, n);
h = 0.2;                                     % face depth: owner stage 1.2 - higher bed 1.0
slope = (1.2 - 0) / 10;
manning = h^(5/3) / n * sqrt(slope);
critical = h * sqrt(g * h);
cap_error = abs(capped - critical);
manning_error = abs(uncapped - manning);
passed = manning > critical && cap_error < 1e-12 && manning_error < 1e-12;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-DIFF-CRIT-001", "Mesh diffusive flux capped at critical flow", ...
    manning, critical, capped, cap_error, status, true, 'VariableNames', {'case_id', ...
    'case_name','manning_q_m2_s','critical_q_m2_s','capped_q_m2_s','cap_error_m2_s', ...
    'status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-DIFF-CRIT-001 failed.');
end
fprintf('VAL-MESH-DIFF-CRIT-001 passed.\n');
