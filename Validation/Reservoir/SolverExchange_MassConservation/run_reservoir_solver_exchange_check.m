% VAL-RES-002: reservoir releases inside the routing solvers.
%
% Each surface solver advances one 1 s step on a 5 x 11 grid of 10 m cells.
% Six pits sit among cells outside the domain (no data), so no water crosses a
% face and only the reservoir controls move water:
%   control 1: source pit (3,3), 2 m deep
%       outlet 1: k = 2, sill 0.5 m, b = 1.5 -> pit (3,5)
%       outlet 2: k = 1, sill 1.0 m, b = 1.5 -> pit (3,7)
%   control 2: source pit (3,9), 1.5 m deep
%       outlet 1: k = 3, sill 0.2 m, b = 1.5 -> out of the model (0, 0)
%   control 3: source pit (5,3), 1.2 m deep, one outlet
%       outlet 1: k = 1.5, sill 0.3 m, b = 1.5 -> pit (5,5)
% Reference (reservoir_outlet_transfers in HydroPol2D-Python): the outlets
% release in order, Q = k max(h - h0, 0)^b, capped at the volume above the
% sill, with h updated after outlet 1. The release leaves the source and
% enters the destination; the external release is booked as outlet flow.
% The Neal run gives the source pit a 2 m wide, 1 m deep channel.
% Before the fix the default, Neal, full-momentum and wave solvers added the
% release to the source (water created), outlet 1 reached its destination
% twice when outlet 2 existed, CA lost outlet 1 without an outlet 2, and a
% destination of (0, 0) was an invalid index.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

addpath(functions_dir, '-begin');

solvers = ["local_inertial", "neal", "full_momentum", "diffusive", "kinematic", "ca"];
Results = table();
for s = solvers
    Results = [Results; run_solver(s)]; %#ok<AGROW>
end

writetable(Results, fullfile(out_dir, 'Reservoir_Solver_Exchange.csv'));
passed = all(Results.passed);
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-RES-002", "Reservoir releases inside the routing solvers", ...
    status, true, 'VariableNames', {'case_id','case_name','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));

disp(Results);
if ~passed
    error('VAL-RES-002 failed.');
end
fprintf('VAL-RES-002 passed.\n');

function Row = run_solver(solver)
ny = 5; nx = 11;
dx = 10;
area = dx^2;
dt_s = 1;
time_step = dt_s / 60;  % [min]

z = NaN(ny, nx);
pits = [3 3; 3 5; 3 7; 3 9; 5 3; 5 5];
for p = 1:size(pits, 1)
    z(pits(p,1), pits(p,2)) = 0;
end
h0 = zeros(ny, nx);
h0(3,3) = 2.0;
h0(3,9) = 1.5;
h0(5,3) = 1.2;

R.x = [3; 9; 3];  R.y = [3; 3; 5];
R.k1 = [2; 3; 1.5]; R.h1 = [0.5; 0.2; 0.3]; R.k2 = [1.5; 1.5; 1.5];
R.xds1 = [5; 0; 5]; R.yds1 = [3; 0; 5];
R.k3 = [1; 0; 0]; R.h2 = [1.0; 0; 0]; R.k4 = [1.5; 1.5; 1.5];
R.xds2 = [7; NaN; NaN]; R.yds2 = [3; NaN; NaN];

River_Width = zeros(ny, nx);
River_Depth = zeros(ny, nx);
flag_subgrid = 0; flag_overbanks = 0;
if solver == "neal"
    River_Width(3,3) = 2;
    River_Depth(3,3) = 1;
    flag_subgrid = 1; flag_overbanks = 1;
end

% Reference volumes and releases.
V0 = area * h0;
if solver == "neal"
    V0(3,3) = hp2d_neal_cell_volume(h0(3,3), 2, 1, dx);
end
[V_ref, external_ref] = reference_releases(V0, h0, R, dt_s, solver == "neal", dx);

% Solver inputs.
idx_nan = isnan(z);
d_tot = 1000 * h0;
d_tot(idx_nan) = NaN;
d_p = d_tot;
n = 0.03 * ones(ny, nx);
outlet_index = false(ny, nx);
row_outlet = []; col_outlet = [];
outflow = zeros(ny, nx, 5);
args_res = {R.x, R.y, R.k1, R.h1, R.k2, R.k3, R.h2, R.k4, R.yds1, R.xds1, R.yds2, R.xds2, 1};

switch solver
    case {"local_inertial", "neal"}
        [~,~,~,~,outlet_flow,d_t,~,~,~,~,~,~,~,C_a] = Local_Inertial_Model_D4(1, args_res{:}, ...
            z, d_tot, d_p, n, n.^2, area, time_step, dx, outlet_index, 1, 0.01, ...
            row_outlet, col_outlet, 1e-3, outflow(:,:,1:3), idx_nan, 0, flag_subgrid, ...
            n, n, River_Width, River_Depth, [], [], [], [], area, struct(), ...
            flag_overbanks, 0, []);
    case "full_momentum"
        [~,~,~,~,outlet_flow,d_t,~,~,~,~,~,~,~,C_a] = Full_Momentum_Model_D4(1, args_res{:}, ...
            z, d_tot, d_p, n, n.^2, area, time_step, dx, outlet_index, 1, 0.01, ...
            row_outlet, col_outlet, 1e-3, outflow, idx_nan, 0, 0, ...
            n, n, River_Width, River_Depth, [], [], [], [], area, struct(), ...
            0, 0, [], struct());
    case "diffusive"
        [~,~,~,~,outlet_flow,d_t,~,~,~,~,~,~,~,C_a] = Diffusive_Wave_Model_Implicit(1, args_res{:}, ...
            z, d_tot, d_p, n, area, time_step, dx, outlet_index, 1, 0.01, ...
            row_outlet, col_outlet, 1e-3, outflow(:,:,1:3), idx_nan, 0, 0, ...
            n, n, River_Width, River_Depth, [], [], [], [], area);
    case "kinematic"
        [~,~,~,~,outlet_flow,d_t,~,~,~,~,~,~,~,C_a] = Kinematic_Wave_Model(1, args_res{:}, ...
            z, d_tot, d_p, n, area, time_step, dx, outlet_index, 1, 0.01, ...
            row_outlet, col_outlet, 1e-3, outflow(:,:,1:3), idx_nan, 0, 0, ...
            n, n, River_Width, River_Depth, [], [], [], [], area);
    case "ca"
        [~,~,~,~,outlet_flow,d_t,~] = CA_Routing(args_res{:}, ...
            z, d_tot, n, area, time_step, zeros(ny, nx), dx, zeros(ny, nx), ...
            outlet_index, 1, 0.01, row_outlet, col_outlet, idx_nan, 0);
        C_a = area;
end

% Model volumes after the step.
V_model = area * d_t / 1000;
V_model(idx_nan) = 0;
if solver == "neal"
    V_model(3,3) = hp2d_neal_cell_volume(d_t(3,3) / 1000, 2, 1, dx);
end
if isscalar(C_a); C_a = C_a * ones(ny, nx); end %#ok<NASGU>
external_model = sum(outlet_flow(:), 'omitnan') / 1000 / 3600 * area * dt_s;

cell_error = max(abs(V_model(:) - V_ref(:)));
external_error = abs(external_model - external_ref);
ledger_residual = sum(V_model(:)) + external_model - sum(V0(:));
passed = cell_error < 1e-9 && external_error < 1e-9 && abs(ledger_residual) < 1e-9;

Row = table(solver, V_ref(3,5), V_model(3,5), V_ref(3,7), V_model(3,7), ...
    V_ref(5,5), V_model(5,5), V_ref(3,3), V_model(3,3), external_ref, external_model, cell_error, ...
    ledger_residual, passed, ...
    'VariableNames', {'solver','dest1_ref_m3','dest1_model_m3','dest2_ref_m3', ...
    'dest2_model_m3','single_dest_ref_m3','single_dest_model_m3', ...
    'source_ref_m3','source_model_m3','external_ref_m3', ...
    'external_model_m3','max_cell_error_m3','ledger_residual_m3','passed'});
end

function [V, external] = reference_releases(V, h0, R, dt_s, neal, dx)
external = 0;
for ii = 1:numel(R.y)
    for outlet = 1:2
        if outlet == 1
            k = R.k1(ii); sill = R.h1(ii); b = R.k2(ii); row = R.yds1(ii); col = R.xds1(ii);
        else
            k = R.k3(ii); sill = R.h2(ii); b = R.k4(ii); row = R.yds2(ii); col = R.xds2(ii);
        end
        if isnan(row); continue; end
        r = R.y(ii); c = R.x(ii);
        if neal && r == 3 && c == 3
            depth = hp2d_neal_depth_from_volume(V(r,c), 2, 1, dx);
            above_sill = V(r,c) - hp2d_neal_cell_volume(sill, 2, 1, dx);
        else
            depth = V(r,c) / dx^2;
            above_sill = (depth - sill) * dx^2;
        end
        dV = min(k * max(depth - sill, 0)^b * dt_s, max(above_sill, 0));
        V(r,c) = V(r,c) - dV;
        if row == 0 && col == 0
            external = external + dV;
        else
            V(row,col) = V(row,col) + dV;
        end
    end
end
%#ok<*INUSL>
end
