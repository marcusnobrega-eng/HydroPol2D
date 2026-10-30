% VAL-WQ-SUBSTEP-001: wash-off sub-steps sum to the step; outflow never
% exceeds the mass held.
%
% build_up_wash_off split a step it could not take in one go into
% floor(T/tmin)+1 pieces, the last of length (T/tmin - floor(T/tmin)) tmin.
% With T = 60 s and tmin = 25 s that integrates 2 x 20 s + 10 s = 50 s, not
% 60 s. It had no limiter either, so a cell could release more than it held
% and the negative mass was booked as "lost". The model time step was also
% capped at tmin_wq. HydroPol2D-Python uses n = ceil(max(W_tot/B) dt) equal
% sub-steps, scales each cell's outflow to the mass it holds, and
% sub-cycles inside the hydraulic step (modules/water_quality.py).
% The check runs one 120 s step on a 3 x 4 grid with flows in all four
% directions and out of every side. It compares the masses and the export
% with Python's mass_based_washoff_step on the same numbers
% (python_reference.csv, made by make_python_reference.py from
% washoff_inputs.csv), and closes the pollutant mass.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

ny = 3; nx = 4; area = 100; dt_s = 120; C_4 = 1.3; min_Bt = 0.5;   % g/m2
[jj, ii] = meshgrid(1:nx, 1:ny);
mass = 0.5 + 0.3 * mod(ii * 3 + jj * 5, 7);                        % kg
mass(2, 2) = 0.02;                                                % below min_Bt (0.05 kg)
C_3 = 3e4 + 1e4 * mod(ii + 2 * jj, 3);                         % kg/h per (m3/s)^C_4 per kg: many sub-steps
qx = 0.03 * sin(1.7 * (1:ny)' + 0.9 * (0:nx));                     % m3/s, (ny, nx+1), + east
qy = 0.02 * cos(1.1 * (0:ny)' + 1.3 * (1:nx));                     % m3/s, (ny+1, nx), + south
to_mm_h = @(q) q / area * 1000 * 3600;
qout_left = to_mm_h(max(-qx(:, 1:end-1), 0));
qout_right = to_mm_h(max(qx(:, 2:end), 0));
qout_up = to_mm_h(max(-qy(1:end-1, :), 0));
qout_down = to_mm_h(max(qy(2:end, :), 0));

fid = fopen(fullfile(case_dir, 'washoff_inputs.csv'), 'w');
fprintf(fid, 'name,values\n');
write_row(fid, 'mass_kg', reshape(mass', 1, []));
write_row(fid, 'qx_m3_s', reshape(qx', 1, []));
write_row(fid, 'qy_m3_s', reshape(qy', 1, []));
write_row(fid, 'C_3', reshape(C_3', 1, []));
write_row(fid, 'C_4', C_4); write_row(fid, 'min_Bt', min_Bt);
write_row(fid, 'cell_area_m2', area); write_row(fid, 'dt_s', dt_s);
fclose(fid);

idx_nan_5 = false(ny, nx, 5);
[B, ~, ~, tmin_wq, ~, mass_lost] = build_up_wash_off(C_3, C_4, qout_left, qout_right, ...
    qout_up, qout_down, zeros(ny, nx), mass, dt_s / 60, nx, ny, area, zeros(ny, nx), ...
    idx_nan_5, 0, 0, zeros(ny, nx), 0, 1e9, min_Bt);

reference = read_rows(fullfile(case_dir, 'python_reference.csv'));
python_mass = reshape(reference.mass_kg, nx, ny)';
python_export = reference.export_kg;
mass_error = max(abs(B - python_mass), [], 'all');
export_error = abs((sum(mass, 'all') - sum(B, 'all')) - python_export);
passed = mass_error < 1e-12 && export_error < 1e-12 && mass_lost == 0 && all(B(:) >= 0) && ...
    tmin_wq == dt_s;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-WQ-SUBSTEP-001", "Wash-off sub-steps and limiter", mass_error, ...
    export_error, mass_lost, tmin_wq, status, true, ...
    'VariableNames', {'case_id','case_name','max_mass_error_vs_python_kg', ...
    'export_error_vs_python_kg','mass_lost_kg','tmin_wq_s','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-WQ-SUBSTEP-001 failed.');
end
fprintf('VAL-WQ-SUBSTEP-001 passed.\n');

function rows = read_rows(path)
% name,values lines; values are space-separated numbers.
lines = splitlines(strtrim(string(fileread(path))));
rows = struct();
for k = 2:numel(lines)
    parts = split(lines(k), ',');
    rows.(char(parts(1))) = str2double(split(strtrim(parts(2))))';
end
end

function write_row(fid, name, values)
fprintf(fid, '%s,%s\n', name, strjoin(compose('%.17g', values), ' '));
end
