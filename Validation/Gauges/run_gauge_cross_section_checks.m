function Summary = run_gauge_cross_section_checks()
%RUN_GAUGE_CROSS_SECTION_CHECKS Gauge discharge cross-sections: geometry, flux sum, full run.
%   1) Geometry on synthetic planes: auto direction, centred lines, even widths, edges.
%   2) Discharge from synthetic per-direction outflows (D4 and D8), with back flow.
%   3) Excel V-tilted run with four gauges: a full-width section one row above the outlet
%      carries 30/31 of the outlet volume (within 1 %), a full-width section carries at
%      least what the 5-cell channel section does, and an even width becomes odd.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(fullfile(root, 'HydroPol2D_Functions'), fullfile(here, '..', 'Groundwater', 'Layered'));
names = {}; details = {};
    function record(name, detail)
        names{end + 1, 1} = name; details{end + 1, 1} = detail; %#ok<AGROW>
        fprintf('[PASS] %s: %s\n', name, detail);
    end

% ---- 1) geometry ----------------------------------------------------------------
[cc, rr] = meshgrid(1:21, 1:15);
valid = true(15, 21);
east_plane = -cc;            % downhill to the east
south_plane = -rr;           % downhill to the south (rows grow southward)
S = hp2d_gauge_cross_sections(8, 11, 5, {'auto'}, east_plane, valid);
assert(strcmp(S.downstream, 'E') && isequal(S.rows', 6:10) && all(S.cols == 11), 'east plane: N-S line of 5');
S = hp2d_gauge_cross_sections(8, 11, 3, {'auto'}, south_plane, valid);
assert(strcmp(S.downstream, 'S') && all(S.rows == 8) && isequal(S.cols', 10:12), 'south plane: E-W line of 3');
warning('off', 'HydroPol2D:GaugeSection');
S = hp2d_gauge_cross_sections(8, 11, 4, {'N'}, south_plane, valid);
warning('on', 'HydroPol2D:GaugeSection');
assert(strcmp(S.downstream, 'N') && S.width == 5 && numel(S.cols) == 5 && ~S.auto, 'override and even width');
S = hp2d_gauge_cross_sections(2, 11, 7, {'auto'}, east_plane, valid);
assert(isequal(S.rows', 1:5), 'line clipped at the grid edge');
holes = valid; holes(6, 11) = false;
S = hp2d_gauge_cross_sections(8, 11, 5, {'auto'}, east_plane, holes);
assert(isequal(S.rows', [7 8 9 10]), 'cells outside the domain dropped');
record('geometry', 'auto E/S, centred lines, override, even width -> odd, edge clipping, holes');

% ---- 2) discharge ------------------------------------------------------------------
z = zeros(15, 21);
q = struct('left', z, 'right', 10 + z, 'up', z, 'down', z);
A = 400;                                         % 20 m cells
S = hp2d_gauge_cross_sections(8, 11, 5, {'E'}, east_plane, valid);
Q = hp2d_gauge_section_discharge(S, q, A);
assert(abs(Q - 5*10*A/3.6e6) < 1e-12, 'uniform eastward flow');
q.left(:, 12) = 4;                               % back flow from the next column
Q = hp2d_gauge_section_discharge(S, q, A);
assert(abs(Q - 5*(10 - 4)*A/3.6e6) < 1e-12, 'back flow subtracted');
qs = struct('left', z, 'right', z, 'up', z, 'down', 6 + z, 'ne', z, 'se', 1 + z, 'sw', 2 + z, 'nw', z);
S = hp2d_gauge_cross_sections(8, 11, 3, {'S'}, south_plane, valid);
Q = hp2d_gauge_section_discharge(S, qs, A);
assert(abs(Q - 3*(6 + 1 + 2)*A/3.6e6) < 1e-12, 'D8 diagonals crossing the line');
Sw = hp2d_gauge_cross_sections(8, 11, 3, {'W'}, south_plane, valid);
qw = struct('left', 3 + z, 'right', z, 'up', z, 'down', z);
assert(hp2d_gauge_section_discharge(Sw, qw, A) > 0, 'westward flow is positive downstream');
% Signed face fluxes (D4 inertial): F east through every east face, stored twice.
F = 7 + z; F(:, 12) = -2;                        % column 12 sends 2 mm/h back west
qf = struct('left', -[z(:, 1), F(:, 1:end - 1)], 'right', F, 'up', z, 'down', z, 'signed_faces', true);
S = hp2d_gauge_cross_sections(8, 11, 5, {'E'}, east_plane, valid);
assert(abs(hp2d_gauge_section_discharge(S, qf, A) - 5*7*A/3.6e6) < 1e-12, 'signed faces: no double count');
S = hp2d_gauge_cross_sections(8, 12, 5, {'E'}, east_plane, valid);
assert(abs(hp2d_gauge_section_discharge(S, qf, A) + 5*2*A/3.6e6) < 1e-12, 'signed faces: reverse flow is negative');
Fs = 7 + z;                                      % south through every south face
qn = struct('left', z, 'right', z, 'up', -[z(1, :); Fs(1:end - 1, :)], 'down', Fs, 'signed_faces', true);
S = hp2d_gauge_cross_sections(8, 11, 3, {'S'}, south_plane, valid);
assert(abs(hp2d_gauge_section_discharge(S, qn, A) - 3*7*A/3.6e6) < 1e-12, 'signed faces: southward');
record('discharge', 'uniform flow, back flow, D8 diagonals, signed faces, sign by direction');

% ---- 3) full Excel run -----------------------------------------------------------------
R = run_layered_vtilted_case('gauges');
% The local-inertial surface run closes to about 0.14 mm over the 508,400 m2 catchment,
% the same with gauges on or off (checked 2026-10-20): the sections only read the fluxes.
assert(abs(R.volume_error_m3)/508400*1000 < 0.5, 'volume error %.3e m3', R.volume_error_m3);
dirs = {R.gauge_sections.downstream};
assert(all(strcmp(dirs(1:3), 'S')), 'auto direction along the channel is south (got %s)', strjoin(dirs, ','));
assert(R.gauge_sections(4).width == 5, 'even width 4 must become 5');
t = R.time_hydrograph_min*60;
Qs = R.gauge_section_discharge_m3_s;
vol_outlet = trapz(t, R.outlet_hydrograph_m3_s);
vol_section = trapz(t, Qs(:, 2));
% The section (row 30) drains every row but the outlet row: with all 31 x 41 cells active
% and the same runoff everywhere, it carries 30/31 of the outlet volume.
ny = size(R.head_m, 1);
share = R.gauge_sections(2).rows(1)/ny;
rel = abs(vol_section - share*vol_outlet)/max(share*vol_outlet, eps);
assert(rel < 0.01, 'section above the outlet: volume %.4g vs %.4g m3 expected (%.2f %%)', ...
    vol_section, share*vol_outlet, 100*rel);
assert(trapz(t, Qs(:, 1)) >= trapz(t, Qs(:, 3)) - 1e-9, 'full width must carry at least the channel section');
record('V-tilted Excel run', sprintf(['volume error %.1f m3; section above the outlet %.4g m3 vs ', ...
    '%.4g m3 expected (30/31 of the outlet, %.2f %%); peaks %.3f vs %.3f m3/s'], R.volume_error_m3, ...
    vol_section, share*vol_outlet, 100*rel, max(Qs(:, 2)), max(R.outlet_hydrograph_m3_s)));
Summary = table(string(names), string(details), 'VariableNames', {'check', 'detail'});
fprintf('Gauge cross-section checks: %d of 3 passed\n', height(Summary));
end
