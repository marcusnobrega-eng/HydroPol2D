function metrics = run_structured_lookup_roughness_validation()
%RUN_STRUCTURED_LOOKUP_ROUGHNESS_VALIDATION Guard the g*n^2 table convention.

g = 9.81;
n_expected = 0.04;
resolution = 10;
depth = 1;

S = struct();
S.sfincs_exact = true;
S.roughness_storage = 'g_n_squared';
S.g = g;
S.dz = 1;
S.maxDepth = 1;
S.area_x = reshape([resolution * depth, resolution * depth], 1, 1, []);
S.width_x = reshape([resolution, resolution], 1, 1, []);
S.wetfrac_x = reshape([1, 1], 1, 1, []);
S.hrep_x = reshape([depth, depth], 1, 1, []);
S.phi_x = reshape([1, 1], 1, 1, []);
S.perimeter_x = reshape([resolution + 2 * depth, resolution + 2 * depth], 1, 1, []);
S.Rh_x = S.area_x ./ S.perimeter_x;
S.nrep_x = reshape([g * n_expected^2, g * n_expected^2], 1, 1, []);
S.n_x = S.nrep_x;
S.K_x = [];
S.invert_x = 0;

face = hp2d_subgrid_face_state(S, depth, 'x', resolution);
n_actual = double(face.n);
assert(abs(n_actual - n_expected) < 1e-12, ...
    'Structured lookup returned %.12g instead of Manning n %.12g.', ...
    n_actual, n_expected);

rh = resolution * depth / (resolution + 2 * depth);
k_expected = resolution * depth * rh^(2/3) / n_expected;
assert(abs(double(face.K) - k_expected) < 1e-10, ...
    'Fallback conveyance does not use the physical Manning coefficient.');

metrics = struct('manning_expected', n_expected, ...
    'manning_returned', n_actual, 'conveyance_expected', k_expected, ...
    'conveyance_returned', double(face.K), 'passed', true);
fprintf('Structured lookup roughness validation passed: n = %.4f.\n', n_actual);
end
