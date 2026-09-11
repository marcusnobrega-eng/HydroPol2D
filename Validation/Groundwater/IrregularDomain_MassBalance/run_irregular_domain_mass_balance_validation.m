clear; clc;

model_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(model_root, 'HydroPol2D_Functions'));

nrows = 31;
ncols = 37;
[x, y] = meshgrid(1:ncols, 1:nrows);
domain = ((x - 19) / 16).^2 + ((y - 16) / 13).^2 <= 1;
domain(13:18, 16:21) = false;

z0 = zeros(nrows, ncols);
h = 5 + 0.01 * x - 0.005 * y;
K = 1e-6 * ones(nrows, ncols);
Sy = 0.2 * ones(nrows, ncols);
R = zeros(nrows, ncols);
h_soil = 100 * ones(nrows, ncols);

z0(~domain) = nan;
h(~domain) = nan;
K(~domain) = nan;
Sy(~domain) = nan;
h_soil(~domain) = nan;

empty_mask = false(nrows, ncols);
nan_head = nan(nrows, ncols);
storage_before = sum(Sy(domain) .* (h(domain) - z0(domain))) * 2000^2;

[h_after, ~, ~, q_exf, ~, balance_error] = Boussinesq_2D_explicit( ...
    3600, 2000, 2000, h, z0, Sy, R, K, empty_mask, 0, nan_head, ...
    nan_head, 0.5, h_soil, domain, empty_mask, nan_head, empty_mask, ...
    nan_head, empty_mask);

storage_after = sum(Sy(domain) .* (h_after(domain) - z0(domain))) * 2000^2;
relative_error = abs(storage_after - storage_before) / storage_before;

assert(all(isfinite(h_after(domain))), 'Valid-domain heads became non-finite.');
assert(all(isnan(h_after(~domain))), 'Heads outside the domain must remain NaN.');
assert(all(h_after(domain) >= z0(domain)), 'Groundwater head fell below bedrock.');
assert(all(q_exf(domain) == 0), 'Unexpected exfiltration occurred.');
assert(relative_error < 1e-10, 'Closed-domain storage was not conserved.');
assert(abs(balance_error) / storage_before < 1e-10, 'Reported mass balance did not close.');

fprintf('PASS: irregular-domain groundwater balance; relative error = %.3e\n', relative_error);
