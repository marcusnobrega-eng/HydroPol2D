function [head, diagnostics] = Voronoi_Boussinesq_Step(mesh, head, aquifer_bottom, hydraulic_conductivity, specific_yield, dt)
%VORONOI_BOUSSINESQ_STEP Conservative unconfined groundwater edge flux.

n = mesh.n_cells;
head = expand(head, n); aquifer_bottom = expand(aquifer_bottom, n);
hydraulic_conductivity = expand(hydraulic_conductivity, n); specific_yield = expand(specific_yield, n);
assert(all(hydraulic_conductivity >= 0 & specific_yield > 0));
thickness = max(head - aquifer_bottom, 0);
storage = specific_yield .* thickness .* mesh.cell_area(:);
transmissivity = hydraulic_conductivity .* thickness;
internal = mesh.edge_neighbor > 0;
o = mesh.edge_owner(internal); d = mesh.edge_neighbor(internal);
To = transmissivity(o); Td = transmissivity(d);
Tface = 2 .* To .* Td ./ max(To + Td, eps);
Q = -Tface .* mesh.edge_length(internal) .* (head(d) - head(o)) ./ mesh.edge_distance(internal);
out_cell = o; out_cell(Q < 0) = d(Q < 0);
requested = accumarray(out_cell, abs(Q) .* dt, [n 1], @sum, 0);
scale = min(1, storage ./ max(requested, eps));
Q = Q .* scale(out_cell);
delta = accumarray(o, -Q .* dt, [n 1], @sum, 0) + accumarray(d, Q .* dt, [n 1], @sum, 0);
storage = max(storage + delta, 0);
head = aquifer_bottom + storage ./ (specific_yield .* mesh.cell_area(:));
diagnostics.max_flux_m3_s = max(abs(Q), [], 'omitnan');
diagnostics.mass_change_m3 = sum(delta);
end

function value = expand(value, n)
if isscalar(value), value = repmat(double(value), n, 1); else, value = double(value(:)); end
assert(numel(value) == n && all(isfinite(value)));
end
