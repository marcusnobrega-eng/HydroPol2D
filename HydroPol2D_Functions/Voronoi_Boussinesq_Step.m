function [head, diagnostics] = Voronoi_Boussinesq_Step(mesh, head, aquifer_bottom, hydraulic_conductivity, specific_yield, dt)
%VORONOI_BOUSSINESQ_STEP Conservative unconfined groundwater edge flux.

n = mesh.n_cells;
head = expand(head, n); aquifer_bottom = expand(aquifer_bottom, n);
hydraulic_conductivity = expand(hydraulic_conductivity, n); specific_yield = expand(specific_yield, n);
assert(all(hydraulic_conductivity >= 0 & specific_yield > 0));
thickness = max(head - aquifer_bottom, 0);
storage = specific_yield .* thickness .* mesh.cell_area(:);
internal = mesh.edge_neighbor > 0;
o = mesh.edge_owner(internal); d = mesh.edge_neighbor(internal);

% Follow the Python reference: arithmetic face properties and a
% predictor-corrector evaluation of the unconfined Boussinesq flux.
Q1 = edge_flux(head,aquifer_bottom,hydraulic_conductivity,mesh,o,d,internal);
delta1 = conservative_delta(Q1,o,d,n,dt,storage);
predicted_storage = max(storage + delta1,0);
predicted_head = aquifer_bottom + predicted_storage ./ ...
    (specific_yield .* mesh.cell_area(:));
Q2 = edge_flux(predicted_head,aquifer_bottom,hydraulic_conductivity,mesh,o,d,internal);
Q = 0.5 .* (Q1 + Q2);
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

function Q = edge_flux(head,bottom,conductivity,mesh,o,d,internal)
thickness=max(head-bottom,0);
face_thickness=0.5.*(thickness(o)+thickness(d));
face_conductivity=0.5.*(conductivity(o)+conductivity(d));
Q=-face_conductivity.*face_thickness.*mesh.edge_length(internal).* ...
    (head(d)-head(o))./mesh.edge_distance(internal);
end

function delta = conservative_delta(Q,o,d,n,dt,storage)
out_cell=o; out_cell(Q<0)=d(Q<0);
requested=accumarray(out_cell,abs(Q).*dt,[n 1],@sum,0);
scale=min(1,storage./max(requested,eps));
limited=Q.*scale(out_cell);
delta=accumarray(o,-limited.*dt,[n 1],@sum,0)+ ...
    accumarray(d,limited.*dt,[n 1],@sum,0);
end

function value = expand(value, n)
if isscalar(value), value = repmat(double(value), n, 1); else, value = double(value(:)); end
assert(numel(value) == n && all(isfinite(value)));
end
