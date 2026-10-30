function [surface_volume, modules, diag] = hp2d_voronoi_modules_step(mesh, modules, surface_volume, edge_q, dt_s)
%HP2D_VORONOI_MODULES_STEP Advance the mesh modules after surface routing.
%   In the grid's order (HydroPol2D-Python mesh_modules.after_surface_step):
%   reservoir releases, then pollutant wash-off on the step's edge
%   discharges, then the human risk index on the new depths.
%   edge_q is the discharge per width [m2/s], positive from owner to
%   neighbour. diag.reservoir_outflow_m3 is the water that left the model
%   through a (0, 0) reservoir destination.

diag = struct('reservoir_outflow_m3', 0, 'pollutant_export_kg', 0);

if ~isempty(modules.reservoirs)
    R = modules.reservoirs;
    volume = max(surface_volume, 0);
    column = ones(size(R.source));
    xds1 = double(R.ds1 ~= 0); xds1(isnan(R.ds1)) = NaN;
    xds2 = double(R.ds2 ~= 0); xds2(isnan(R.ds2)) = NaN;
    [exchange, external, release] = hp2d_reservoir_releases(volume, column, R.source, ...
        R.k1, R.h1, R.k2, R.k3, R.h2, R.k4, R.ds1, xds1, R.ds2, xds2, dt_s, ...
        struct('type', 'area', 'area', mesh.surface_area(:)));
    surface_volume = max(volume + exchange, 0);
    diag.reservoir_outflow_m3 = sum(external);
    R.discharge_m3_s = release / dt_s;
    R.cumulative_outflow_m3 = R.cumulative_outflow_m3 + diag.reservoir_outflow_m3;
    modules.reservoirs = R;
end

if ~isempty(modules.water_quality)
    W = modules.water_quality;
    [W.mass_kg, exported, W.outlet_concentration_mg_l, W.cell_concentration_mg_l] = ...
        washoff_edges(W, edge_q(:) .* mesh.edge_length(:), mesh.edge_owner(:), ...
        mesh.edge_neighbor(:), dt_s);
    W.cumulative_export_kg = W.cumulative_export_kg + exported;
    diag.pollutant_export_kg = exported;
    modules.water_quality = W;
end

if ~isempty(modules.human_risk)
    H = modules.human_risk;
    P = H.parameters;
    depth = max(surface_volume ./ mesh.surface_area(:), 0);
    hu = H.discharge_x * edge_q(:);
    hv = H.discharge_y * edge_q(:);
    speed = zeros(size(depth));
    wet = depth > 1e-3;   % the grid's wet-cell threshold for the risk velocity
    speed(wet) = min(hypot(hu(wet), hv(wet)) ./ depth(wet), H.maximum_velocity_m_s);
    % flag_human_instability = 1 (Human_Instability_Module.m)
    F_person = P.weight_person * P.gravity;
    F_buoy = P.width1_person * P.width2_person * depth * P.ro_water * P.gravity;
    available_friction = P.mu * max(F_person - F_buoy, 0);
    hydro_force = max(0.5 * P.Cd * P.ro_water * P.width1_person * depth .* speed.^2, 0);
    risk = zeros(size(depth));
    supported = available_friction > 0;
    risk(supported) = hydro_force(supported) ./ available_friction(supported);
    risk(~supported & depth > 0) = 1;
    H.risk = min(max(risk, 0), 1);
    H.max_risk = max(H.risk, H.max_risk);
    modules.human_risk = H;
end
end

function [mass, exported, outlet_concentration, cell_concentration] = washoff_edges(W, q, owner, neighbor, dt_s)
% W_e = C Q_e^p B on every edge a cell drains through, kg/s (C = C_3/3600).
% Sub-steps n = ceil(max(gross rate / B) dt) of dt/n, and each cell's
% outflow is scaled so it never exceeds the mass it holds
% (mass_based_washoff_edge_step in HydroPol2D-Python). A boundary edge
% (neighbour 0) exports the pollutant out of the model.
n = numel(W.mass_kg);
mass = max(W.mass_kg, 0);
forward = q >= 0;
donor = owner; donor(~forward) = neighbor(~forward);
receiver = neighbor; receiver(~forward) = owner(~forward);
draining = donor > 0;   % an edge fed from outside brings water, not pollutant
donor(~draining) = 1;
q_out = abs(q); q_out(~draining) = 0;
leaving = draining & receiver == 0;
receiver(receiver == 0) = 1;
q_total = accumarray(donor, q_out, [n 1]);

[~, total] = edge_rates(mass);
ratio = zeros(n, 1);
ratio(mass > 0) = total(mass > 0) ./ max(mass(mass > 0), 1e-300);
substeps = max(1, ceil(max(ratio) * dt_s));
sub_dt = dt_s / substeps;
exported = 0;
load = zeros(n, 1);
for k = 1:substeps
    [rate, total] = edge_rates(mass);
    scale = ones(n, 1);
    over = total * sub_dt > mass;
    scale(over) = mass(over) ./ max(total(over) * sub_dt, 1e-30);
    rate = rate .* scale(donor);
    load = accumarray(donor, rate, [n 1]);
    incoming = accumarray(receiver, rate .* ~leaving, [n 1]);
    exported = exported + sum(rate(leaving)) * sub_dt;
    mass = max(mass + (incoming - load) * sub_dt, 0);
end
outlet_water = sum(q_out(leaving)) * dt_s;
outlet_concentration = 0;
if outlet_water > 0, outlet_concentration = 1000 * exported / outlet_water; end
cell_concentration = zeros(n, 1);
cell_concentration(q_total > 0) = 1000 * load(q_total > 0) ./ q_total(q_total > 0);

    function [rate, total] = edge_rates(current)
        active = current .* (current >= W.threshold_kg);
        rate = W.coefficient_kg_s(donor) .* active(donor) .* q_out .^ W.exponent(donor);
        rate(~draining) = 0;
        total = accumarray(donor, rate, [n 1]);
    end
end
