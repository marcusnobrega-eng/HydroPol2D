function [S, F, host_m3] = hp2d_green_infrastructure_step(G, S, host_m3, rain_rate_mm_h, ...
    et_rate_mm_day, apply_et_crop_coefficient, dt_s)
%HP2D_GREEN_INFRASTRUCTURE_STEP Advance the green-infrastructure (GI) compartments one step.
%   Port of green_infrastructure_step in HydroPol2D-Python
%   (modules/green_infrastructure.py). Each compartment is the part of one GI
%   facility inside one grid cell, with three stores [m3]: surface ponding,
%   engineered media water and gravel water. In order:
%     1. rain on the compartment area;
%     2. capture from the host cell, Q = C W h sqrt(2 g dH), scaled so the
%        compartments of a host cell never take more than it holds;
%     3. weir overflow back to the host cell, Q = Cw W (excess depth)^1.5;
%     4. media infiltration, Ks A (1 + (h + psi(Se) dtheta) / L_media),
%        capped by the ponded water and the media deficit;
%     5. Mualem percolation to the gravel, Ks Kr(Se) A;
%     6. ET from the ponding, then from the media above residual;
%     7. lateral flow between compartments of a facility (D4): local-inertial
%        surface flow and Darcy gravel flow, in Courant substeps;
%     8. native exfiltration request, min(gravel, Ks_native A dt);
%     9. power-law underdrain, Q = Cd max(h_gravel - sill, 0)^b, to a host
%        surface cell or out of the model.
%   Every transfer is capped by what its donor holds, so the step conserves
%   mass. Unlike Python (review C14), a donor with several lateral edges
%   cannot give more than it holds in total, and a gravel store cannot be
%   filled above its capacity by lateral inflow: Python caps each edge on
%   its own and then clips the store, which creates or deletes water.
%
%   G: prepared geometry and parameters (hp2d_prepare_green_infrastructure).
%   S: state struct with surface_m3, media_m3, gravel_m3 [n x 1] and
%      edge_q_m3_s [n_edges x 1] (surface-flow memory).
%   host_m3: water volume of the host cells G.host_cells [n_host x 1].
%   rain_rate_mm_h, et_rate_mm_day: forcing per compartment [n x 1].
%   apply_et_crop_coefficient: true when et_rate is a reference ET (the
%      design crop coefficient applies); false for a prescribed actual ET.
%   F: fluxes of the step [m3 per compartment unless noted]: rain, capture,
%      infiltration, overflow, percolation, et, native_request, underdrain,
%      facility_underdrain [n_facilities x 1], external_m3 (scalar, leaves the
%      model), surface_to_host [n_host x 1] (underdrain released to surface
%      targets), gravel_edge_q_m3_s, substeps.

g = 9.80665;
EPS = 1.0e-12;

area = G.area_m2;
bottom = G.bottom_elevation_m;
terrain = G.terrain_elevation_m;
media_depth = G.media_depth_m;
media_porosity = G.media_porosity;
media_residual = G.media_residual_water_content;
gravel_depth = G.gravel_depth_m;
gravel_void = G.gravel_void_fraction;
media_capacity = area .* media_depth .* media_porosity;
media_residual_volume = area .* media_depth .* media_residual;
gravel_capacity = area .* gravel_depth .* gravel_void;
n_host = numel(host_m3);
host_m3 = double(host_m3(:));
host = G.compartment_host(:);

surface = max(double(S.surface_m3(:)), 0);
media = min(max(double(S.media_m3(:)), 0), media_capacity);
gravel = min(max(double(S.gravel_m3(:)), 0), gravel_capacity);

% 1. Rain on the compartment.
rain = max(double(rain_rate_mm_h(:)), 0) / 1000 / 3600 .* area * dt_s;
surface = surface + rain;

% 2. Capture from the host cell (head above the GI water or the terrain).
host_depth = host_m3(host) / G.cell_area_m2;
gi_depth = surface ./ area;
host_head = terrain + host_depth;
gi_head = bottom + gi_depth;
head_drop = max(host_head - max(gi_head, terrain), 0);
capture_q = G.surface_exchange_coefficient .* G.inlet_width_m .* max(host_depth, 0) ...
    .* sqrt(2 * g * head_drop);
proposed_capture = capture_q * dt_s;
requested_by_host = accumarray(host, proposed_capture, [n_host 1]);
capture_scale = min(1, host_m3 ./ max(requested_by_host, EPS));
capture = proposed_capture .* capture_scale(host);
host_m3 = host_m3 - accumarray(host, capture, [n_host 1]);
surface = surface + capture;

% 3. Weir overflow above the crest, back to the host cell.
crest_capacity = area .* G.overflow_crest_depth_m;
excess_depth = max((surface - crest_capacity) ./ area, 0);
overflow_q = G.overflow_coefficient .* G.inlet_width_m .* excess_depth .^ 1.5;
overflow = min(overflow_q * dt_s, max(surface - crest_capacity, 0));
surface = surface - overflow;
host_m3 = host_m3 + accumarray(host, overflow, [n_host 1]);

% 4. Media infiltration. As in Python, the ponded depth is the one after the
% rain and before capture (gi_depth above).
media_deficit = max(media_capacity - media, 0);
Se = min(max((media - media_residual_volume) ./ ...
    max(media_capacity - media_residual_volume, EPS), 1.0e-6), 1);
media_n = G.media_n_vg;
media_m = 1 - 1 ./ media_n;
suction_head = (1 ./ G.media_alpha_vg_m_inverse) .* (Se .^ (-1 ./ media_m) - 1) .^ (1 ./ media_n);
moisture_deficit = max((media_porosity - media_residual) .* (1 - Se), 0);
infiltration_capacity = G.media_ksat_m_s .* area .* ...
    (1 + (gi_depth + suction_head .* moisture_deficit) ./ max(media_depth, 1e-6));
infiltration = min(surface, infiltration_capacity * dt_s);
infiltration = min(infiltration, media_deficit);
surface = surface - infiltration;
media = media + infiltration;

% 5. Mualem percolation to the gravel.
Se = min(max((media - media_residual_volume) ./ ...
    max(media_capacity - media_residual_volume, EPS), 1.0e-6), 1);
conductivity_term = min(max(1 - Se .^ (1 ./ media_m), 0), 1);
relative_conductivity = Se .^ 0.5 .* (1 - conductivity_term .^ media_m) .^ 2;
percolation = G.media_ksat_m_s .* relative_conductivity .* area * dt_s;
percolation = min(percolation, max(media - media_residual_volume, 0));
percolation = min(percolation, max(gravel_capacity - gravel, 0));
media = media - percolation;
gravel = gravel + percolation;

% 6. ET: ponded water first, then media water above residual.
if apply_et_crop_coefficient
    et_factor = G.et_crop_coefficient;
else
    et_factor = ones(size(area));
end
et_demand = max(double(et_rate_mm_day(:)), 0) / 1000 / 86400 .* et_factor .* area * dt_s;
surface_et = min(surface, et_demand);
media_et = min(max(media - media_residual_volume, 0), max(et_demand - surface_et, 0));
surface = surface - surface_et;
media = media - media_et;
et_loss = surface_et + media_et;

% 7. Lateral flow between the compartments of each facility.
edge_q = double(S.edge_q_m3_s(:));
gravel_q = zeros(size(edge_q));
substeps = 0;
if G.connection_count > 0
    s = G.edge_start;
    e = G.edge_end;
    width = G.edge_width_m;
    len = G.edge_length_m;
    n_face = 0.5 * (G.surface_manning_n(s) + G.surface_manning_n(e));
    gravel_k = 0.5 * (G.gravel_ksat_m_s(s) + G.gravel_ksat_m_s(e));
    cross_section = width .* min(gravel_depth(s) .* gravel_void(s), gravel_depth(e) .* gravel_void(e));
    maximum_depth = max(max(surface ./ area, 0));
    stable_dt = G.courant * min(len) / max(sqrt(g * max(maximum_depth, 1.0e-6)), 1.0e-6);
    substeps = min(G.maximum_substeps, max(1, ceil(dt_s / max(stable_dt, G.minimum_substep_s))));
    sub_dt = dt_s / substeps;
    for sub = 1:substeps %#ok<FXUP>
        depth = surface ./ area;
        water_surface = bottom + depth;
        edge_depth = max((depth(s) + depth(e)) * 0.5, 0);
        flow_area = width .* edge_depth;
        radius = flow_area ./ max(width + 2 * edge_depth, EPS);
        slope = (water_surface(s) - water_surface(e)) ./ len;
        numerator = edge_q + g * flow_area * sub_dt .* slope;
        denominator = 1 + g * n_face .^ 2 .* abs(edge_q) * sub_dt ./ ...
            max(flow_area .* radius .^ (4 / 3), EPS);
        proposed = numerator ./ denominator;
        proposed(~(flow_area > EPS)) = 0;
        [surface, edge_q] = limited_edge_transfer(surface, s, e, proposed, sub_dt, []);

        gravel_fraction = gravel ./ max(gravel_capacity, EPS);
        gravel_head = bottom - media_depth + gravel_fraction .* gravel_depth;
        proposed_gravel = gravel_k .* cross_section .* (gravel_head(s) - gravel_head(e)) ./ len;
        [gravel, gravel_q] = limited_edge_transfer(gravel, s, e, proposed_gravel, sub_dt, ...
            gravel_capacity);
    end
end

% 8. Native exfiltration request into the soil below (the caller books
% what the soil accepts and returns the rest to the gravel).
native_request = min(gravel, G.native_ksat_m_s .* area * dt_s);
gravel = gravel - native_request;

% 9. Underdrain.
gravel_head_depth = gravel ./ max(area .* gravel_void, EPS);
compartment_target_kind = G.target_kind(G.facility_index);
underdrain_q = G.underdrain_coefficient .* ...
    max(gravel_head_depth - G.underdrain_sill_depth_m, 0) .^ G.underdrain_exponent;
underdrain_q(compartment_target_kind == 0) = 0;
underdrain = min(underdrain_q * dt_s, gravel);
gravel = gravel - underdrain;
facility_underdrain = accumarray(G.facility_index(:), underdrain, [G.facility_count 1]);
surface_to_host = zeros(n_host, 1);
external_m3 = 0;
for facility = 1:G.facility_count
    released = facility_underdrain(facility);
    switch G.target_kind(facility)
        case 1   % host surface cell
            target = G.target_host(facility);
            host_m3(target) = host_m3(target) + released;
            surface_to_host(target) = surface_to_host(target) + released;
        case 3   % out of the model
            external_m3 = external_m3 + released;
    end
end

S.surface_m3 = surface;
S.media_m3 = media;
S.gravel_m3 = gravel;
S.edge_q_m3_s = edge_q;

F = struct();
F.rain = rain;
F.capture = capture;
F.overflow = overflow;
F.infiltration = infiltration;
F.percolation = percolation;
F.et = et_loss;
F.native_request = native_request;
F.underdrain = underdrain;
F.facility_underdrain = facility_underdrain;
F.external_m3 = external_m3;
F.surface_to_host = surface_to_host;
F.gravel_edge_q_m3_s = gravel_q;
F.substeps = substeps;
end

function [V, q] = limited_edge_transfer(V, s, e, proposed_q, dt_s, capacity)
%LIMITED_EDGE_TRANSFER Move proposed_q * dt along each edge, capped so no store
%   gives more than it holds (all its outflows together) and, with capacity,
%   no store receives more than its free capacity. Edges flow from s to e
%   when proposed_q >= 0 (Python's _limited_edge_transfer caps each edge on
%   its own).
n = numel(V);
reverse = proposed_q < 0;
donor = s;
donor(reverse) = e(reverse);
receiver = e;
receiver(reverse) = s(reverse);
amount = abs(proposed_q) * dt_s;
outflow = accumarray(donor, amount, [n 1]);
donor_scale = ones(n, 1);
has = outflow > 0;
donor_scale(has) = min(1, max(V(has), 0) ./ outflow(has));
amount = amount .* donor_scale(donor);
if ~isempty(capacity)
    inflow = accumarray(receiver, amount, [n 1]);
    receiver_scale = ones(n, 1);
    has = inflow > 0;
    receiver_scale(has) = min(1, max(capacity(has) - V(has), 0) ./ inflow(has));
    amount = amount .* receiver_scale(receiver);
end
transfer = amount;
transfer(reverse) = -transfer(reverse);
V = V - accumarray(s, transfer, [n 1]) + accumarray(e, transfer, [n 1]);
q = transfer / dt_s;
end
