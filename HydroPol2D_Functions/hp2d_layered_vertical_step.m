function S = hp2d_layered_vertical_step(column, S, P, dt_s)
%HP2D_LAYERED_VERTICAL_STEP Column processes of one layered-subsurface step (cells x slots).
%   S (state, [N x 1] unless noted): water_mm [N x L], head_m, surface_depth_mm,
%   open_water_evaporation_mm (demand over the step), soil_et_demand_mm (demand over
%   the step, before the pervious fraction).
%   P (parameters): valid, pervious, root_depth_m, ksat_near_surface_mm_h, l_vg,
%   ltop_m, dh_max_m, ksurf_weight, feddes ([] = off), drainage_substeps, and the flags
%   infiltration, evapotranspiration, groundwater.
%   Order (HydroPol2D-Python advance_layered_subsurface): infiltration into the soil
%   slots; open-water evaporation on wet cells; soil ET (Feddes, or with Feddes off the
%   rest of the demand from the saturated part of the root zone through the table);
%   drain cascade; Darcy flux at the table. The caller applies recharge - capillary
%   with hp2d_water_table_update, then river exchange and lateral flow.
%   Outputs added to S: infiltration_mm, surface_evaporation_mm, soil_et_mm,
%   recharge_mm, capillary_mm.
valid = logical(P.valid(:));
N = numel(valid);
zeros_n = zeros(N, 1, 'like', S.head_m);
water = S.water_mm;
head = S.head_m(:);
surface = S.surface_depth_mm(:);
pervious = P.pervious(:) + zeros_n;
pervious(~valid) = 0;
root_depth = P.root_depth_m(:) + zeros_n;

S.infiltration_mm = zeros_n;
if P.infiltration
    V = hp2d_soil_view(column, head, water, root_depth);
    rate = hp2d_layered_infiltration_potential(V, surface, column.theta_r(:, 1), ...
        column.theta_s(:, 1), column.alpha_m_inverse(:, 1), column.n_vg(:, 1), ...
        P.ksat_near_surface_mm_h(:) + zeros_n, dt_s, P.l_vg, P.ltop_m, P.dh_max_m, ...
        P.ksurf_weight, valid);
    rate(~isfinite(rate)) = 0;
    request = min(rate*dt_s/3600.*pervious, surface);
    [water3, accepted] = hp2d_add_slot_infiltration(V, request, valid);
    water(valid, 1:3) = water3(valid, :);
    accepted(~valid) = 0;
    S.infiltration_mm = accepted;
    surface = max(surface - accepted, 0);
end

S.surface_evaporation_mm = zeros_n;
S.soil_et_mm = zeros_n;
if P.evapotranspiration
    wet = valid & surface > 0;
    evaporation = min(surface, S.open_water_evaporation_mm(:) + zeros_n);
    evaporation(~wet) = 0;
    S.surface_evaporation_mm = evaporation;
    surface = max(surface - evaporation, 0);
    V = hp2d_soil_view(column, head, water, root_depth);
    demand_step = S.soil_et_demand_mm(:) + zeros_n;
    demand = demand_step.*pervious;
    [water3, actual] = hp2d_slot_soil_et(V, demand, valid & ~wet, valid, P.feddes, ...
        column.theta_r(:, 2), column.theta_s(:, 2), column.alpha_m_inverse(:, 2), ...
        column.n_vg(:, 2), demand_step*(86400/dt_s));
    water(valid, 1:3) = water3(valid, :);
    actual(~valid) = 0;
    S.soil_et_mm = actual;
    if isempty(P.feddes) && P.groundwater
        % Roots below the water table take the rest of the demand from the aquifer.
        depth = hp2d_water_table_depth(column, head);
        rooted = max(root_depth, 0);
        saturated_root = min(max(rooted - depth, 0), rooted);
        fraction = zeros_n;
        has = rooted > 0;
        fraction(has) = saturated_root(has)./max(rooted(has), 1e-12);
        remaining = max(demand - S.soil_et_mm, 0);
        remaining(~(valid & ~wet)) = 0;
        take = min(remaining.*fraction, hp2d_drainable_storage_mm(column, head));
        [head, water, ~, unmet] = hp2d_water_table_update(column, head, water, -take);
        S.soil_et_mm = S.soil_et_mm + take - unmet;
    end
end

S.recharge_mm = zeros_n;
S.capillary_mm = zeros_n;
if P.groundwater
    water = hp2d_drain_cascade(column, head, water, dt_s, P.drainage_substeps);
    [water, recharge, capillary] = hp2d_water_table_flux(column, head, water, dt_s);
    recharge(~valid) = 0;
    capillary(~valid) = 0;
    S.recharge_mm = recharge;
    S.capillary_mm = capillary;
end
S.water_mm = water;
S.head_m = head;
S.surface_depth_mm = surface;
end
