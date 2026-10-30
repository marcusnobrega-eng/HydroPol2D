function E = hp2d_detention_pond_exchange(surface_m3, storage_m3, pump_active, cumulative_m3, P, dt_s)
%HP2D_DETENTION_POND_EXCHANGE Critical capture into the ponds, then all their devices.
%   surface_m3: surface-water volume of the exchange cells P.exchange_cells [m3];
%   storage_m3: pond storage [m3]; pump_active: pump switch of every device;
%   cumulative_m3: cumulative device volume [m3]; P: hp2d_prepare_detention_ponds.
%   pond_exchange_step in HydroPol2D-Python (routing/detention.py):
%     1. Critical capture. Each face between a pond and an outside cell takes
%        q = w sqrt(g) H^1.5, H = max(z + h - max(pond stage, sill), 0), from that
%        cell into the pond; a donor shared by several faces is scaled to its volume.
%     2. Devices, all at once from the stages after the capture, with
%        dH = max(stage - level, 0) - max(destination stage - level, 0):
%          orifice   Q = Cd A sqrt(2 g |dH|) sign(dH)   (reverses when dH < 0)
%          spillway  Q = Cd L |dH|^1.5 sign(dH)
%          pump      on at stage >= start, off at stage <= stop, else as before;
%                    Q = Qmax sqrt(max(1 - H_lift / H_shutoff, 0)),
%                    H_lift = max(destination - pond stage, 0) + extra loss
%        Releases from a pond are scaled to its storage and reverse flows to the
%        destination cell's volume, so the exchange conserves mass. A destination of
%        0 (workbook (0, 0)) leaves the model: a free outfall that never reverses;
%        E.external_m3 is that volume, per device, for the outlet-flow booking.

if ~(dt_s > 0)
    error('HydroPol2D:DetentionPonds', 'Pond exchange dt_s must be positive.');
end
g = 9.81;
n = P.n_ponds;
A = P.cell_area_m2;
m = numel(P.exchange_cells);
surface = max(double(surface_m3(:)), 0);
storage = max(double(storage_m3(:)), 0);
pump_active = logical(pump_active(:));
cumulative = double(cumulative_m3(:));
terrain = P.exchange_terrain_m;

% 1. Critical capture.
stage = hp2d_detention_pond_stage(storage, P);
fl = P.face_local;
fp = P.face_pond;
neighbour_stage = terrain(fl) + surface(fl) / A;
head = max(neighbour_stage - max(stage(fp), P.face_sill_m), 0);
request = P.face_width_m .* sqrt(g) .* head .^ 1.5 .* dt_s;
requested_by_cell = accumarray(fl, request, [m 1]);
donor_scale = ones(m, 1);
has = requested_by_cell > 0;
donor_scale(has) = min(1, surface(has) ./ max(requested_by_cell(has), 1e-30));
captured = request .* donor_scale(fl);
surface = surface - accumarray(fl, captured, [m 1]);
lateral = accumarray(fp, captured, [n 1]);
storage = storage + lateral;

% 2. Devices.
stage = hp2d_detention_pond_stage(storage, P);
dp = P.dev_pond;
dl = P.dev_local;
inside = dl > 0;
pond_stage = stage(dp);
destination_stage = -inf(P.n_devices, 1);
destination_stage(inside) = terrain(dl(inside)) + surface(dl(inside)) / A;
level = P.dev_level_m;
delta_head = max(pond_stage - level, 0) - max(destination_stage - level, 0);
orifice_q = P.dev_coefficient .* P.dev_geometry .* sqrt(2 * g * abs(delta_head)) .* sign(delta_head);
spillway_q = P.dev_coefficient .* P.dev_geometry .* abs(delta_head) .^ 1.5 .* sign(delta_head);
pump_now = pump_active;
pump_now(pond_stage <= P.dev_stop_stage_m) = false;
pump_now(pond_stage >= P.dev_start_stage_m) = true;
lift = max(destination_stage - pond_stage, 0) + P.dev_head_loss_m;
pump_q = zeros(P.n_devices, 1);
pump_q(pump_now) = P.dev_maximum_flow_m3_s(pump_now) .* ...
    sqrt(max(1 - lift(pump_now) ./ P.dev_shutoff_head_m(pump_now), 0));
requested_q = pump_q;
requested_q(P.dev_kind == 0) = orifice_q(P.dev_kind == 0);
requested_q(P.dev_kind == 1) = spillway_q(P.dev_kind == 1);
requested_q(~P.dev_enabled) = 0;
is_pump = P.dev_kind == 2;
active_after = pump_active;
active_after(is_pump) = pump_now(is_pump);

requested = requested_q * dt_s;
pond_outgoing = accumarray(dp, max(requested, 0), [n 1]);
cell_outgoing = accumarray(dl(inside), max(-requested(inside), 0), [m 1]);
pond_scale = ones(n, 1);
has = pond_outgoing > 0;
pond_scale(has) = min(1, storage(has) ./ max(pond_outgoing(has), 1e-30));
cell_scale = ones(m, 1);
has = cell_outgoing > 0;
cell_scale(has) = min(1, surface(has) ./ max(cell_outgoing(has), 1e-30));
actual = requested .* pond_scale(dp);
reverse = requested < 0;
actual(reverse) = requested(reverse) .* cell_scale(dl(reverse));
storage = storage - accumarray(dp, actual, [n 1]);
surface = surface + accumarray(dl(inside), actual(inside), [m 1]);
external = zeros(P.n_devices, 1);
external(~inside) = actual(~inside);

E = struct('surface_m3', max(surface, 0), 'storage_m3', max(storage, 0), ...
    'pump_active', active_after, 'cumulative_m3', cumulative + actual, ...
    'lateral_inflow_m3', lateral, 'device_discharge_m3_s', actual / dt_s, ...
    'external_m3', external);
end
