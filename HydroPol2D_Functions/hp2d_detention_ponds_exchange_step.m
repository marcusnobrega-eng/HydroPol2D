function [Detention_Ponds, d_t, outlet_flow] = hp2d_detention_ponds_exchange_step(Detention_Ponds, d_t, outlet_flow, dt_s)
%HP2D_DETENTION_PONDS_EXCHANGE_STEP Critical capture and devices after the routing step.
%   d_t: surface depth [mm] after this step's face fluxes; outlet_flow: outlet sink
%   [mm/h over the coarse cell]. Runs hp2d_detention_pond_exchange on the face donors
%   and device destinations and writes their new depths back. A device discharging
%   out of the model is booked as outlet flow at its inlet cell (or the pond's first
%   cell), as hp2d_reservoir_book_external does, so the outlet record and the mass
%   ledgers count it.

P = Detention_Ponds;
A = P.cell_area_m2;
cells = P.exchange_cells;
depth_before = double(gather(d_t(cells)));
depth_before(~isfinite(depth_before)) = 0;
volume_before = max(depth_before, 0) / 1000 * A;

E = hp2d_detention_pond_exchange(volume_before, P.storage_m3, P.pump_active, ...
    P.cumulative_device_m3, P, dt_s);

changed = E.surface_m3 ~= volume_before;
if any(changed)
    d_t(cells(changed)) = E.surface_m3(changed) / A * 1000;
end
P.storage_m3 = E.storage_m3;
P.pump_active = E.pump_active;
P.cumulative_device_m3 = E.cumulative_m3;
P.device_discharge_m3_s = E.device_discharge_m3_s;
P.cumulative_lateral_m3 = P.cumulative_lateral_m3 + E.lateral_inflow_m3;
leaving = E.external_m3 > 0;
if any(leaving)
    P.cumulative_external_m3 = P.cumulative_external_m3 + ...
        accumarray(P.dev_pond(leaving), E.external_m3(leaving), [P.n_ponds 1]);
    [book, ~, which] = unique(P.dev_book_cell(leaving));
    booked = accumarray(which, E.external_m3(leaving));
    outlet_flow(book) = outlet_flow(book) + 1000 * 3600 * booked / (A * dt_s);
end
Detention_Ponds = P;
end
