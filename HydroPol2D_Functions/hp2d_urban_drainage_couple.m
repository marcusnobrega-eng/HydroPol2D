function [d_t, outlet_flow, Drainage_States] = hp2d_urban_drainage_couple(d_t, outlet_flow, Drainage_States, dt_s, cell_area_m2)
%HP2D_URBAN_DRAINAGE_COUPLE Exchange water between the surface and the
%   storm-drain network over one surface step, after the routing solver, as
%   HydroPol2D-Python does after local_inertial_step.
%
%   d_t: surface depth [mm]; outlet_flow: outlet flow [mm/h over the cell].
%   Only the inlet cells change depth. Water leaving through an outfall is
%   added to the outlet flow at the outfall's booking cell, so the outlet
%   record and the mass ledger count it (as hp2d_reservoir_book_external
%   does for reservoir releases). Water entering through a fixed-stage
%   outfall is kept in step_inflow_m3 for the ledger. Water in the nodes and
%   pipes is State.node_volume_m3.

N = Drainage_States.Network;
dt_s = double(gather(dt_s));
cell_area_m2 = double(gather(cell_area_m2));
cells = N.inlet_unique_cell;
surface_m3 = double(gather(d_t(cells))) / 1000 * cell_area_m2;
[Drainage_States.State, surface_m3, step] = hp2d_urban_drainage_advance( ...
    N, Drainage_States.State, surface_m3, dt_s);
d_t(cells) = surface_m3 / cell_area_m2 * 1000;

booked = step.node_outflow_m3 > 0;
if any(booked)
    cell = Drainage_States.booking_cell(booked);
    rate_mm_h = 1000 * 3600 * step.node_outflow_m3(booked) / (cell_area_m2 * dt_s);
    added = accumarray(cell, rate_mm_h, [numel(outlet_flow) 1]);
    target = find(added);
    outlet_flow(target) = outlet_flow(target) + added(target);
end
Drainage_States.step_outflow_m3 = step.outflow_m3;
Drainage_States.step_inflow_m3 = step.inflow_m3;
Drainage_States.last_substeps = step.substeps;
end
