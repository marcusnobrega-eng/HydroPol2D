# Reservoir releases inside the routing solvers (VAL-RES-002)

Purpose: check the reservoir controls where the model applies them, inside
each surface solver. VAL-RES-001 integrates the rating curve in its own
script and never calls a solver.

All solvers now call `hp2d_reservoir_releases.m`, the rule of
`reservoir_outlet_transfers` in HydroPol2D-Python:

- Each control drains its source cell through up to two outlets, in order:
  Q = k max(h - h0, 0)^b, capped at the water above the sill.
- h is the source depth after this step's face fluxes and after outlet 1.
- The release leaves the source and enters the destination cell.
- A destination of (0, 0) in the workbook sends the release out of the
  model. It is booked as outlet flow at the source cell, so the mass ledger
  counts it. NaN switches the outlet off.
- Neal cells measure the sill above the channel bed.

The case is a 5 x 11 grid of 10 m cells. Six pits sit among cells outside the
domain, so no water crosses a face. One 1 s step per solver: local inertial,
Neal, full momentum, diffusive, kinematic and CA.

- Control 1: two outlets to two pits.
- Control 2: one outlet out of the model.
- Control 3: one outlet to a pit.

The check computes the reference volumes itself.

Before the fix:

- The default, Neal, full-momentum and wave solvers added the release to the
  source, which created water.
- Outlet 1 reached its destination twice when outlet 2 existed.
- CA lost outlet 1 when a control had no outlet 2.
- Neal delivered nothing to the destinations.
- A (0, 0) destination was an invalid index.

The lookup-subgrid solver (`Local_Inertial_Model_D4_Subgrid.m`) uses the same
helper with its SFINCS tables. It is not in this check, because building those
tables needs a fine DEM.

## Run

```matlab
run('Validation/Reservoir/SolverExchange_MassConservation/run_reservoir_solver_exchange_check.m')
```

Acceptance: for every solver, each cell volume, the external volume and the
ledger residual within 1e-9 m3 of the reference.
