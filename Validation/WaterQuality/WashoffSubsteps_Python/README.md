# Wash-off sub-steps and limiter (VAL-WQ-SUBSTEP-001)

Purpose: `build_up_wash_off.m` integrates the whole step, never washes off more
pollutant than a cell holds, and does not slow the model down.

Three MATLAB problems, from the audit:

- **M15.** A step the wash-off could not take in one go was split into
  floor(T/tmin)+1 pieces, and the last piece lasted
  (T/tmin − floor(T/tmin))·tmin. With T = 60 s and tmin = 25 s that
  integrates 2 × 20 s + 10 s = 50 s, not 60 s.
- **L20.** There was no limiter, so a cell could release more than it held.
  The negative mass was then booked as "lost".
- **L21.** The model time step was capped at the wash-off's stable step.

HydroPol2D-Python (`modules/water_quality.py`) splits the step into
n = ceil(max(W_tot/B)·dt) equal sub-steps. In each sub-step it scales a
cell's outflow to the mass that cell holds. The pollutant sub-cycles inside
the hydraulic step. MATLAB now does the same and keeps its own units and
options:

- `C_3` per hour, with W in kg/h.
- The rating-curve wash-off.
- The outlet direction.

`refreshing_timestep.m` no longer limits the step for water quality.

The check runs one 120 s step on a 3 x 4 grid, with flows in all four
directions and out of every side. C_3 is high enough to need many
sub-steps, and one cell is below `min_Bt`. The masses and the total export
are compared with Python's `mass_based_washoff_step` on the same numbers.
MATLAB writes the inputs to `washoff_inputs.csv`, and `make_python_reference.py`
turns them into `python_reference.csv` (HydroPol2D-Python b5f224b). Python
reads C_3 per second, so it gets C_3/3600.

## Run

```matlab
run('Validation/WaterQuality/WashoffSubsteps_Python/run_washoff_substep_check.m')
```

Acceptance:

- Masses and export within 1e-12 kg of Python.
- No mass booked as lost, and no negative mass.
- `tmin_wq` equals the step.

| | Old code | New code |
|---|---|---|
| Mass error | 0.036 kg | 6.6e-17 kg |
| Export error | 0.13 kg | 1.8e-15 kg |
| Mass booked as lost | 0.034 kg | none |

The analytical wash-off check (VAL-WQ-001A–C, BuildUpWashOff_MassBalance) still
passes.
