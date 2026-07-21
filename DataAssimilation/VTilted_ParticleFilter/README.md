# V-Tilted Augmented Particle Filter Prototype

This folder contains an isolated data-assimilation prototype. It does not edit
the production HydroPol2D source or the shared V-tilted static rasters.

## Run The Full HydroPol2D Backend

From MATLAB:

```matlab
cd('/Users/mngomes/Documents/HydroPol2D_CodexUpdate/HydroPol2D_Model/DataAssimilation/VTilted_ParticleFilter')
test_pf_core_math
run_vtilted_augmented_pf_full_hydropol2d('SmokeTest', true, 'NParticles', 1, 'NWindows', 1, 'Routing', 'full_momentum')
run_vtilted_augmented_pf_full_hydropol2d('SmokeTest', true, 'NParticles', 1, 'NWindows', 1, 'Routing', 'local_inertial')
```

`run_vtilted_augmented_pf_full_hydropol2d.m` is the authoritative prototype
for data assimilation because it calls HydroPol2D preprocessing and the real
HydroPol2D time loop for the synthetic truth and every particle.

`run_vtilted_augmented_pf.m` remains only as a fast mathematical sandbox for
particle-filter plotting and objective-function experiments. It should not be
used to claim HydroPol2D data-assimilation behavior.

## What The Prototype Does

- Uses augmented particles with explicit `theta + state` memory.
- Starts with two LULC classes and two SOIL classes.
- Includes groundwater parameters from the beginning.
- Solves the full HydroPol2D model in the full backend: hydrologic processes,
  groundwater/recharge, and the selected routing solver.
- Runs the default synthetic case for 360 min with constant rainfall.
- Uses a 30 min assimilation window with mixed observation clocks:
  discharge every 5 min, soil moisture every 15 min, and groundwater depth
  every 30 min.
- Treats discharge observations as interval means over their sample interval;
  soil moisture and groundwater depth are sampled as point states.
- Updates parameters only in v1.
- Carries and resamples state explicitly.
- Computes particle weights with the paper method:
  `sum(beta * normal_pdf(error, sigma))`, then normalizes by the sum over
  particles.
- Converts discharge from `m3/s` to area-normalized `mm/h` before weighting.
- Uses default paper-style group weights `beta_Q = beta_GW = beta_SM = 1/3`.
- Carries particle weights between windows, resamples only when effective
  sample size falls below the configured threshold, preserves one unperturbed
  copy of each selected particle, and perturbs only extra duplicated copies.
- Uses `pf_apply_state_update_identity.m` as the placeholder for future state
  correction.
- Keeps all DA window-control edits in the copied loop
  `HydroPol2D_Main_While_DA_Window.m`, not in the production model loop.

## Outputs

The compact sandbox writes:

- `Outputs/PF_Results.mat`
- `Outputs/Parameter_Evolution.csv`
- `Outputs/Objective_Function.csv`
- `Outputs/Effective_Sample_Size.csv`
- `Outputs/Observation_Fit.csv`
- `Outputs/Observation_Ensemble_Summary.csv`
- `Outputs/Posterior_Summary.csv`
- `Outputs/Posterior_Summary_Physical.csv`
- `Outputs/State_Summary.csv`
- `Figures/parameter_evolution.png`
- `Figures/initial_parameter_spread.png`
- `Figures/objective_function.png`
- `Figures/observation_fit.png`
- `Figures/observation_ensemble_fit.png`
- `Figures/particle_weights.png`
- `Figures/state_summary.png`

The full HydroPol2D backend writes:

- `Outputs/FullHydroPol2D_PF_Results.mat`
- `Outputs/FullHydroPol2D_Synthetic_Observations.csv`
- `Outputs/FullHydroPol2D_Observation_Fit.csv`
- `Outputs/FullHydroPol2D_Parameter_Evolution.csv`
- `Outputs/FullHydroPol2D_Effective_Sample_Size.csv`

## Isolation Contract

All generated files stay under this folder. The normal HydroPol2D workflow is
unchanged.
