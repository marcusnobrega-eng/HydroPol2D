# Staged overbank benchmark

This case tests the transition from channel flow to internal-bench flooding and
then to external-floodplain inundation. It compares three local-inertial models
of the same 1 km reach:

- a 5 m explicit representation of the full cross section;
- an ordinary 100 m grid sampled bilinearly from the 5 m rasters;
- a 100 m Neal subgrid model with a 50 m embedded channel.

The channel bed is at 0 m. Two 25 m benches begin at 0.5 m, and the external
floodplain begins at 1.0 m. Channel Manning roughness is 0.035; bench and
external-floodplain roughness is 0.10. The longitudinal slope is 0.005 m/m.

The inflow is a Nash-shaped hydrograph with a peak of 200 m3/s at 45 min. The
simulation lasts 240 min so the rising limb, both flooding transitions, and the
recession are represented. All models use the original Bates local-inertial
scheme and the same 0.01 min time step.

Run in MATLAB:

```matlab
cd('/Users/mngomes/Documents/HydroPol2D_CodexUpdate/HydroPol2D_Model/Validation/Subgrid/Neal2012_StagedOverbank_100m')
run_neal_staged_overbank_100m
```

Numerical results are written to `Outputs/Validation`, and figures are written
to `Figures`.

## Results

The analytical uniform-flow thresholds are approximately 31.4 m3/s at the
0.5 m bench elevation and 110.8 m3/s at the 1.0 m external-floodplain
elevation. The 200 m3/s inflow therefore crosses both transitions.

| Metric | Ordinary 100 m | Neal 100 m |
|---|---:|---:|
| Outlet hydrograph RMSE (m3/s) | 15.476 | 2.723 |
| Outlet hydrograph NSE | 0.945 | 0.998 |
| Internal hydrograph RMSE (m3/s) | 12.053 | 3.043 |
| Internal hydrograph NSE | 0.967 | 0.998 |
| Bench-depth RMSE (m) | 0.180 | 0.018 |
| Maximum-depth RMSE (m) | 0.280 | 0.040 |
| Wet-area CSI at 0.01 m | 0.334 | 0.998 |

All three simulations conserve the imposed water volume to numerical
precision. The ordinary 100 m grid does not reach the external floodplain at
the internal gauge because resampling removes the two 25 m benches and turns
the central coarse row into an overly wide channel.

The Neal model reproduces the rising limb, peak stage, overbank activation,
outlet hydrograph, and maximum inundation pattern much more closely. The first
wetting of the internal benches occurs at 21 min in both the 5 m and Neal
models. External-floodplain wetting occurs at 34 min in the 5 m model and at
35 min in the Neal model.

The recession remains unresolved at 100 m. At the internal section, the Neal
benches fall below 0.01 m at 126 min, compared with 196 min in the 5 m model.
The external floodplain falls below the same threshold at 152 and 237 min,
respectively. The case therefore passes mass conservation, hydrograph, peak
extent, and activation-timing checks, but fails the recession-timing check.
It is reported as a partial validation rather than a full pass.
