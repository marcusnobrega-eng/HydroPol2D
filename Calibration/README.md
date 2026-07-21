# HydroPol2D Calibration Framework

This folder contains calibration wrappers around HydroPol2D. The normal
model source remains unchanged: calibration candidates are passed through
existing bypass/override inputs.

## Generic GA

The generic single-objective genetic algorithm is implemented in:

```matlab
Calibration/GenericGA/run_hydropol2d_ga_calibration.m
```

It expects a configuration struct with:

- `parameter_specs`: parameter names, targets, bounds, transforms, and enabled flags.
- `cases`: rainfall/event cases to evaluate.
- `model_runner`: a case-specific function handle that runs HydroPol2D.
- `objective`: weights for RMSE, volume error, peak error, timing error, and bias.

Generic helper files:

- `hydropol2d_default_parameter_catalog.m`: default LULC, SOIL, groundwater, and generic input-data parameter ranges.
- `hydropol2d_build_parameter_specs.m`: expands the catalog to however many LULC/SOIL classes your catchment has.
- `hydropol2d_enable_parameter_specs.m`: enables exact parameter names or wildcard patterns.
- `hydropol2d_finalize_parameter_specs.m`: validates bounds and creates internal GA coordinates.
- `hydropol2d_calibration_config_template.m`: commented template for a new catchment.

The key idea is that class count is not hard coded. A catchment config only
needs to provide:

```matlab
Cal.base_overrides.LULC.Index = [...];
Cal.base_overrides.LULC.roughness = [...];
Cal.base_overrides.LULC.root_depth_m = [...];

Cal.base_overrides.SOIL.Index = [...];
Cal.base_overrides.SOIL.theta_i = [...];
Cal.base_overrides.SOIL.Ks_multiplier_near_surface = [...];
...
```

If there are 3 LULC classes and 5 soil classes, specs are built for those.
If there are 20 LULC classes and 30 soil classes, specs are built for those
instead.

Enable parameters using exact names:

```matlab
enableNames = [
    "lulc_roughness_multiplier_all"
    "soil_theta_i_multiplier_all"];
Cal.parameter_specs = hydropol2d_enable_parameter_specs(Cal.parameter_specs, enableNames);
```

Or wildcard patterns:

```matlab
enableNames = [
    "lulc_*_roughness"
    "soil_*_theta_i"];
```

Use global multipliers early in calibration. Use class-specific parameters
only when observations can identify them.

## Runner Contract

The generic engine does not choose a forcing dataset or run HydroPol2D on
its own. The catchment-specific `Cal.model_runner` must run one candidate
for one case and return either a finite scalar `Eval.objective`, or an
`Eval.metrics` structure containing:

```matlab
RMSE_m3_s
MAE_m3_s
Bias_m3_s
NSE
Observed_Volume_m3
Modeled_Volume_m3
Volume_Error_pct
Observed_Peak_m3_s
Modeled_Peak_m3_s
Observed_Peak_Time_min
Modeled_Peak_Time_min
Peak_Timing_Error_min
```

When `Eval.objective` is absent, the engine minimizes the weighted,
dimensionless objective

```text
J = w_rmse RMSE/Qobs,peak
  + w_volume |Vsim - Vobs|/Vobs
  + w_peak |Qsim,peak - Qobs,peak|/Qobs,peak
  + w_timing |dtpeak|/Tscale
  + w_bias |Bias|/Qobs,peak.
```

The weights and `Tscale` are defined in `Cal.objective`. A runner may
instead provide a finite `Eval.objective` when a different, documented
objective is required.

## Outputs

Each calibration writes its population history and restart state to
`Cal.output_root`:

- `GA_History.csv` and `Generation_###.csv` record every evaluation.
- `Best_Parameters.csv` and `BestSoFar.mat` record the current optimum.
- `GA_Checkpoint.mat` and `GA_Results.mat` preserve restartable results.
- `GA_Diagnostics.png` plots objective evolution by generation.

## Repeatable Engine Test

Run the following command from the repository root before using a new
catchment runner:

```matlab
addpath('Calibration/GenericGA')
test_generic_ga_known_optimum
```

This deterministic two-parameter inverse problem starts away from the
known optimum, exercises the metric-based objective, and verifies that the
GA recovers the target parameter values.

## Adding A New Catchment

1. Copy `GenericGA/hydropol2d_calibration_config_template.m` into a new
   catchment calibration folder.
2. Fill in `case_root`, `Cal.base_overrides`, `Cal.cases`, and
   `Cal.model_runner`.
3. Keep the generic catalog if the default ranges are acceptable, or edit
   `Catalog.LULC.Fields`, `Catalog.SOIL.Fields`, and `Catalog.InputData`
   before building specs.
4. Enable a small identifiable subset first.
5. Run a smoke test with `1` generation and `2` individuals before launching
   a full calibration.
