# SNISB Dam-Case Runners

These MATLAB functions run HydroPol2D for dam-case folders prepared from the SNISB workflow. They are intended for batch execution after the terrain, land-cover, soil, inflow, and outlet inputs have been prepared for each dam.

## Expected Case Inputs

Each dam folder is expected to contain:

```text
dam_dir/
  rasters/domain_clipped/
    dem_fabdem_30m_domain.tif
    lulc_mapbiomas_30m_domain.tif
    soil_texture_usda_30m_domain.tif
  hydrograph/inlet_cells_hydrograph.csv
  outlet/outlet_cells.csv
```

The outlet file is preferred. When it is absent, the runner uses HydroPol2D's fallback outlet placement.

## Runners

Add this folder to the MATLAB path, then choose the runner that matches the execution pattern:

```matlab
addpath(fullfile(pwd, 'Applications', 'SNISB'))

% Run one prepared dam folder.
run_snisb_hydropol2d_case('/path/to/dam_dir')

% Run a dam selected from a CSV list, for example in a SLURM array task.
run_snisb_hydropol2d_from_list('/path/to/dam_list.csv', 1)

% Run a batch of prepared dam folders.
run_snisb_hydropol2d_batch
```

`run_snisb_hydropol2d_case` accepts name-value options for routing, duration, output interval, rainfall, infiltration, warmup depth, and boundary treatment. Its configuration bridge is [`input_data_bypass_snisb.m`](input_data_bypass_snisb.m).

Outputs are written within each `dam_dir/hydropol2d/` folder so that independent cases can be run and retained separately.
