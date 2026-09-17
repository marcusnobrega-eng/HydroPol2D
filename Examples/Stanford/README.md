# Stanford 30 m applied example

This portable example runs the Stanford watershed with either Excel workbooks
or MATLAB configuration code. Both entry points use the same 30 m rasters,
soil and land-cover tables, rainfall event, numerical settings, and 180-minute
simulation period.

## Canonical setup

- grid resolution: **30 m**;
- rainfall: **100 mm/h from 0 to 60 minutes**;
- simulation duration: **180 minutes**;
- map and hydrograph interval: **15 minutes**;
- local-inertial routing with adaptive time stepping;
- CPU and double precision by default;
- raster paths relative to `Input_Data_Sheets/General_Data.xlsx`.

The rainfall ends at 60 minutes, while routing continues to 180 minutes so
the hydrograph response and recession are retained.

## Run with Excel inputs

Open MATLAB in the HydroPol2D repository root and run:

```matlab
run('Examples/Stanford/run_excel.m')
```

The launcher reads:

- `Input_Data_Sheets/General_Data.xlsx`;
- `Input_Data_Sheets/LULC_parameters.xlsx`;
- `Input_Data_Sheets/SOIL_parameters.xlsx`;
- `Input_Data_Sheets/Rainfall_Intensity_Data.xlsx`;
- the rasters under `Static/`.

Results are written to:

```text
Examples/Stanford/Outputs/Reruns_CurrentModel/Stanford_30m_Excel
```

The output directory is ignored by Git.

## Run with MATLAB configuration code

From the repository root:

```matlab
addpath('Examples/Stanford')
run_code
```

This route uses `Config/input_data_bypass_script.m` and
`Config/input_paths_bypass.m`. Results are written to:

```text
Examples/Stanford/Outputs/Reruns_CurrentModel/Stanford_30m_Code
```

Advanced overrides remain available through the lower-level runner. For
example:

```matlab
run_stanford_local_case( ...
    'SimulationMinutes', 360, ...
    'RecordTimeMapsMinutes', 30, ...
    'OutputTag', 'Stanford_30m_6h')
```

## Folder layout

```text
Stanford/
|-- README.md
|-- run_excel.m
|-- run_code.m
|-- run_stanford_local_case.m
|-- Config/
|   |-- input_data_bypass_script.m
|   `-- input_paths_bypass.m
|-- Input_Data_Sheets/
|   |-- General_Data.xlsx
|   |-- LULC_parameters.xlsx
|   |-- SOIL_parameters.xlsx
|   `-- Rainfall_Intensity_Data.xlsx
|-- Static/
|   |-- DEM.tif
|   |-- LULC.tif
|   |-- SOIL.tif
|   `-- optional supporting rasters
|-- Forcing/Rainfall/
|   `-- Rainfall_Intensity_Event_Clean.csv
`-- Outputs/                         # generated locally; not versioned
```

## Editing the case

For Excel runs, edit the workbooks under this example's
`Input_Data_Sheets/` directory. Do not edit the shared templates at the
repository root for a case-specific experiment. Keep the relative raster
paths (`../Static/...`) when copying the complete `Stanford` directory.

For code-driven studies, edit the two files under `Config/` or pass supported
overrides to `run_stanford_local_case`.

`Static_10m_LiDAR`, `Static_10m_LiDAR_ModelDEM`, and
`Static_10m_Upsampled` are research variants. They are not selected by the
default 30 m launchers and should not be presented as the introductory case.
