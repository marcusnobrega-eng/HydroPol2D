# V-Tilted validation example

This portable example runs the canonical V-Tilted case through either Excel
workbooks or MATLAB configuration code. Both entry points use the same rasters,
forcing, outlet, solver settings, and prepared Voronoi mesh:

- rainfall intensity: **10.8 mm/h**;
- rainfall duration: **90 minutes**;
- simulation duration: **180 minutes**, so the hydrograph propagation and
  recession are recorded after the storm ends;
- Manning roughness: **0.04**;
- maximum time step: **1 s**;
- outlet: **60 m**, represented by three 20 m cells or the matching Voronoi
  boundary edges;
- outlet slope: **0.004 m/m**.

`Input_Data_Sheets/General_Data.xlsx` remains the main simulation workbook.
The `flag_voronoi` field in its `Flags` worksheet selects the mesh:

| `flag_voronoi` | Mesh path |
|---|---|
| `0` | 20 m raster grid from `Static/DEM.tif` |
| `1` | Prepared UGRID mesh referenced by the `Voronoi case file` field |

The distributed workbook defaults to `flag_voronoi = 0`. Switching the flag
does not change the storm or create a second `General_Data.xlsx`.
`Input_Data_Sheets/Voronoi_Settings.xlsx` holds the prepared case path,
mesh-compatibility metadata, descriptions, and Voronoi output controls.

## Run with Excel inputs

Open MATLAB in the HydroPol2D repository and run:

```matlab
run('Examples/VTilted_Excel/run_excel.m')
```

The launcher reads `flag_voronoi` and writes results to
`Examples/VTilted_Excel/Outputs/regular` or
`Examples/VTilted_Excel/Outputs/voronoi`. The output directory is ignored by Git.

The distributed workbook defaults to the regular grid. To run the prepared
Voronoi mesh, find `flag_voronoi` in the `Flags` worksheet and change its value
from `0` to `1`. Change it back to `0` for the raster grid. Do not create a
second `General_Data.xlsx` file.

## Run with MATLAB configuration code

From the repository root:

```matlab
addpath('Examples/VTilted_Excel')
run_code('regular')
run_code('voronoi')
```

The code path reads
`Config/vtilted_input_data.m` and `Config/vtilted_input_paths.m`. It does not
read model parameters from Excel. Results are written to `Outputs/code_regular`
or `Outputs/code_voronoi`. The launcher restores the caller's environment
variables after each run.

## Rebuild the validation inputs

The repository already contains the prepared inputs. To regenerate the
canonical rasters and prepared Voronoi case, run:

```matlab
run('Examples/VTilted_Excel/prepare_standard_case.m')
```

This produces a 41 x 31 regular grid at 20 m resolution and a matching
uniform 20 m Voronoi/UGRID case with 1,271 cells.

## Verified results

All four combinations completed to 180 minutes with MATLAB R2025b on
September 16, 2026:

- **Regular 20 m grid:** Excel and MATLAB-code inputs produced identical 181
  hydrograph samples. The peak was 1.524886608 m3/s at 89.991 min, discharge
  at the end of the run was 0.120215490 m3/s, and outlet runoff was
  15.186905861 mm.
- **Uniform 20 m Voronoi:** Excel and MATLAB-code inputs produced identical
  numerical diagnostics. Recorded rainfall volume was 8,236.08 m3, maximum
  mass-residual fraction was 7.12e-15, maximum depth was 0.155864 m, and the
  run passed its internal checks.
- **Canonical D4/Voronoi validation:** the 180-minute uniform-20 m comparison
  passed, with NSE 0.99881, peak-discharge error 0.000686%, outlet-volume error
  1.1547%, and final-depth RMSE 0.001526 m.

The separate variable-resolution Voronoi validation remains useful for mesh
development, but it currently fails the mesh-to-mesh convergence gate for
outlet volume and final depth. It is therefore not the default Excel
validation fixture and must not be presented as proof of regular/Voronoi
equivalence.

## What the workbooks control

The main workbook contains the simulation period, time-step controls, input
paths, process flags, and boundary method. The companion files must remain in
the same `Input_Data_Sheets` directory:

- `Voronoi_Settings.xlsx`
- `LULC_parameters.xlsx`
- `SOIL_parameters.xlsx`
- `Rainfall_Intensity_Data.xlsx`

Raster paths may be absolute or relative to the main workbook. This example
uses `../Static/...`, so the whole example directory can be moved without
rewriting paths.

With `flag_voronoi = 0`, the workbook reads the DEM, LULC, and soil rasters and
runs the structured-grid solver. With `flag_voronoi = 1`, the `Voronoi case
file` cell in `Voronoi_Settings.xlsx` selects a prepared case containing the
validated UGRID mesh and raster-overlap reference. The compatibility settings
must match the metadata stored in that case. Excel selects the prepared mesh;
it does not regenerate the mesh during a simulation.

The rainfall workbook stores 10.8 mm/h at 90 minutes and 0 mm/h at 180
minutes. HydroPol2D treats those times as interval ends. Rain therefore falls
from 0 to 90 minutes, while routing continues without rain from 90 to 180
minutes.

## Run another Excel case without editing the launcher

Set these environment variables in MATLAB, then run the normal launcher:

```matlab
setenv('HYDROPOL_RUN_MODE','excel')
setenv('HYDROPOL_INPUT_EXCEL_FILE','/absolute/path/to/General_Data.xlsx')
setenv('HYDROPOL_EXPORT_ROOT_DIR','/absolute/path/to/Outputs/my_run')
setenv('HYDROPOL_SKIP_POSTPROCESS','')
run('HydroPol2D_V115.m')
```

Set `HYDROPOL_SKIP_POSTPROCESS` to `1` only for a solver-only diagnostic run.
The launcher automatically adds the selected workbook directory to the MATLAB
path so its companion spreadsheets are found.

Excel and code mode use the same numerical solver. Their verified regular-grid
hydrographs are identical, so users can select the input style that best fits
interactive work or version-controlled studies.
