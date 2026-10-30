# Rain gauges on meshes (VAL-MESH-GAUGE-001)

Purpose: a MATLAB mesh run can take rain-gauge series, as HydroPol2D-Python
does (audit M17).

Python spreads rain-gauge series over a mesh with normalised
inverse-distance weights (`forcing/gauge_rainfall.py`, `voronoi.py`):

- the weights are taken over the gauges that report;
- a point on a reporting gauge takes that gauge's value;
- the rain stops one interval after the last row.

MATLAB meshes took no gauges. A prepared mesh case can now give
`forcing.rain_gauges`, with fields `x`, `y`, `time_s` (row times, each the
start of its interval), `intensity_mm_h` (rows × gauges, NaN = not
reporting) and an optional `power` (default 2).

`hp2d_gauge_rainfall` sets up the weights once. With an overlap they go on
the overlap raster cells, which are then averaged over each mesh cell as
Python samples its HRU pieces. Without an overlap they go on the mesh cell
centres. Before the first row the rate is 0, as for every start-stamped
series in MATLAB. Python holds its first row there, which makes a
difference only if the series starts after the run does. Set
`config.forcing_interval_s` to the gauge interval so steps end at the row
times.

The check has two parts:

- A. IDW weights: three gauges, four points (one on a gauge), a gauge that
  stops reporting, and the rate outside the rows.
- B. A 2-hour run on a 2 x 2 mesh with walls and two gauges. The water
  stored at the end must equal the gauge rain.

## Run

```matlab
run('Validation/Voronoi/MeshRainGauges/run_mesh_rain_gauges_check.m')
```

Acceptance: IDW error below 1e-12 mm/h, rate 0 outside the rows, stored
water equal to the gauge rain (1e-12).
