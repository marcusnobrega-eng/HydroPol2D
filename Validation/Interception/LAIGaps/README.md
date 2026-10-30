# Active cells without an LAI value keep their rain (VAL-INT-LAIGAP-001)

Purpose: a gap in the LAI raster must not swallow rain (audit L10).

A NaN LAI in an active cell made `interceptionModel`'s canopy store NaN.
The cell's throughfall was then 0, and its rain never reached the ground.
HydroPol2D-Python fills such cells with the configured LAI fallback
(`static.lai`, default 0) and warns (`core/preprocess.py`).
`HydroPol2D_preprocessing.m` now does the same after its LAI smoothing, with
`hp2d_fill_lai_gaps`. The fallback is `lai_fallback` (Flags sheet or bypass
flags), default 0.

The case is a 2 x 2 LAI map with one gap in an active cell and one cell
outside the domain, under 5 mm of rain on an empty canopy.

## Run

```matlab
run('Validation/Interception/LAIGaps/run_lai_gap_check.m')
```

Acceptance: one cell filled. The gap cell gets all 5 mm as throughfall
(0 mm unfilled), and the outside cell stays NaN.
