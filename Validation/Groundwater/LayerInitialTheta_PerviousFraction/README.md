# Per-layer initial moisture with the pervious fraction (VAL-SOIL-THETA0-001)

Purpose: the layered soil starts with the water the inputs describe.

HydroPol2D-Python starts each soil layer with f_p (θ − θr) × thickness.
θ is that layer's initial-moisture raster or, when the raster is absent,
the cell's initial θ (`modules/soil.py`, `core/preprocess.py`). MATLAB
(audit M8):

- used the three layer rasters (`InitialNearSurfaceTheta`,
  `InitialRootZoneTheta`, `InitialTransmissionTheta`) only when all three
  existed;
- then wrote (θ − θr) × thickness × 1000 with no f_p, so a half-paved cell
  started with twice its water, often full.

Now each raster is optional. `HydroPol2D_preprocessing.m` stores the layers
that have one in `Soil_Properties.initial_layer_theta`.
`initialize_layered_soil_storage.m` uses each layer's θ, or θi, in the f_p
form it already used for θi, bounded by the layer capacity.

The check builds the storages of six cells, with pervious fractions 1 to 0,
when only the root-zone raster is given. It compares them with Python's
`initialize_layered_soil_storage` on the same numbers: MATLAB writes
`layer_inputs.csv`, and `make_python_reference.py` turns it into
`python_reference.csv` (HydroPol2D-Python 66cd12f).

## Run

```matlab
run('Validation/Groundwater/LayerInitialTheta_PerviousFraction/run_layer_initial_theta_check.m')
```

Acceptance: storages within 1e-9 mm of Python. Result: 0 mm. On the old code
the root-zone raster was ignored when it came alone, 45 mm off.
