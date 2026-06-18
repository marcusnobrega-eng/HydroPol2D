# P1-INFIL-LAYERS-001: Layered Soil Profile Fallbacks

Purpose: verify preprocessing logic for the proposed near-surface, root-zone, transmission-zone, and groundwater-zone discretization.

Expected behavior:
- root depths shallower than the near-surface target collapse into the near-surface layer;
- root depths deeper than bedrock are truncated by available soil depth;
- shallow groundwater truncates vadose layers at the water table;
- groundwater at the land surface produces zero vadose thickness;
- layer capacities are non-negative and sum to total vadose capacity.

Run:

```matlab
run('HydroPol2D_Model/Validation/Infiltration/LayeredSoil_Profile/run_layered_soil_profile_test.m')
```

Required diagnostics:
- `Outputs/Validation/Layered_Profile_Diagnostics.csv`
- `Outputs/Validation/Pass_Fail.csv`

## P1-INFIL-LAYERS-002: Layered Soil Dynamics

Purpose: verify active water-balance operations on the near-surface,
root-zone, and transmission-zone storages.

Expected behavior:
- infiltration fills layers top-down;
- evapotranspiration extracts only from root-accessible storage;
- roots shallower than the near-surface layer extract from the accessible
  near-surface fraction;
- layer percolation and groundwater recharge conserve storage;
- capacity reductions create saturation-excess water.

Run:

```matlab
run('HydroPol2D_Model/Validation/Infiltration/LayeredSoil_Profile/run_layered_soil_dynamics_test.m')
```

Required diagnostics:
- `Outputs/Validation/Layered_Dynamics_Diagnostics.csv`
- `Outputs/Validation/Layered_Dynamics_Pass_Fail.csv`
