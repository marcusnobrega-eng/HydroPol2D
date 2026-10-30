# Mesh soil layers: per-layer Ksat and the soil's l (VAL-MESH-SOIL-001)

Purpose: the legacy mesh soil drains like HydroPol2D-Python's (audit M19).

Python's legacy mesh soil (`voronoi.py`) drains the near-surface, root-zone
and transmission layers at their own saturated conductivities. Capillary rise
uses the root-zone conductivity, and Mualem's pore connectivity is the soil's
`l_vg`. MATLAB used one `ksat_m_s` for every layer and fixed l at 0.5 in the
drainage.

`Voronoi_Hydrology_Initialize` now takes the optional HRU fields
`ksat_root_m_s` and `ksat_transmission_m_s`, which default to `ksat_m_s`.
`Voronoi_Hydrology_Step` drains each layer at its own Ksat with `l_vg`. The
relative conductivity is clipped to [0, 1], as in Python. Capillary rise uses
`ksat_root_m_s`. Cases that give neither field and keep `l_vg` at its
default of 0.5 are unchanged.

The check uses one 100 m2 cell. The root zone holds water at saturation 0.6
and the other layers are empty. The transmission layer is nearly
impermeable, and there is no rain, ET or groundwater. After one 600 s step
the water that left the root zone must equal f_p · K_root · kr(Se; l) · dt,
with K_root = 2e-5 m/s and l = −1.2. The near-surface Ksat is 1e-7 m/s.

## Run

```matlab
run('Validation/Voronoi/MeshSoilLayers/run_mesh_soil_layers_check.m')
```

Acceptance: error below 1e-12 m. Result: 2.5e-18 m. On the old code the root
zone drained 9.7e-8 m (Ksat 1e-7, l = 0.5) instead of 4.6e-5 m.
