# Mesh diffusive flux capped at critical flow (VAL-MESH-DIFF-CRIT-001)

Purpose: on steep, smooth ground the mesh diffusive wave passes no more than
a critical section can.

When `critical_flow` is on, HydroPol2D-Python caps the mesh diffusive-wave
face discharge at the critical discharge h·√(g·h) (`voronoi.py`). MATLAB had
no cap (audit L24). `Voronoi_Diffusive_Face_Flux` and
`Voronoi_Diffusive_Step` now take `critical_flow`. `HydroPol2D_Voronoi_Run`
passes the case's setting. Without it the flux is unchanged.

The case is two 10 m cells with 0.2 m of water on a 1 m drop and Manning n
0.01. Manning's discharge per width (2.37 m²/s) is about 8 times the critical
discharge (0.28 m²/s).

## Run

```matlab
run('Validation/Voronoi/DiffusiveCriticalCap/run_diffusive_critical_cap_check.m')
```

Acceptance: the capped flux equals h·√(g·h) within 1e-12, and the uncapped
flux equals Manning's.
