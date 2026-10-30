# Cell velocity uses all four faces (VAL-HYDRO-VELFACE-001)

Purpose: check that `refreshing_timestep.m` gives each cell the velocity on
its west face.

The west face of cell j is the east face of cell j−1. MATLAB built the west
face from the east face of cell j itself:
`vel_right(:,2:end)` where it should be `vel_right(:,1:end-1)`. So the west
face dropped out of four results:

- the cell velocity,
- the velocity raster,
- the wave-speed time-step limit,
- human risk mode 1.

HydroPol2D-Python uses all four faces (`routing/timestep.py`). The south face
was already right.

The case is a 3 x 4 grid of 10 m cells at 0.2 m depth. Only the east faces of
column 1 carry 0.5 m3/s, into column 2, so column 2 must see 0.25 m/s on its
west face.

## Run

```matlab
run('Validation/Hydrodynamics/FaceVelocity_WestFace/run_west_face_velocity_check.m')
```

Acceptance: the velocity raster, the x component and the west-face velocity
are exact (errors below 1e-12 m/s). On the old code column 2 read 0 m/s, an
error of 0.25 m/s.
