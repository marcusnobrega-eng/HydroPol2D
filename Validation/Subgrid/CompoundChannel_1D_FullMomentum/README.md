# P1-SUBGRID-002: Compound-channel 1D full-momentum stress test

Purpose: compare a 10 m explicit full-momentum compound-channel simulation against a 30 m full-momentum lookup-subgrid simulation.

Geometry:
- Longitudinal slope: `0.5%`
- Inbank width: `10 m`
- Bank height: `1 m`
- Overbanks: `40 m` each side
- Overbank lateral slope: `1%`
- Manning roughness: `0.035`

Forcing:
- Triangular upstream inflow hydrograph with `35 m3/s` peak.
- Downstream right boundary opened with normal-flow outlet slope `0.005`.

Current-code result: diagnostic failed. The 30 m lookup-subgrid run matches peak magnitude and cumulative outlet volume reasonably, but drains too aggressively and under-predicts final storage. This indicates that full-momentum lookup subgrid needs face-area/effective-depth conveyance and outlet coupling beyond the present volume-conservative representative-depth approximation.

Key metrics:
- Outlet hydrograph NSE: `0.7932`
- Outlet RMSE: `5.026 m3/s`
- Peak error: `-3.44%`
- Peak timing error: `-0.76 min`
- Cumulative outlet volume error: `3.96%`
- Final storage error: `-97.29%`
