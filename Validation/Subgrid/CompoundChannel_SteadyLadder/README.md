# P1-SUBGRID-003: Compound-channel steady-flow ladder

Purpose: compare a 10 m ordinary local-inertial numerical reference against a 30 m lookup-subgrid local-inertial model using a longer, quasi-steady compound-channel benchmark. This isolates subgrid/coarsening effects from equation-set differences.

Geometry:
- Length: `900 m`
- Longitudinal slope: `0.5%`
- Inbank width: `10 m`
- Bank height: `1 m`
- Overbanks: `40 m` each side
- Overbank lateral slope: `1%`
- Manning roughness: `0.035`

Flow ladder:
- `Q = [5, 10, 15, 18, 22, 30, 35] m3/s`
- Internal gauges at `x = 600 m` and `x = 750 m`
- Each flow is initialized independently from its own normal-depth estimate and run until storage changes by `<0.1%` over a 10 min window or until `180 min`.

Current-code result: diagnostic failed. The lookup-subgrid branch now uses a SFINCS-style local-inertial closure based on grid-average face depth, wet fraction, and effective roughness rather than hydraulic-radius conveyance. The strict ladder still passes `0/7` flows. The lowest inbank flow, `Q = 5 m3/s`, has very small gauge-stage error but misses the profile and storage thresholds. At `Q = 10 m3/s`, both runs remain finite and converged, but the 30 m subgrid model under-stores water and misses the internal stage/profile target. From `Q = 15 m3/s` upward, the 10 m ordinary local-inertial reference does not produce a valid steady benchmark within the maximum runtime, so those rungs are retained as diagnostic instability/coarse-representation evidence rather than report-ready subgrid accuracy evidence.

Interpretation: the previous 30 m local-inertial versus 10 m full-momentum comparison was unfair. After correcting the benchmark and replacing the hydraulic-radius face closure with the SFINCS-style effective-depth/roughness closure, the formula/unit tests remain sound, but the compound-channel ladder still is not validated. The result indicates that the remaining problem is not just the face friction closure. The next repairs should focus on local-inertial stability/mass behavior for compound overbank flow, grid/channel alignment, face connectivity through compound sections, outlet-face definition, and whether an explicit 1D/channel-network component is needed for narrow inbank flow at coarse resolution.

Key current metrics:
- Passed flows: `0/7`
- `Q = 5 m3/s`: near-pass; depth RMSE `0.006 m`, profile RMSE `0.057 m`, storage error `-7.34%`.
- `Q = 10 m3/s`: finite/converged but failed; depth RMSE `0.672 m`, profile RMSE `0.648 m`, storage error `-56.32%`.
- `Q >= 15 m3/s`: diagnostic only; the 10 m local-inertial reference does not provide a valid steady truth source.
