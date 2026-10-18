# Layered subsurface validation

Python -> MATLAB parity for the layered soil-over-regolith column (HydroPol2D-Python
`modules/subsurface.py`, `core/layered_vertical.py`, `voronoi_layered.py`).

| Script | What it checks |
| --- | --- |
| `run_layered_kernel_parity.m` | Every column kernel against `layered_kernels.json` (written by Python `tools/export_layered_fixtures.py`): geometry, MMF table update, drainage, Darcy flux, Feddes, shared-table bisection, lateral advance. 193 checks, worst 2.8e-16 |
| `run_layered_step_parity.m` | The layered vertical step against `layered_steps.json` (Python `tools/export_layered_step_fixture.py`, 2 x 2 V-tilted, 4 scenarios, 312 steps). 3120 checks, worst 2.1e-17 |
| `run_layered_vtilted_case.m` | Full-model V-tilted runs (`local`, `lateral`, `feddes`, `legacy`) from a temporary copy of `Examples/VTilted_Excel` |
| `run_feddes_legacy_consistency.m` | Feddes in the legacy three-layer ET equals the layered Feddes ET for the same state; Feddes never takes more than the unstressed rule |
| `run_layered_mesh_case.m` | The 20 m V-tilted Voronoi mesh with 2 HRUs: legacy, layered and layered + Feddes; every step closes mass within 1e-8 of the rain volume |

Regenerate the fixtures from the HydroPol2D-Python repository:

```bash
python tools/export_layered_fixtures.py <HydroPol2D>/Validation/Groundwater/Layered/layered_kernels.json
python tools/export_layered_step_fixture.py <HydroPol2D>/Validation/Groundwater/Layered/layered_steps.json
```
