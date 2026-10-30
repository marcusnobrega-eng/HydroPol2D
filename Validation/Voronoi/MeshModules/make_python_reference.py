"""Reference numbers for VAL-MESH-MODULES-001 from the HydroPol2D-Python kernels.

Reads kernel_cells.csv, kernel_edges.csv and kernel_outlets.csv (written by
run_mesh_modules_check.m) and writes python_reference_cells.csv and
python_reference_scalars.csv. Run with the HydroPol2D-Python venv, e.g.

    PYTHONPATH=<HydroPol2D-Python>/src <venv>/bin/python make_python_reference.py

Made with HydroPol2D-Python release/studio-0.9.0rc2 at b5f224b.
"""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402

from hydropol2d.mesh_modules.risk import cell_discharge_operators  # noqa: E402
from hydropol2d.modules.human_risk import force_balance_risk  # noqa: E402
from hydropol2d.modules.water_quality import mass_based_washoff_edge_step  # noqa: E402
from hydropol2d.routing.reservoir import reservoir_outlet_transfers  # noqa: E402

HERE = Path(__file__).resolve().parent
DT_S = 30.0


def read(name: str) -> dict[str, np.ndarray]:
    with (HERE / name).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {key: np.asarray([float(row[key]) for row in rows]) for key in rows[0]}


cells = read("kernel_cells.csv")
edges = read("kernel_edges.csv")
outlets = read("kernel_outlets.csv")

# MATLAB ids are 1-based and use 0 for "outside"; Python uses 0-based and -1.
source = outlets["source"].astype(np.int64) - 1
external = outlets["external"].astype(bool)
destination = np.where(external, 0, outlets["destination"].astype(np.int64) - 1)
volume, external_m3, _transferred, control_release = reservoir_outlet_transfers(
    jnp.asarray(np.maximum(cells["volume_m3"], 0.0)),
    source_cells=jnp.asarray(source),
    destination_cells=jnp.asarray(destination),
    external=jnp.asarray(external),
    coefficient=jnp.asarray(outlets["coefficient"]),
    sill_depth_m=jnp.asarray(outlets["sill_depth_m"]),
    exponent=jnp.asarray(outlets["exponent"]),
    control_index=jnp.asarray(outlets["control"].astype(np.int64) - 1),
    source_area_m2=jnp.asarray(cells["area_m2"][source]),
    control_count=2,
    dt_s=DT_S,
)
volume = np.maximum(np.asarray(volume), 0.0)

owner = edges["owner"].astype(np.int64) - 1
neighbor = edges["neighbor"].astype(np.int64) - 1
quality = mass_based_washoff_edge_step(
    pollutant_mass_kg=cells["mass_kg"],
    edge_discharge_m3_s=edges["q_m2_s"] * edges["length_m"],
    edge_owner=owner,
    edge_neighbor=neighbor,
    washoff_coefficient=cells["coefficient_kg_s"],
    washoff_exponent=cells["exponent"],
    minimum_buildup_kg=cells["threshold_kg"],
    dt_s=DT_S,
)

mesh = SimpleNamespace(
    n_cells=cells["cell"].size,
    n_edges=edges["edge"].size,
    edge_owner=owner,
    edge_neighbor=neighbor,
    edge_length=edges["length_m"],
    edge_midpoint_x=edges["midpoint_x_m"],
    edge_midpoint_y=edges["midpoint_y_m"],
    cell_area=cells["area_m2"],
    cell_x=cells["x_m"],
    cell_y=cells["y_m"],
)
operator_x, operator_y = cell_discharge_operators(mesh)
depth = volume / cells["area_m2"]
hu = operator_x @ edges["q_m2_s"]
hv = operator_y @ edges["q_m2_s"]
speed = np.divide(np.hypot(hu, hv), depth, out=np.zeros_like(depth), where=depth > 1.0e-3)
risk = force_balance_risk(depth_m=depth, velocity_m_s=speed)

with (HERE / "python_reference_cells.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["cell", "volume_m3", "mass_kg", "cell_concentration_mg_l", "risk_index"])
    for index in range(mesh.n_cells):
        writer.writerow([
            index + 1,
            repr(float(volume[index])),
            repr(float(quality.pollutant_mass_kg[index])),
            repr(float(quality.cell_concentration_mg_l[index])),
            repr(float(risk.risk_index[index])),
        ])
with (HERE / "python_reference_scalars.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "external_m3", "control_1_release_m3", "control_2_release_m3",
        "exported_kg", "outlet_concentration_mg_l",
    ])
    writer.writerow([
        repr(float(external_m3)),
        repr(float(control_release[0])),
        repr(float(control_release[1])),
        repr(float(quality.total_export_kg)),
        repr(float(quality.outlet_concentration_mg_l)),
    ])
print("wrote python_reference_cells.csv and python_reference_scalars.csv")
