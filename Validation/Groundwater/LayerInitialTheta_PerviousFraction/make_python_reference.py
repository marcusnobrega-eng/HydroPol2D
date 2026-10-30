"""Reference for VAL-SOIL-THETA0-001 from HydroPol2D-Python's initialize_layered_soil_storage.

Reads layer_inputs.csv (written by run_layer_initial_theta_check.m) and writes
python_reference.csv. Run with the HydroPol2D-Python venv:

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

from hydropol2d.modules import soil  # noqa: E402

HERE = Path(__file__).resolve().parent
with (HERE / "layer_inputs.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
col = {key: np.asarray([float(row[key]) for row in rows]) for key in rows[0]}
profile = SimpleNamespace(
    pervious_fraction=jnp.asarray(col["pervious_fraction"]),
    near_surface_thickness_m=jnp.asarray(col["near_thickness_m"]),
    root_zone_thickness_m=jnp.asarray(col["root_thickness_m"]),
    transmission_thickness_m=jnp.asarray(col["transmission_thickness_m"]),
    near_surface_capacity_mm=jnp.asarray(col["near_capacity_mm"]),
    root_zone_capacity_mm=jnp.asarray(col["root_capacity_mm"]),
    transmission_capacity_mm=jnp.asarray(col["transmission_capacity_mm"]),
)
theta = jnp.stack((
    jnp.asarray(col["theta_near"]), jnp.asarray(col["theta_root"]), jnp.asarray(col["theta_transmission"])
))
state = soil.initialize_layered_soil_storage(
    profile=profile, theta_initial=theta, theta_r=jnp.asarray(col["theta_r"]),
    valid_mask=jnp.ones(col["theta_r"].shape, dtype=bool),
)
with (HERE / "python_reference.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["near_storage_mm", "root_storage_mm", "transmission_storage_mm"])
    for values in zip(
        np.asarray(state.near_surface_storage_mm),
        np.asarray(state.root_zone_storage_mm),
        np.asarray(state.transmission_storage_mm),
    ):
        writer.writerow([repr(float(value)) for value in values])
print("wrote python_reference.csv")
