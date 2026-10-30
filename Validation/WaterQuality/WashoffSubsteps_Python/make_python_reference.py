"""Reference for VAL-WQ-SUBSTEP-001 from HydroPol2D-Python's grid wash-off kernel.

Reads washoff_inputs.csv (written by run_washoff_substep_check.m) and writes
python_reference.csv. Run with the HydroPol2D-Python venv:

    PYTHONPATH=<HydroPol2D-Python>/src <venv>/bin/python make_python_reference.py

Made with HydroPol2D-Python release/studio-0.9.0rc2 at 66cd12f. MATLAB's C_3 is
per hour (W in kg/h); Python's coefficient is per second, so it gets C_3/3600.
"""

from __future__ import annotations

import csv
from pathlib import Path

import jax

jax.config.update("jax_enable_x64", True)

import numpy as np  # noqa: E402

from hydropol2d.modules.water_quality import mass_based_washoff_step  # noqa: E402

HERE = Path(__file__).resolve().parent
NY, NX = 3, 4

with (HERE / "washoff_inputs.csv").open(newline="") as handle:
    rows = {row["name"]: np.asarray([float(v) for v in row["values"].split()]) for row in csv.DictReader(handle)}

mass = rows["mass_kg"].reshape(NY, NX)
qx = rows["qx_m3_s"].reshape(NY, NX + 1)
qy = rows["qy_m3_s"].reshape(NY + 1, NX)
step = mass_based_washoff_step(
    pollutant_mass_kg=mass,
    qx_m3_s=qx,
    qy_m3_s=qy,
    valid_mask=np.ones((NY, NX), dtype=bool),
    washoff_coefficient=rows["C_3"].reshape(NY, NX) / 3600.0,
    washoff_exponent=float(rows["C_4"][0]),
    minimum_buildup_kg_m2=float(rows["min_Bt"][0]) / 1000.0,
    cell_area_m2=float(rows["cell_area_m2"][0]),
    dt_s=float(rows["dt_s"][0]),
)
with (HERE / "python_reference.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["name", "values"])
    writer.writerow(["mass_kg", " ".join(repr(float(v)) for v in np.asarray(step.pollutant_mass_kg).ravel())])
    writer.writerow(["export_kg", repr(float(step.total_export_kg))])
print("wrote python_reference.csv")
