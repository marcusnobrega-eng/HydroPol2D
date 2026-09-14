#!/usr/bin/env python3
"""Fit Lin hydraulic geometry, or record use of the existing global fallback."""

from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n")


def apply_global_fallback(final_config: Path, defaults: dict) -> None:
    config = json.loads(final_config.read_text())
    config["river-geometry-source"] = "power_law"
    for key in list(config):
        if key.startswith("spatial-"):
            del config[key]
    hbd = defaults["hydrobathydem"]
    config.update(
        {
            "beta-1": hbd["beta_1"],
            "beta-2": hbd["beta_2"],
            "alfa-1": hbd["alfa_1"],
            "alfa-2": hbd["alfa_2"],
        }
    )
    write_json(final_config, config)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--catalog", type=Path, default=HERE / "india_highres_cases.json")
    args = parser.parse_args()

    defaults = json.loads(args.catalog.read_text())["defaults"]
    hbd_defaults = defaults["hydrobathydem"]
    root = args.case_dir / "Build" / "HydroBathyDEM"
    first = root / "FirstPass"
    calibration = root / "Lin2020" / "calibration"
    final_config = root / "hydrobathydem_final.json"
    decision = calibration / "hydraulic_geometry_decision.json"
    threshold = int(hbd_defaults["lin_calibration_threshold_km2"])
    coefficients = [
        calibration / f"D4_beta_1_width_{threshold}km2.tif",
        calibration / f"D4_beta_2_width_{threshold}km2.tif",
        calibration / f"D4_alfa_1_depth_{threshold}km2.tif",
        calibration / f"D4_alfa_2_depth_{threshold}km2.tif",
    ]

    if decision.is_file():
        previous = json.loads(decision.read_text())
        if previous["mode"] == "global_power_law_fallback":
            apply_global_fallback(final_config, defaults)
            print("[SKIP] Reapplying the recorded global hydraulic-geometry fallback.")
            return
        if previous["mode"] == "spatial_coefficients_or_power_law" and all(path.is_file() for path in coefficients):
            print("[SKIP] Existing spatial hydraulic-geometry calibration is complete.")
            return

    calibration.mkdir(parents=True, exist_ok=True)
    command = [
        "hydrobathydem-calibrate-hydraulics",
        "--fac-area", str(first / "d4" / "D4_Wshed_Properties_fac_area_km2.tif"),
        "--d4-direction", str(first / "d4" / "D4_flow_direction.tif"),
        "--lin-gpkg", str(root / "Lin2020" / "processed" / "lin2020_dem_domain_width_depth.gpkg"),
        "--out-dir", str(calibration),
        "--selected-threshold-km2", str(threshold),
        "--fit-area-source", "d4",
        "--application-min-area-km2", str(defaults["d4_min_river_area_km2"]),
        "--max-H-abg-m", str(hbd_defaults["max_h_abg_m"]),
    ]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="")
    record = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "returncode": result.returncode,
        "existing_global_coefficients": {
            "beta_1": hbd_defaults["beta_1"],
            "beta_2": hbd_defaults["beta_2"],
            "alfa_1": hbd_defaults["alfa_1"],
            "alfa_2": hbd_defaults["alfa_2"],
        },
    }
    if result.returncode == 0 and all(path.is_file() for path in coefficients):
        record["mode"] = "spatial_coefficients_or_power_law"
        record["reason"] = "Lin et al. spatial calibration passed the HydroBathyDEM quality criteria."
    elif "Global Lin hydraulic-geometry fit failed" in result.stdout:
        apply_global_fallback(final_config, defaults)
        record["mode"] = "global_power_law_fallback"
        record["reason"] = (
            "The available Lin et al. samples did not pass the predefined calibration criteria; "
            "the existing HydroPol2D global width-area and depth-area relationships were retained."
        )
        record["calibration_output_tail"] = result.stdout[-4000:]
        print("[FALLBACK] Lin calibration was insufficient; retained existing global coefficients.")
    else:
        raise RuntimeError(f"Hydraulic-geometry calibration failed unexpectedly (exit {result.returncode}).")
    write_json(decision, record)


if __name__ == "__main__":
    main()
