#!/usr/bin/env python3
"""Run the self-contained HydroPol2D release checks.

Each MATLAB driver runs in a fresh process. This prevents workspace and path
state from one validation case changing another case. The script records the
case-level pass/fail tables produced by each driver and fails only when a
report-ready result is not a pass. Diagnostic rows remain visible in the
summary but do not block a release.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
import sys
import time


@dataclass(frozen=True)
class Check:
    name: str
    driver: str
    pass_file: str | None = None
    profile: str = "full"


CHECKS = (
    Check("bundled TopoToolbox runtime", "", profile="fast"),
    Check("canopy interception", "Validation/Canopy_Interception/SingleCell_Storage_Balance/run_canopy_interception_model.m", "Validation/Canopy_Interception/SingleCell_Storage_Balance/Outputs/Validation/Pass_Fail.csv", "fast"),
    Check("snow", "Validation/Snow/ColdWarmPartition_Melt/run_snow_model.m", "Validation/Snow/ColdWarmPartition_Melt/Outputs/Validation/Pass_Fail.csv", "fast"),
    Check("layered-soil dynamics", "Validation/Infiltration/LayeredSoil_Profile/run_layered_soil_dynamics_test.m", "Validation/Infiltration/LayeredSoil_Profile/Outputs/Validation/Layered_Dynamics_Pass_Fail.csv", "fast"),
    Check("V-tilted infiltration equations", "Validation/Infiltration/VTilted_Infiltration/run_vtilted_infiltration_validation_suite.m", "Validation/Infiltration/VTilted_Infiltration/Outputs/Validation/VTilted_Infiltration_Pass_Fail.csv", "full"),
    Check("V-tilted infiltration full-model dynamics", "Validation/Infiltration/VTilted_Infiltration/run_hydropol2d_vtilted_infiltration_cases.m", profile="full"),
    Check("evapotranspiration", "Validation/Evapotranspiration/Reference_ET_Extraction/run_vtilted_et_validation_suite.m", "Validation/Evapotranspiration/Reference_ET_Extraction/Outputs/Validation/VTilted_ET_Pass_Fail.csv", "fast"),
    Check("asynchronous groundwater", "Validation/Groundwater/AsyncRecharge_VTilted/run_hydropol2d_vtilted_groundwater_cases.m", "Validation/Groundwater/AsyncRecharge_VTilted/Outputs/Validation/VTilted_Groundwater_Pass_Fail.csv", "full"),
    Check("analytical groundwater hillslope", "Validation/Groundwater/Hillslope_Analytical/run_hillslope_analytical_validation.m", "Validation/Groundwater/Hillslope_Analytical/Outputs/Validation/Hillslope_Groundwater_Pass_Fail.csv", "fast"),
    Check("full momentum", "Validation/Hydrodynamics/FullMomentum/run_full_momentum_hydrodynamics_validation.m", "Validation/Hydrodynamics/FullMomentum/Outputs/Validation/FullMomentum_Hydrodynamics_Pass_Fail.csv", "full"),
    Check("local inertial and cellular automata", "Validation/Hydrodynamics/LocalInertial_CA/run_local_inertial_ca_hydrodynamics_validation.m", "Validation/Hydrodynamics/LocalInertial_CA/Outputs/Validation/LocalInertial_CA_Hydrodynamics_Pass_Fail.csv", "full"),
    Check("kinematic and diffusive routing", "Validation/Hydrodynamics/Diffusive_Kinematic/run_diffusive_kinematic_hydrodynamics_validation.m", "Validation/Hydrodynamics/Diffusive_Kinematic/Outputs/Validation/Diffusive_Kinematic_Hydrodynamics_Pass_Fail.csv", "full"),
    Check("non-breaking wave", "Validation/Hydrodynamics/NonBreakingWave/run_nonbreaking_wave_hydrodynamics_validation.m", "Validation/Hydrodynamics/NonBreakingWave/Outputs/Validation/NonBreakingWave_Hydrodynamics_Pass_Fail.csv", "full"),
    Check("reservoir rating curve", "Validation/Reservoir/RatingCurve_StorageDepletion/run_reservoir_rating_curve_validation.m", "Validation/Reservoir/RatingCurve_StorageDepletion/Outputs/Validation/Reservoir_RatingCurve_Pass_Fail.csv", "fast"),
    Check("inflow hydrograph boundary", "Validation/BoundaryConditions/InflowHydrograph_RectangularDomain/run_inflow_hydrograph_validation.m", "Validation/BoundaryConditions/InflowHydrograph_RectangularDomain/Outputs/Validation/InflowHydrograph_Boundary_Pass_Fail.csv", "fast"),
    Check("spatial rainfall", "Validation/SpatialRainfall/TinyRaster_TotalRainfall/run_spatial_rainfall_raster_validation.m", "Validation/SpatialRainfall/TinyRaster_TotalRainfall/Outputs/Validation/SpatialRainfall_Raster_Pass_Fail.csv", "fast"),
    Check("water quality", "Validation/WaterQuality/BuildUpWashOff_MassBalance/run_water_quality_washoff_validation.m", "Validation/WaterQuality/BuildUpWashOff_MassBalance/Outputs/Validation/WaterQuality_Washoff_Pass_Fail.csv", "fast"),
    Check("human risk", "Validation/HumanRisk/DepthVelocity_Classification/run_human_risk_instability_validation.m", "Validation/HumanRisk/DepthVelocity_Classification/Outputs/Validation/Pass_Fail.csv", "fast"),
    Check("Neal (2012) simple channel", "Validation/Subgrid/Neal2012_SimpleStraightChannel/run_neal2012_simple_examples.m", "Validation/Subgrid/Neal2012_SimpleStraightChannel/Outputs/Validation/Pass_Fail.csv", "full"),
    Check("Neal (2012) composite steady channel", "Validation/Subgrid/Neal2012_CompositeChannel_100cms/run_neal_composite_channel_100cms.m", "Validation/Subgrid/Neal2012_CompositeChannel_100cms/Outputs/Validation/Pass_Fail.csv", "full"),
    Check("Neal (2012) composite Nash hydrograph", "Validation/Subgrid/Neal2012_CompositeChannel_100cms/run_neal_composite_channel_nash100.m", "Validation/Subgrid/Neal2012_CompositeChannel_100cms/Outputs/Validation/NashTransient/Pass_Fail.csv", "full"),
    Check("Neal (2012) staged overbank channel", "Validation/Subgrid/Neal2012_StagedOverbank_100m/run_neal_staged_overbank_100m.m", "Validation/Subgrid/Neal2012_StagedOverbank_100m/Outputs/Validation/Pass_Fail.csv", "full"),
    Check("Neal (2012) two-dimensional junction", "Validation/Subgrid/Neal2012_2D_CrossJunction/run_neal_2d_cross_junction_validation.m", "Validation/Subgrid/Neal2012_2D_CrossJunction/Outputs/Validation/Pass_Fail.csv", "full"),
)


def matlab_quote(value: Path | str) -> str:
    return str(value).replace("'", "''")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matlab",
        default=shutil.which("matlab"),
        help="MATLAB executable (default: matlab on PATH)",
    )
    parser.add_argument(
        "--profile",
        choices=("fast", "full"),
        default="full",
        help="fast excludes longer domain tests; full runs the release suite",
    )
    parser.add_argument("--keep-going", action="store_true", help="run remaining checks after a failure")
    return parser.parse_args()


def matlab_runtime_command(root: Path) -> str:
    root_q = matlab_quote(root)
    fn_q = matlab_quote(root / "HydroPol2D_Functions")
    return (
        "restoredefaultpath; "
        f"cd('{root_q}'); "
        f"addpath('{fn_q}','-begin'); "
        f"runtime=hydropol2d_add_runtime_paths('{root_q}'); "
        "assert(startsWith(which('GRIDobj'), runtime.topotoolbox_lite_root));"
    )


def read_pass_file(path: Path) -> tuple[int, int, list[str]]:
    """Return required passes, required failures, and diagnostic descriptions."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return 0, 1, ["empty pass/fail table"]

    required_passes = 0
    required_failures = 0
    diagnostics: list[str] = []
    for row in rows:
        case_id = row.get("case_id", "unnamed")
        status = row.get("status", "").strip().lower()
        report_ready = row.get("report_ready", "").strip().lower()
        is_required = report_ready in {"1", "true", "yes"}
        is_pass = status in {"pass", "passed"} or row.get("passed", "").strip().lower() in {"1", "true"}
        if is_required:
            required_passes += int(is_pass)
            required_failures += int(not is_pass)
        elif status and status not in {"pass", "passed"}:
            diagnostics.append(f"{case_id}: {status}")
    return required_passes, required_failures, diagnostics


def run_check(matlab: str, root: Path, check: Check) -> dict[str, str | int | float]:
    started = time.monotonic()
    if not check.driver:
        expression = matlab_runtime_command(root)
    else:
        driver = root / check.driver
        expression = matlab_runtime_command(root) + f" run('{matlab_quote(driver)}');"

    process = subprocess.run(
        [matlab, "-batch", expression],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    elapsed_s = time.monotonic() - started
    result: dict[str, str | int | float] = {
        "check": check.name,
        "driver": check.driver or "bundled runtime smoke test",
        "exit_code": process.returncode,
        "elapsed_s": round(elapsed_s, 2),
        "required_passes": 0,
        "required_failures": 0,
        "diagnostics": "",
        "status": "failed" if process.returncode else "passed",
    }
    if process.returncode:
        result["diagnostics"] = process.stdout[-3000:].replace("\n", " | ")
        return result

    if check.pass_file:
        pass_file = root / check.pass_file
        if not pass_file.is_file():
            result["status"] = "failed"
            result["required_failures"] = 1
            result["diagnostics"] = f"missing {check.pass_file}"
            return result
        passes, failures, diagnostics = read_pass_file(pass_file)
        result["required_passes"] = passes
        result["required_failures"] = failures
        result["diagnostics"] = "; ".join(diagnostics)
        if failures:
            result["status"] = "failed"
        elif diagnostics:
            result["status"] = "passed_with_diagnostics"
    return result


def main() -> int:
    args = parse_args()
    if not args.matlab:
        print("ERROR: provide --matlab or add MATLAB to PATH", file=sys.stderr)
        return 2
    matlab = str(Path(args.matlab).expanduser())
    if not Path(matlab).is_file():
        print(f"ERROR: MATLAB executable not found: {matlab}", file=sys.stderr)
        return 2

    root = Path(__file__).resolve().parents[2]
    output_dir = root / "Validation" / "Release" / "Outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    selected = [check for check in CHECKS if args.profile == "full" or check.profile == "fast"]
    results: list[dict[str, str | int | float]] = []
    for index, check in enumerate(selected, start=1):
        print(f"[{index}/{len(selected)}] {check.name}", flush=True)
        result = run_check(matlab, root, check)
        results.append(result)
        print(f"  {result['status']} ({result['elapsed_s']} s)", flush=True)
        if result["status"] == "failed" and not args.keep_going:
            break

    fieldnames = ["check", "driver", "status", "exit_code", "elapsed_s", "required_passes", "required_failures", "diagnostics"]
    summary = output_dir / "Validation_Release_Summary.csv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    failed = [result for result in results if result["status"] == "failed"]
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    print(f"\nWrote {summary.relative_to(root)} at {stamp}")
    if failed:
        print(f"{len(failed)} release check(s) failed.", file=sys.stderr)
        return 1
    print(f"{len(results)} release checks completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
