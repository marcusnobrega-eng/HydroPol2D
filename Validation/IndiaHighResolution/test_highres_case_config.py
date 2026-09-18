#!/usr/bin/env python3

import importlib.util
import json
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("highres", ROOT / "Config" / "prepare_india_highres_case.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
HYDRAULICS_SPEC = importlib.util.spec_from_file_location(
    "highres_hydraulics", ROOT / "Config" / "calibrate_india_highres_hydraulics.py"
)
HYDRAULICS = importlib.util.module_from_spec(HYDRAULICS_SPEC)
HYDRAULICS_SPEC.loader.exec_module(HYDRAULICS)
PREFLIGHT_SPEC = importlib.util.spec_from_file_location(
    "highres_preflight", ROOT / "Config" / "preflight_india_highres_case.py"
)
PREFLIGHT = importlib.util.module_from_spec(PREFLIGHT_SPEC)
PREFLIGHT_SPEC.loader.exec_module(PREFLIGHT)
SUMMARY_SPEC = importlib.util.spec_from_file_location(
    "highres_summary", ROOT / "Validation" / "IndiaHighResolution" / "summarize_india_highres_output.py"
)
SUMMARY = importlib.util.module_from_spec(SUMMARY_SPEC)
SUMMARY_SPEC.loader.exec_module(SUMMARY)


class HighResolutionCaseConfigTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = json.loads((ROOT / "Config" / "india_highres_cases.json").read_text())

    def test_all_cases_have_thirty_day_warmup(self):
        defaults = self.catalog["defaults"]
        for case in self.catalog["cases"].values():
            start = MODULE.simulation_start(case, defaults)
            event = datetime.fromisoformat(case["event_start"])
            self.assertEqual((event - start).days, 30)

    def test_memory_scales_with_active_cells(self):
        defaults = self.catalog["defaults"]
        self.assertEqual(MODULE.request_memory_gib(120_128, defaults), 16)
        self.assertEqual(MODULE.request_memory_gib(856_630, defaults), 24)
        self.assertEqual(MODULE.request_memory_gib(2_262_704, defaults), 48)

    def test_hydrobathydem_resolution_settings(self):
        defaults = self.catalog["defaults"]
        config = MODULE.hbd_config(Path("dem.tif"), Path("out"), 90, defaults, spatial=False)
        self.assertEqual(config["channel-cell-width-m"], 90)
        self.assertEqual(config["protect-stream-buffer-m"], 270)
        self.assertEqual(config["breach-dist-cells"], 100)
        self.assertEqual(config["river-geometry-source"], "power_law")

    def test_domain_steps_are_resumable(self):
        command = MODULE.resumable_command("build-domain", [Path("a.tif"), Path("b.tif")])
        self.assertEqual(
            command,
            "if ! ( test -s a.tif && test -s b.tif ); then build-domain; fi; test -s a.tif && test -s b.tif",
        )

    def test_hydraulic_fallback_removes_spatial_coefficients(self):
        with tempfile.TemporaryDirectory() as folder:
            config = Path(folder) / "final.json"
            config.write_text(json.dumps({"river-geometry-source": "spatial", "spatial-beta-1-raster": "x"}))
            HYDRAULICS.apply_global_fallback(config, self.catalog["defaults"])
            updated = json.loads(config.read_text())
        self.assertEqual(updated["river-geometry-source"], "power_law")
        self.assertNotIn("spatial-beta-1-raster", updated)
        self.assertEqual(updated["beta-1"], self.catalog["defaults"]["hydrobathydem"]["beta_1"])

    def test_no_unknown_groundwater_source(self):
        source = self.catalog["sources"]["water_table"]
        self.assertIn("Fan", source["product"])
        self.assertNotIn("GW_table.tif", source["product"])

    def test_effective_soil_depth_has_existing_numerical_floor(self):
        self.assertEqual(self.catalog["defaults"]["minimum_effective_soil_depth_m"], 0.01)

    def test_highres_cases_use_case_specific_etp(self):
        launcher = (ROOT / "HydroPol2D_V115.m").read_text()
        runner = (ROOT / "Config" / "run_india_highres_case.sbatch").read_text()
        variable = "HYDROPOL_DISABLE_ERA5_DAILY_NETCDF"
        self.assertIn(f"getenv('{variable}')", launcher)
        self.assertIn(f"export {variable}=1", runner)

    def test_bypass_paths_support_current_and_legacy_launchers(self):
        bypass = (ROOT / "Config" / "input_paths_bypass.m").read_text()
        self.assertIn("OptionalOverrides", bypass)
        self.assertIn("InputPaths.topo_path", bypass)
        self.assertIn("InputPaths.hydropol2d_tools", bypass)

    def test_groundwater_check_allows_one_float32_ulp(self):
        elevations = np.array([1.0, 2100.0], dtype="float64")
        tolerance = PREFLIGHT.float32_raster_tolerance(elevations)
        self.assertGreaterEqual(tolerance[0], 1e-4)
        self.assertGreaterEqual(tolerance[1], np.spacing(np.float32(2100.0)))

    def test_matlab_serial_dates_are_interpreted(self):
        expected = pd.Timestamp("2005-06-24")
        serial = expected.toordinal() + 366
        actual = SUMMARY.model_dates(pd.Series([serial]), expected)[0]
        self.assertEqual(actual.normalize(), expected)

    def test_depth_raster_timestamp(self):
        actual = SUMMARY.depth_timestamp(Path("Flood_Depths_2018_08_15_12_30_00.tif"))
        self.assertEqual(actual, datetime(2018, 8, 15, 12, 30))

    def test_complete_runs_enable_postprocessing(self):
        runner = (ROOT / "Config" / "run_india_highres_case.sbatch").read_text()
        self.assertIn('float(sys.argv[1]) > 0', runner)
        self.assertIn('unset HYDROPOL_SKIP_POSTPROCESS', runner)
        self.assertIn('HP2D_MAX_TIMESTEP_SECONDS:-60', runner)
        self.assertIn('HP2D_MAX_TIMESTEP_SECONDS:-300', runner)

    def test_india_runtime_profile_uses_current_public_flags(self):
        config = (ROOT / "Config" / "input_data_bypass_script.m").read_text()
        profile = config.split("if ~isempty(getenv('HYDROPOL_INDIA_CASE_ROOT'))", 1)[1]
        self.assertIn("flag_spatial_rainfall = 1", profile)
        self.assertIn("flag_groundwater_modeling = 1", profile)
        self.assertIn("flag_initial_soil_moisture = 1", profile)
        self.assertIn("flag_neal_channel = 1", profile)
        self.assertIn("HYDROPOL_RAINFALL_INTERVAL_MIN", profile)
        self.assertIn("HYDROPOL_RAINFALL_FILENAME_EXAMPLE", profile)
        self.assertNotIn("flag_subgrid = 1", profile)
        self.assertNotIn("flag_overbanks = 1", profile)

    def test_postprocessing_uses_saved_map_count_and_interval_end_times(self):
        postprocessing = (ROOT / "HydroPol2D_Functions" / "post_processing.m").read_text()
        animations = (ROOT / "HydroPol2D_Functions" / "Inundation_Maps.m").read_text()
        self.assertIn("n_saved_maps =", postprocessing)
        self.assertIn("map_time_records = running_control.time_records(2:end)", postprocessing)
        self.assertIn("for i = 1:n_saved_maps", postprocessing)
        self.assertNotIn("for i = 1:length(running_control.time_records)", postprocessing)
        self.assertNotIn("length(running_control.time_records)", animations)
        self.assertNotIn("fullfile('Temporary_Files'", animations)
        self.assertEqual(animations.count("max(max(y_grid)) zmin zmax"), 1)
        self.assertIn("flags.flag_overbanks == 1", postprocessing)
        self.assertIn("flags.flag_overbanks == 1", animations)
        self.assertIn("-c:v mpeg4", animations)
        self.assertIn("double(gather(Maps.Hydro.velocity", postprocessing)
        self.assertIn("double(gather(Maps.Hydro.hazard_dv", postprocessing)
        self.assertIn("double(gather(Maps.Hydro.I_t", postprocessing)
        self.assertIn("gather(HydroMaps.GWdepth_save", postprocessing)

    def test_smoke_run_can_exercise_official_postprocessing(self):
        runner = (ROOT / "Config" / "run_india_highres_case.sbatch").read_text()
        self.assertIn("HP2D_POSTPROCESS_SMOKE", runner)

    def test_summary_reads_actual_tables_folder(self):
        summary = (ROOT / "Validation" / "IndiaHighResolution" / "summarize_india_highres_output.py").read_text()
        self.assertIn('"Tables_CSV" / "Rating_Curve_Gauges.csv"', summary)


if __name__ == "__main__":
    unittest.main()
