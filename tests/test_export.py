"""Actual output-format regression; no optimizer, sampling, or writer mocking."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ce3d_export_check import run_export_self_check


class ExportRoundTripTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = run_export_self_check()

    def test_real_formats_for_all_model_variants(self):
        self.assertEqual(self.report["status"], "passed")
        self.assertEqual(self.report["formats"], ["HTML", "PNG", "SVG", "CSV", "JSON", "NPZ"])
        cases = self.report["cases"]
        self.assertEqual({(case["model"], case["satellites"]) for case in cases},
                         {("LA", False), ("LA", True), ("Voigt", False), ("Voigt", True)})
        for case in cases:
            with self.subTest(model=case["model"], satellites=case["satellites"]):
                self.assertEqual(case["parameter_count"], 15 if case["satellites"] else 11)
                self.assertEqual(case["component_count"], 12 if case["satellites"] else 10)
                self.assertEqual(case["exported_sample_rows"], 24)
                self.assertEqual(case["chain_steps"], 8)
                self.assertEqual(case["chain_blocks"], 2)
                self.assertTrue(case["svg_xml_valid"])
                self.assertTrue(case["area_conservation_passed"])
                self.assertTrue(case["result_quantiles_match_fixture"])
                self.assertTrue(case["primary_report_generated"] and case["organized_layout"])
                self.assertTrue({"Results.html", "figures/background.png", "details/prepared_spectrum.csv",
                                 "figures/fit_ce3d.png", "figures/fit_ce3d.svg", "details/fit_result.json",
                                 "details/posterior_samples_ce3d.csv", "tables/posterior_parameters.csv",
                                 "tables/peak_table.csv", "tables/fit_curves.csv", "details/chain/manifest.json",
                                 "details/chain/last_state.npz", "details/chain/random_state.json"}.issubset(case["files_checked"]))

    def test_fixture_scope_and_cleanup(self):
        self.assertFalse(self.report["sampling_performed"])
        self.assertIn("not MCMC samples or measured data", self.report["fixture"])
        self.assertTrue(self.report["temporary_files_removed"])


if __name__ == "__main__":
    unittest.main()
