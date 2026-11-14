"""Frozen source/unit/dependence/missingness tests for DOI meta sensitivity."""
from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import rsds_doi_meta_heterogeneity_missingness_v1 as m


class OriginalDOIMetaAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.input = m.read_pinned()
        cls.result = m.analyze(cls.input)

    def test_pin_and_reproducibility(self):
        self.assertEqual(m.analyze(m.read_pinned()), self.result)
        self.assertEqual(self.result["source_git_blob_sha"],
                         "af0574c6bb6383402a97b5bebadb8cf0408e9740")

    def test_tampered_original_bytes_fail_closed(self):
        with tempfile.TemporaryDirectory() as t:
            f = Path(t) / "changed.json"
            f.write_text('{"fake": "source"}', encoding="utf-8")
            with patch.object(m, "SOURCE", f):
                with self.assertRaisesRegex(ValueError, "Git blob changed"):
                    m.read_pinned()

    def test_source_not_native_or_tag_is_statistical_unit(self):
        self.assertEqual(self.result["original_studies"], 12)
        self.assertEqual(self.result["original_native_rows"], 22)
        for field in ("q", "daylength"):
            e = self.result["endpoints"][field]
            self.assertEqual(len(e["study_details_rho0"]), 12)
            self.assertEqual(e["models"]["assumed_within_DOI_rho0"]["studies"], 9)
            self.assertEqual(e["missingness"]["unsupported_native_rows"], 7)
            self.assertEqual(e["missingness"]["fully_unsupported_source_DOI"], 3)
            self.assertEqual(sum(s["status"]=="HOLD" for s
                                 in e["study_details_rho0"]), 3)

    def test_no_native_zero_fill_or_nonreplication_reversal(self):
        corrupted = copy.deepcopy(self.input)
        row = next(z for z in corrupted["all_original_study_native_taxon_rows_including_HOLD"]
                   if z["q"]["mean"] is None)
        row["q"]["mean"] = 0.0
        with self.assertRaises((ValueError, TypeError)):
            m.analyze(corrupted)
        corrupted = copy.deepcopy(self.input)
        corrupted["external_48_source_negative_shortwave_direction_nonreplication_UNCHANGED"] = "REPLICATED"
        with self.assertRaisesRegex(ValueError, "nonreplication"):
            m.analyze(corrupted)

    def test_duplicate_pseudoreplication_is_blocked(self):
        corrupted = copy.deepcopy(self.input)
        corrupted["all_original_study_native_taxon_rows_including_HOLD"][-1] = (
            copy.deepcopy(corrupted["all_original_study_native_taxon_rows_including_HOLD"][0]))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            m.analyze(corrupted)

    def test_covariance_extremes_monotone_for_each_study(self):
        for field in ("q", "daylength"):
            e = self.result["endpoints"][field]
            for zero, one in zip(e["study_details_rho0"], e["study_details_rho1"]):
                self.assertEqual(zero["doi"], one["doi"])
                self.assertEqual(zero["estimate"], one["estimate"])
                if zero["status"] == "SUPPORTED":
                    self.assertLessEqual(zero["se_proxy"], one["se_proxy"] + 1e-12)
                else:
                    self.assertIsNone(zero["se_proxy"])

    def test_nine_study_q_meta_diagnostic(self):
        e = self.result["endpoints"]["q"]
        a = e["models"]["assumed_within_DOI_rho0"]
        b = e["models"]["assumed_within_DOI_rho1"]
        self.assertAlmostEqual(a["study_equal_mean"], .010262345757767848, places=12)
        self.assertAlmostEqual(a["random_effects_mean"], .009363575063856868, places=9)
        self.assertAlmostEqual(a["tau2"], .010399088485455489, places=9)
        self.assertAlmostEqual(b["random_effects_mean"], .00941087529419987, places=9)
        for o in (a,b):
            self.assertLess(o["model_95_interval"][0], 0)
            self.assertGreater(o["model_95_interval"][1], 0)
            self.assertLess(o["model_95_prediction_interval"][0], 0)
            self.assertGreater(o["model_95_prediction_interval"][1], 0)
        self.assertEqual(e["source_sign_counts"], {"positive":4, "negative":5})

    def test_nine_study_daylight_meta_diagnostic(self):
        e = self.result["endpoints"]["daylength"]
        a = e["models"]["assumed_within_DOI_rho0"]
        self.assertAlmostEqual(a["study_equal_mean"], -.006383247914400234, places=12)
        self.assertAlmostEqual(a["random_effects_mean"], -.0023006808972531136, places=9)
        self.assertAlmostEqual(a["tau2"], .00012078142534279878, places=9)
        self.assertLess(a["model_95_interval"][0],0)
        self.assertGreater(a["model_95_interval"][1],0)
        self.assertEqual(e["source_sign_counts"], {"positive":5, "negative":4})

    def test_tipping_points_are_conditional_not_observed_effects(self):
        q = self.result["endpoints"]["q"]["missingness"]
        d = self.result["endpoints"]["daylength"]["missingness"]
        self.assertAlmostEqual(q["missing_source_equal_common_effect_at_zero_pooled_mean"],
                               -.030787037273303546,places=12)
        self.assertAlmostEqual(d["missing_source_equal_common_effect_at_zero_pooled_mean"],
                               .0191497437432007,places=12)
        qh=q["all_registered_12_studies_hypothetical_native_range"]
        self.assertLess(qh["all_missing_native_at_observed_min"],0)
        self.assertGreater(qh["all_missing_native_at_observed_max"],0)
        self.assertNotEqual(qh["observed_native_min"],0)
        self.assertNotEqual(qh["observed_native_max"],0)


if __name__ == "__main__":
    unittest.main()
