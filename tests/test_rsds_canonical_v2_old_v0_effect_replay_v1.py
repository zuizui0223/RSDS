"""Regression tests: 36 historical source context effects and 18 source means."""
from __future__ import annotations
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.compare_rsds_canonical_v2_old_v0_biology_v1 import (
    FILES, ROOT, run, sign_interval
)


class CanonicalERA5OldNewSourceModelReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=run()

    def test_thirty_six_context_coefficient_replay_machine_precision(self):
        r=self.report
        self.assertEqual(r["original_field_studies"],2)
        self.assertEqual(r["original_source_taxon_groups"],6)
        self.assertEqual(r["original_source_weather_contexts"],12)
        self.assertEqual(r["native_original_context_and_three_axis_comparisons"],36)
        self.assertLess(r["max_abs_new_vs_old_original_context_fitted_coefficient_difference"],1e-12)
        self.assertTrue(r["all_original_36_context_fitted_coefficients_reproduced_to_1e-12"])

    def test_eighteen_first_source_tag_group_mean_coefficients_replay(self):
        r=self.report
        self.assertEqual(r["native_original_source_animal_mean_comparisons"],18)
        self.assertLess(r["max_abs_new_vs_old_original_source_tag_means_difference"],1e-12)
        self.assertTrue(r["all_original_18_source_animal_mean_coefficients_reproduced_to_1e-12"])

    def test_same_tag_paired_q_and_sun_sin_are_not_fictitious_high_low_radiation(self):
        r=self.report["paired_original_same_tag_context_nominal_descriptive_bootstrap_sign_counts"]
        self.assertEqual(r["q"],{"NOMINAL_POSITIVE":2,"NOMINAL_NEGATIVE":0,
                                  "INCLUDES_ZERO":2,"NOT_ESTIMABLE":2})
        self.assertEqual(r["sun_sin"],{"NOMINAL_POSITIVE":1,"NOMINAL_NEGATIVE":0,
                                        "INCLUDES_ZERO":3,"NOT_ESTIMABLE":2})
        self.assertEqual(r["sun_cos"],{"NOMINAL_POSITIVE":0,"NOMINAL_NEGATIVE":0,
                                        "INCLUDES_ZERO":4,"NOT_ESTIMABLE":2})
        self.assertEqual(len(self.report["all_six_source_taxon_comparisons"]),6)
        self.assertTrue(any("not higher/lower irradiance" in k.lower()
                            for k in self.report["interpretation_ceiling"]))

    def test_old_physical_arrays_not_claimed_identical(self):
        r=self.report
        self.assertTrue(r["historic_12_raw_JSON_SHA_source_HOLDS_UNCHANGED"])
        self.assertEqual(r["old_v0_weather_physical_array_identity"],
                         "NOT_IDENTIFIABLE_FROM_UNARCHIVED_HISTORIC_HOURLY_ARRAYS")
        self.assertEqual(r["historical_external_negative_radiation_nonreplication_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")

    def test_source_group_has_no_original_tag_ID_or_location_leak(self):
        s=json.dumps(self.report)
        for word in ('"individual_id"','"animal_id"','"original_tag_id"',
                     '"location-lat"','"location-long"','"timestamp"',
                     '"gps_fixes"','"source_gridyear"'):
            self.assertNotIn(word,s)

    def test_mismatching_original_first_effect_blob_fails_closed(self):
        import scripts.compare_rsds_canonical_v2_old_v0_biology_v1 as mod
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for rel in FILES:
                tgt=p/rel
                tgt.parent.mkdir(parents=True,exist_ok=True)
                tgt.write_bytes((ROOT/rel).read_bytes())
            f=p/list(FILES)[1]
            x=json.loads(f.read_text())
            x["original_studies"]=3
            f.write_text(json.dumps(x))
            # FILES paths are rooted at global ROOT in load_pinned(), so
            # test direct SHA implementation separately with patched root.
            with patch.object(mod,"ROOT",p):
                with self.assertRaisesRegex(ValueError,"Prior original data result changed"):
                    mod.run()

    def test_known_ci_sign_checks(self):
        self.assertEqual(sign_interval([.01,.15]),"NOMINAL_POSITIVE")
        self.assertEqual(sign_interval([-.15,-.01]),"NOMINAL_NEGATIVE")
        self.assertEqual(sign_interval([-.1,.05]),"INCLUDES_ZERO")
        self.assertEqual(sign_interval(None),"NOT_ESTIMABLE")


if __name__=="__main__":
    unittest.main()
