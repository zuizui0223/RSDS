"""Known-truth and original-source immutable season/position audit checks."""
from __future__ import annotations
import unittest

from scripts import rsds_six_same_species_q_position_season_audit_v1 as m


class SourcePositionVsSeason(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=m.audit()

    def test_source_sample_and_missing_not_zero(self):
        x=self.r
        self.assertEqual(x["original_previous_source_studies"],63)
        self.assertEqual(x["original_new_source_study_slots"],6)
        self.assertEqual(x["original_source_file_admissible"],4)
        self.assertEqual(x["original_source_file_HOLD"],2)
        self.assertEqual(x["original_q_supported_new"],3)
        self.assertEqual(x["original_daylength_supported_new"],3)
        self.assertEqual(x["new_source_bootstrap_primary_q_sign_verdicts"],{
            "NOMINAL_SUPPORTED_SAME_SIGN":3,
            "UNAVAILABLE_OR_ORIGINAL_SOURCE_HOLD":3,
        })
        for z in x["all_original_new_source_and_prior_q_source_daylength_rows"]:
            if not z["new_q"]["available"]:
                self.assertIsNone(z["new_q"]["mean"])
                self.assertEqual(z["new_q"]["supported"],"UNAVAILABLE")
                self.assertFalse(z["new_q_daylength"]["available"])

    def test_negative_mean_q_persists_but_slope_can_oppositely_sign(self):
        z={x["species"]:x for x in self.r["all_original_new_source_and_prior_q_source_daylength_rows"] if x["new_q"]["available"]}
        self.assertEqual(set(z),{"Anser albifrons","Odocoileus hemionus","Vulpes vulpes"})
        for sp in ("Anser albifrons","Odocoileus hemionus"):
            self.assertEqual(z[sp]["new_q"]["supported"],"NEGATIVE")
            self.assertGreaterEqual(len(z[sp]["prior_season_opposite_new_sign_supported_source_DOIs"]),1)
            self.assertTrue(z[sp]["new_daylength_opposite_prior_nominal_sign"])
        self.assertEqual(z["Anser albifrons"]["new_daylength"]["supported"],"POSITIVE")
        self.assertEqual(z["Odocoileus hemionus"]["new_daylength"]["supported"],"NEGATIVE")
        self.assertEqual(z["Vulpes vulpes"]["new_q"]["supported"],"POSITIVE")
        self.assertEqual(z["Vulpes vulpes"]["new_daylength"]["supported"],"INCLUDES_ZERO")
        self.assertFalse(z["Vulpes vulpes"]["new_daylength_opposite_prior_nominal_sign"])
        self.assertEqual(self.r["species_with_opposite_nominal_daylength_slope_support_across_archived_sources"],["Anser albifrons","Odocoileus hemionus"])

    def test_prior_solar_GPS_exposure_cannot_be_called_untouched(self):
        x=self.r
        self.assertEqual(x["new_same_q_sign_with_known_prior_GPS_exposure_count"],1)
        self.assertEqual(x["new_same_q_sign_with_prior_GPS_exposure_unresolved_count"],2)
        self.assertTrue(x["source_to_source_tag_and_population_independence_unverified"])
        self.assertTrue(x["post_outcome_comparison_not_a_registered_daylength_reversal_test"])
        self.assertTrue(x["source_to_source_slope_difference_interval_not_estimated"])
        self.assertTrue(x["no_common_seasonal_calendar_control"])
        self.assertTrue(x["no_new_GPS_or_original_provider_metadata_access"])

    def test_nominal_CI_boundary(self):
        self.assertEqual(m.check_effect({
            "available":True,"q":0.02,"source_tag_bootstrap95":[-0.01,0.05],
            "original_source_tagged_animals":6,"point_sign":"POSITIVE",
            "nominal_source_CI_sign":"INCLUDES_ZERO",
        })["supported"],"INCLUDES_ZERO")
        with self.assertRaises(ValueError):
            m.check_effect({
                "available":True,"q":0.02,"source_tag_bootstrap95":[-0.01,0.05],
                "original_source_tagged_animals":6,"point_sign":"POSITIVE",
                "nominal_source_CI_sign":"POSITIVE",
            })
        self.assertFalse(m.check_effect({
            "available":False,"q":None,"nominal_source_CI_sign":"UNAVAILABLE"
        })["available"])


if __name__=="__main__":
    unittest.main()
