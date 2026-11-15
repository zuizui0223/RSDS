"""Known-truth source-native time niche position vs seasonal slope safety checks."""
import unittest
from scripts import rsds_same_species_position_season_slope_v1 as m

class NicheAxisTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output=m.analyze()

    def test_distinct_source_q_position_and_slope(self):
        z=self.output
        self.assertEqual(z["earlier_source_DOI"],63)
        self.assertEqual(z["previous_source_native_taxon_rows"],88)
        self.assertEqual(z["new_original_source_study_IDs"],6)
        self.assertEqual(z["new_original_source_q_respondents"],3)
        self.assertEqual(z["new_original_source_daylength_respondents"],3)
        self.assertEqual(z["new_q_same_all_prior_supported_direction"],3)
        self.assertEqual(z["new_supported_daylength_opposite_to_any_prior_supported_source"],2)

    def test_new_source_species_and_supported_opposite_slope(self):
        r={x["scientific_name"]:x for x in self.output["original_new_source_axis_rows_including_HOLD"]
           if x["new_native_q"]["available"]}
        self.assertEqual(set(r),{"Anser albifrons","Odocoileus hemionus","Vulpes vulpes"})
        self.assertTrue(all(x["source_niche_position_same_nominal_supported_direction"] for x in r.values()))
        self.assertTrue(r["Anser albifrons"]["new_slope_supported_opposite_to_at_least_one_old_source"])
        self.assertTrue(r["Odocoileus hemionus"]["new_slope_supported_opposite_to_at_least_one_old_source"])
        self.assertFalse(r["Vulpes vulpes"]["new_slope_supported_opposite_to_at_least_one_old_source"])

    def test_barnacle_goose_unavailable_not_a_biological_negative(self):
        rows=self.output["original_new_source_axis_rows_including_HOLD"]
        barn=next(x for x in rows if x["scientific_name"]=="Branta leucopsis")
        self.assertFalse(barn["new_native_q"]["available"])
        self.assertFalse(barn["source_niche_position_same_nominal_supported_direction"])

    def test_no_cross_study_confidence_or_independence_reclassification(self):
        z=self.output
        self.assertTrue(z["prior_GPS_source_independence_not_verified"])
        self.assertTrue(z["original_tagged_animal_nominal95_intervals_not_crosssource_tests"])
        self.assertTrue(z["post_outcome_and_species_target_adaptive"])
        self.assertFalse(z["original_raw_GPS_reopened"])
        self.assertTrue(z["old_48_source_external_solar_direction_not_replicated"])
        self.assertTrue(z["old_12_source_ERA5_raw_SHA_HOLD"])

    def test_bootstrap_direction_requires_ci_excluding_zero(self):
        self.assertIsNone(m.supported_sign({"available":True,"nominal_source_CI_sign":"INCLUDES_ZERO"}))
        self.assertEqual(m.supported_sign({"available":True,"nominal_source_CI_sign":"NEGATIVE"}),"NEGATIVE")
        self.assertIsNone(m.supported_sign({"available":False,"nominal_source_CI_sign":"NEGATIVE"}))

if __name__=="__main__":
    unittest.main()
