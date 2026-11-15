"""Source-only, no-original-GPS six pre-event fixed follow-up admission tests."""
import unittest
from scripts import rsds_six_same_species_followup_source_only_v1 as m

class TargetSourceOnlyFreeze(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=m.evaluate()

    def test_all_six_and_no_resampling(self):
        x=self.result
        self.assertEqual(x["original_source_metadata_eligible_six"],6)
        self.assertEqual(x["original_previous_q_source_studies_EXCLUDED"],63)
        self.assertEqual(len(x["original_target_metadata_candidates_eligible"]),6)
        self.assertEqual(len({r["original_repository_study_id"] for r in x["original_target_metadata_candidates_eligible"]}),6)
        self.assertEqual(x["eligible_by_target_taxon"],{
            "Anas platyrhynchos":1,"Anser albifrons":2,
            "Branta leucopsis":1,"Odocoileus hemionus":1,"Vulpes vulpes":1
        })
        self.assertEqual(x["target_names_with_no_new_source"],["Ovis canadensis"])

    def test_no_original_animal_GPS_or_provider_files_before_SHA_freeze(self):
        x=self.result
        self.assertTrue(x["no_original_provider_file_metadata_or_GPS_opened"])
        self.assertTrue(all(z["original_animal_GPS_events_opened"]==0
                            and z["provider_original_file_metadata_opened"] is False
                            for z in x["original_target_metadata_candidates_eligible"]))
        self.assertFalse(x["source_eligibility_targeted_after_original_q_outcomes"] is False)
        self.assertFalse(x["individual_population_independence_verified"])
        self.assertFalse(x["other_RSDS_ODSP_prior_original_GPS_source_overlap_verified"])
        self.assertFalse(x["no_global_species_confirmation_or_meta_probability_sample"] is False)

    def test_source_only_rejects_other_target_candidates_instead_of_substituting(self):
        x=self.result
        self.assertEqual(x["other_noneligible_target_metadata_count"],
                         len(x["all_other_original_target_metadata_candidates_noneligible"]))
        self.assertTrue(all(z["source_only_reasons"] for z in
                            x["all_other_original_target_metadata_candidates_noneligible"]))
        self.assertEqual(x["original_source_metadata_eligible_six"],6)
        self.assertTrue(x["old_12_ERA5_raw_SHA_HOLD"])
        self.assertFalse(x["old_solar_48_source_external_signed_direction_replicated"])

if __name__=="__main__":
    unittest.main()
