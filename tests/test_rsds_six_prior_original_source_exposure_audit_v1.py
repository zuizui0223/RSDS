"""Known-truth original source exposure: prior GPS != untouched new q."""
import unittest

from scripts import rsds_six_prior_original_source_exposure_audit_v1 as m


class SourceExposureAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=m.evaluate()

    def test_existing_source_exposure_not_reclassified_as_untouched(self):
        a={x["original_DOI"]:x for x in
           self.result["original_followup_source_exposure_rows"]}
        self.assertEqual(len(a),6)
        self.assertEqual(a["10.5441/001/1.72hh609t"]["prior_data_exposure_status"],
                         "KNOWN_PREVIOUS_ORIGINAL_GPS_DATA_ANALYZED_FOR_SOLAR_DIRECTION")
        self.assertEqual(a["10.5441/001/1.8dc0v84m"]["prior_data_exposure_status"],
                         "KNOWN_PREVIOUS_SOURCE_METADATA_HOLD_NO_GPS_FROM_THAT_ROUTE")
        self.assertEqual(self.result["known_previous_GPS_analysis_sources"],1)
        self.assertEqual(self.result["known_previous_original_provider_metadata_HOLD_sources"],1)
        self.assertEqual(self.result["no_match_in_eight_checked_previous_outputs_sources"],4)

    def test_never_claim_global_independence_from_eight_outputs(self):
        self.assertEqual(self.result["number_of_original_prior_archive_surfaces_inspected"],8)
        self.assertFalse(self.result["global_prior_GPS_exposure_audit_completed"])
        self.assertFalse(self.result["same_animal_independence_verified"])
        self.assertTrue(self.result["no_original_source_GPS_animal_events_reopened"])
        self.assertTrue(all(x["new_original_GPS_rows_reopened"]==0 for x in
                            self.result["original_followup_source_exposure_rows"]))
        self.assertTrue(all(x["prior_other_uninspected_outputs_and_tag_overlap_unresolved"]
                            for x in self.result["original_followup_source_exposure_rows"]))
        self.assertFalse(self.result["prior_original_48_source_solar_direction_replicated"])
        self.assertTrue(self.result["prior_12_ERA5_raw_SHA_HOLD"])

    def test_exact_prior_source_result_case_has_no_false_independence(self):
        fox={"dois":["doi:10.5441/001/1.72hh609t"],
             "source_state":"original_MD5_verified_source_analyzed",
             "original_provider_checksum_verified":{"provider_checksum_verified":True},
             "source_original_coverage":{"original_total_rows":95413}}
        old=m.matching_original_source_nodes({"all_preselected_source_originals":[fox]},
                                             "10.5441/001/1.72hh609t")
        self.assertEqual(len(old),1)
        self.assertTrue(old[0]["provider_original_checksum_verified"])
        self.assertEqual(old[0]["original_source_coverage_rows_reported"],95413)

if __name__=="__main__":
    unittest.main()
