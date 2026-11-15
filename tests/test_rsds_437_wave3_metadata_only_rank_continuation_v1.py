"""Known-truth source-accounting tests for Wave3 metadata-only source selection."""
import unittest

from scripts.rsds_437_wave3_metadata_only_rank_continuation_v1 import freeze, select


class Wave3MetadataOnlySourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = freeze()
        cls.result = select()

    def test_public_source_denominators_and_rank_only(self):
        r = self.result
        self.assertEqual(r["source_frame_packages"], 437)
        self.assertEqual(r["source_frame_GPS_packages"], 259)
        self.assertEqual(r["source_frame_metadata_priority_candidates"], 120)
        self.assertEqual(r["prior_selected_source_studies"], 35)
        self.assertEqual(r["wave3_new_original_source_studies"], 16)
        self.assertEqual(r["selected_original_source_studies_after_wave3"], 51)
        self.assertEqual(r["remaining_metadata_priority_single_native_sources"], 69)
        self.assertEqual(r["wave3_selected_single_native_rank_range_inclusive"], [25, 40])
        self.assertEqual(len(r["all_wave3_selected_source_metadata_in_SHA_rank_order"]), 16)

    def test_no_resampling_or_new_animal_outcome(self):
        r = self.result
        rows = r["all_wave3_selected_source_metadata_in_SHA_rank_order"]
        self.assertEqual([z["source_rank_1_based"] for z in rows], list(range(25, 41)))
        self.assertEqual(len({z["package_uuid"] for z in rows}), 16)
        self.assertEqual(len({z["original_repository_study_id"] for z in rows}), 16)
        self.assertEqual(len({z["original_DOI_labels"][0] for z in rows}), 16)
        self.assertTrue(all(z["new_provider_file_metadata_opened"] is False for z in rows))
        self.assertTrue(all(z["new_GPS_animal_events_opened"] == 0 for z in rows))
        self.assertTrue(all(len(z["scientific_name_labels"]) == 1 for z in rows))
        self.assertFalse(r["wave3_new_original_provider_file_metadata_opened"])
        self.assertEqual(r["wave3_new_original_animal_GPS_events_opened"], 0)

    def test_fixed_historical_negative_and_holds(self):
        r = self.result
        self.assertTrue(r["previous_48_source_solar_direction_external_nonreplication_UNCHANGED"])
        self.assertTrue(r["previous_12_source_ERA5_raw_SHA_HOLD_UNCHANGED"])
        self.assertFalse(r["global_biological_probability_sample"])
        self.assertTrue(r["quota_selected_after_two_prior_waves"])
        self.assertIn("Homo sapiens", self.contract["human_privacy_rule"])
        self.assertTrue(all(
            row["human_privacy_file_and_GPS_HOLD"] ==
            ("Homo sapiens" in row["scientific_name_labels"])
            for row in r["all_wave3_selected_source_metadata_in_SHA_rank_order"]
        ))


if __name__ == "__main__":
    unittest.main()
