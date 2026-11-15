"""Frozen six already-discovered species and no-new-outcome source selection checks."""
import unittest
from unittest.mock import patch

from scripts import rsds_six_species_untouched_original_source_selection_v1 as m


class SixOriginalSpeciesMetadataFirstSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=m.select()

    def test_exact_six_source_and_six_targets(self):
        r=self.result
        self.assertEqual(r["previous_study_DOI_count"],63)
        self.assertEqual(r["catalog_original_packages"],437)
        self.assertEqual(len(r["fixed_six_previous_native_q_repeat_species"]),6)
        self.assertEqual(r["new_source_candidates_with_same_taxa_even_if_metadata_ineligible"],16)
        self.assertEqual(r["new_original_metadata_qualified_source_studies"],6)
        self.assertEqual(r["remaining_metadata_unqualified_or_inapplicable_candidates"],10)
        self.assertEqual(r["new_qualified_source_counts_by_six_species"],{
            "Anas platyrhynchos":1,
            "Anser albifrons":2,
            "Branta leucopsis":1,
            "Odocoileus hemionus":1,
            "Ovis canadensis":0,
            "Vulpes vulpes":1,
        })

    def test_original_DOI_study_id_and_first_output_order(self):
        items=self.result["new_original_source_metadata_selection_rows"]
        self.assertEqual(len(items),6)
        self.assertEqual([x["source_index"] for x in items],list(range(6)))
        self.assertEqual(len({x["original_source_study_id"] for x in items}),6)
        self.assertEqual(len({m.normalize_doi(x["original_source_DOI"]) for x in items}),6)
        self.assertEqual(len({x["original_package_UUID"] for x in items}),6)
        self.assertTrue(all(x["source_status"]=="FROZEN_NEW_SOURCE_ID_FOR_FILE_PREFLIGHT" for x in items))
        self.assertTrue(all(x["original_license_labels"]==["CC0 1.0 Universal"] for x in items))
        self.assertTrue(all("GPS" in x["original_sensor_labels"] or "GNSS" in x["original_sensor_labels"] for x in items))
        self.assertTrue(all(x["reported_animal_count_metadata_only"]>=10 for x in items))
        self.assertTrue(all(x["reported_location_count_metadata_only"]>=10000 for x in items))
        self.assertTrue(all(x["original_GPS_animal_events_opened"]==0 for x in items))

    def test_all_six_species_kept_including_zero_opportunity(self):
        self.assertIn("Ovis canadensis",self.result["fixed_six_previous_native_q_repeat_species"])
        self.assertEqual(self.result["new_qualified_source_counts_by_six_species"]["Ovis canadensis"],0)

    def test_outcome_adaptive_target_taxa_but_metadata_only_source_selection(self):
        r=self.result
        self.assertTrue(r["taxon_selection_was_after_prior_biological_outcomes"])
        self.assertTrue(r["no_source_was_selected_based_on_new_animal_GPS_or_q_sign"])
        self.assertFalse(r["new_provider_file_metadata_opened"])
        self.assertEqual(r["new_GPS_rows_opened"],0)
        self.assertTrue(r["old_48_source_sunlight_direction_heldout_not_replicated"])
        self.assertTrue(r["old_12_ERA5_raw_SHA_hold"])
        self.assertTrue(r["source_population_independence_not_verified"])


if __name__=="__main__":
    unittest.main()
