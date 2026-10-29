"""Pin and audit four new original multi-native source results and goose study."""
from __future__ import annotations
import json
import math
import unittest
from unittest.mock import patch

import scripts.validate_rsds_four_new_multinative_temporal_first_v1 as audit


class OriginalMultinativeStudyResultsVerified(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=audit.infer()

    def test_original_four_field_studies_ten_taxa_with_no_source_replacement(self):
        d=self.data
        self.assertEqual(d["original_source_studies_selected_before_new_q"],4)
        self.assertEqual(d["registered_native_study_taxon_rows"],10)
        self.assertEqual(d["source_original_MD5_verified"],4)
        self.assertEqual(d["original_multinative_studies_with_any_all_taxa_shared_q_calendar"],1)
        self.assertEqual(d["original_studies_with_no_all_native_shared_calendar"],2)

    def test_new_goose_three_registered_species_exact_120_source_dates(self):
        d=self.data
        self.assertEqual(d["new_coobserved_geese_original_study_DOI"],"10.5441/001/1.fk899541")
        self.assertEqual(d["new_geese_original_same_source_calendar_days"],120)
        self.assertEqual(d["new_geese_source_year_month_blocks"],5)
        self.assertEqual({x["native_taxon"]:x["original_same_calendar_q_tag_count"] for x in
                          d["goose_native_q_and_study_reference_daylength"]},{
                              "Branta leucopsis":62,"Anser albifrons":11,"Anser anser":10})

    def test_three_original_goose_q_signs_and_opposing_temporal_slope(self):
        rows={r["native_taxon"]:r for r in self.data["goose_native_q_and_study_reference_daylength"]}
        self.assertTrue(all(r["q_nominal_interval_sign"]=="NEGATIVE" for r in rows.values()))
        self.assertEqual({k:r["daylength_nominal_interval_sign"] for k,r in rows.items()},{
            "Branta leucopsis":"NEGATIVE",
            "Anser albifrons":"INCLUDES_ZERO",
            "Anser anser":"POSITIVE"
        })
        self.assertAlmostEqual(rows["Branta leucopsis"]["observed_NIGHT_minus_DAY_GPS_q"],-.1840335077564767)
        self.assertAlmostEqual(rows["Anser anser"]["q_per_hour_study_reference_daylength_slope"],.020657883034855887)

    def test_all_three_native_goose_pairs_and_one_nominal_month_interval(self):
        d=self.data
        self.assertEqual(d["goose_native_pair_q_nominal_interval_counts"],{
            "POSITIVE":0,"NEGATIVE":2,"INCLUDES_ZERO":1
        })
        self.assertEqual(d["goose_within_month_pair_daylength_nominal_interval_counts"],{
            "POSITIVE":0,"NEGATIVE":1,"INCLUDES_ZERO":2
        })
        rows=d["all_three_frozen_original_goose_source_native_pairs"]
        self.assertEqual(len(rows),3)
        specific=next(x for x in rows if x["taxa"]==["Anser albifrons","Anser anser"])
        self.assertEqual(specific["source_year_month_blocks"],5)
        self.assertAlmostEqual(specific["within_month_study_reference_daylength_slope"],
                               -.03594158661342082,places=12)
        self.assertLess(specific["descriptive_source_month_block_bootstrap95"][1],0)

    def test_tortoise_schema_hold_mallard_source_q_and_vulture_no_support(self):
        d=self.data
        self.assertEqual(d["Galapagos_exact_registered_taxon_label_mapping"],
                        "UNRESOLVED_ZERO_EXACT_REGISTERED_TAXON_ROWS_NOT_BIOLOGICAL_ZERO")
        self.assertEqual(d["duck_mallard_full_study_q_eligible_tagged_animals"],21)
        self.assertEqual(d["duck_mallard_full_study_negative_q"],
                         "NEGATIVE_UNADJUSTED_ORIGINAL_TAG_BOOTSTRAP")
        self.assertEqual(d["African_vulture_source_night_day_q"],
                         "UNSUPPORTED_NO_APPROPRIATE_SOURCE_TAGS")

    def test_previous_external_direction_and_raw_ERA5_HOLDS_are_immutable(self):
        d=self.data
        self.assertTrue(d["previous_12_ERA5_raw_JSON_HOLDS_UNCHANGED"])
        self.assertEqual(d["previous_48_source_signed_shortwave_direction_EXTERNAL_NONREPLICATED_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")

    def test_prior_source_evidence_changes_fail_closed(self):
        with patch.object(audit,"SOURCE_BLOB","0"*40):
            with self.assertRaisesRegex(ValueError,"first GPS science blob changed"):
                audit.infer()
        with patch.object(audit,"METADATA_BLOB","0"*40):
            with self.assertRaisesRegex(ValueError,"source metadata or first GPS science blob changed"):
                audit.infer()

    def test_source_group_release_has_no_raw_original_tag_or_coordinates(self):
        serialized=json.dumps(self.data)
        for key in ('"animal_id"','"individual_id"','"tag_id"','"gps_fix"',
                    '"location-lat"','"location-long"','"latitude"','"longitude"',
                    '"calendar_date"','"event_timestamp"','"source_GPS_track"'):
            self.assertNotIn(key,serialized)


if __name__=="__main__":
    unittest.main()
