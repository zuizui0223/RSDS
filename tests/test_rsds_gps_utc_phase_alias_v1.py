"""Regression tests for source-conditioned GPS native UTC-hour phase alias."""
import json
import unittest
from pathlib import Path
from unittest.mock import patch

import scripts.summarize_rsds_gps_utc_phase_alias_v1 as model


class NativeOriginalGPSPhaseAliasingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=model.analyze()

    def test_frozen_eight_original_studies_and_twelve_native_taxon_rows(self):
        d=self.result
        self.assertEqual(d["original_field_studies"],8)
        self.assertEqual(d["source_taxon_rows"],12)
        self.assertEqual(d["source_cohorts_with_native_UTC_phase_fix_schedule_imbalance_below_20pct"],4)
        self.assertEqual(d["source_cohorts_with_near_balanced_UTC_phase_fix_schedule"],7)
        self.assertEqual(d["source_cohorts_metadata_HOLD"],1)
        self.assertEqual(len(d["all_original_source_taxon_phase_records"]),12)

    def test_utah_raw_hours_strongly_odds_not_every_other_hour_subsampling(self):
        rows=[x for x in self.result["all_original_source_taxon_phase_records"]
              if x["original_study_doi"]=="10.5441/001/1.711"]
        self.assertEqual(len(rows),4)
        self.assertTrue(all(x["source_original_fix_schedule_phase_class"]==
                            "PARITY_ALIASED_ORIGINAL_COLLAR_SCHEDULE" for x in rows))
        for x in rows:
            e,o=x["phase_effects"]
            self.assertLess(e["original_UTC_hour_fix_share"],0.2)
            self.assertGreater(o["original_UTC_hour_fix_share"],0.8)
            self.assertEqual(e["original_UTC_hour_selected_fixes"]+
                             o["original_UTC_hour_selected_fixes"],
                             x["total_selected_source_original_hourly_bin_fixes"])
            self.assertTrue(x["fixes_even_plus_odd_exact_partition"])

    def test_balanced_wyoming_two_species_show_opposite_nominal_phase_bias(self):
        opposite=self.result[
            "native_taxa_with_opposite_nominal_signed_q_cadence_effects_between_even_odd_UTC_phases"
        ]
        self.assertEqual({x["taxon"] for x in opposite},
                         {"Odocoileus hemionus","Ovis canadensis"})
        self.assertEqual({x["study_doi"] for x in opposite},
                         {"10.5441/001/1.703"})
        for row in self.result["all_original_source_taxon_phase_records"]:
            if row["original_study_doi"]=="10.5441/001/1.703":
                e,o=row["phase_effects"]
                self.assertEqual(e["q_nominal_animal_bootstrap_sign"],"POSITIVE")
                self.assertEqual(o["q_nominal_animal_bootstrap_sign"],"NEGATIVE")

    def test_no_fake_turkey_effect_or_great_egret_values(self):
        by={x["taxon"]:x for x in self.result["all_original_source_taxon_phase_records"]}
        turkey=by["Meleagris gallopavo"]
        self.assertEqual(turkey["phase_effects"][0]["matched_original_q_animal_count"],0)
        self.assertIsNone(turkey["phase_effects"][1]["q_two_hour_minus_onehour"])
        egret=by["Ardea alba"]
        self.assertEqual(egret["source_original_fix_schedule_phase_class"],"UNAVAILABLE")
        self.assertIsNone(egret["supported_q_effects"])
        self.assertEqual(self.result["source_cohorts_metadata_HOLD"],1)

    def test_one_original_study_gets_one_source_mean_not_four_ecosystems(self):
        by=self.result["source_original_study_equal_phase_q_bias"]
        self.assertEqual(set(by),set(model.PHASES))
        self.assertEqual(len(by["EVEN_UTC_HOURS"]),8)
        self.assertEqual(len(by["ODD_UTC_HOURS"]),8)
        e=next(x for x in by["EVEN_UTC_HOURS"] if x["original_study_doi"]=="10.5441/001/1.711")
        self.assertEqual(e["native_taxa_with_original_paired_q_effect"],4)
        self.assertAlmostEqual(e["native_taxon_study_mean_q_measurement_bias"],
                               0.0008635034457465542,places=12)

    def test_no_old_radiation_revision_or_individual_gps_export(self):
        d=self.result
        self.assertEqual(d["historic_48_source_sunlight_direction_external_NONREPLICATED_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")
        self.assertTrue(d["historical_12_ERA5_source_version_HOLDS_UNCHANGED"])
        s=json.dumps(d)
        for marker in ('"animal_id"','"tag_id"','"individual_id"',
                       '"timestamp"','"location-lat"','"location-long"',
                       '"latitude"','"longitude"','"original_fix"'):
            self.assertNotIn(marker,s)

    def test_previous_result_sha_change_fails(self):
        original=model.EXPECTED_SHA
        with patch.object(model,"EXPECTED_SHA","0"*40):
            with self.assertRaisesRegex(ValueError,"changed"):
                model.source()


if __name__=="__main__":
    unittest.main()
