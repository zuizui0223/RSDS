"""Regression test of original Illinois+Wyoming shared-date within-month daylength comparisons."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.validate_rsds_shared_month_photoperiod_v1 import (
    PIN,ROOT,source_interpretation
)


class RSDSMonthSeasonTaxonEffectValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=source_interpretation()

    def test_all_two_study_six_registered_original_species_pairs_and_schedules(self):
        r=self.r
        self.assertEqual(r["original_source_field_studies"],2)
        self.assertEqual(r["source_native_taxon_and_schedule_comparisons"],6)
        self.assertEqual(r["all_original_unadjusted_native_daylength_slopes_replayed"],6)
        self.assertEqual(r["month_conditioned_comparisons_nominal_bootstrap_interval_signs"],{
            "POSITIVE":1,"NEGATIVE":1,"INCLUDES_ZERO":4
        })

    def test_original_illinois_169_days_and_wyoming_193_dates_unchanged(self):
        rows=self.r["source_case_results"]
        for x in rows:
            self.assertEqual(x["same_original_source_calendar_days"],
                             169 if x["region"]=="Illinois" else 193)
            self.assertEqual(x["original_source_month_blocks"],
                             9 if x["region"]=="Illinois" else 8)
            self.assertEqual(x["month_conditioned_year_month_blocks"],8)
        self.assertEqual({x["taxon_first"] for x in rows if x["region"]=="Illinois"},
                         {"Canis latrans","Lynx rufus"})

    def test_two_persistent_signs_are_weak_and_not_a_global_photoeffect(self):
        x={(i["region"],i["taxon_first"],i["taxon_second"],i["GPS_timestamp_schedule"]):i
           for i in self.r["source_case_results"]}
        ill=x[("Illinois","Lynx rufus","Odocoileus virginianus","original")]
        wyo=x[("Wyoming","Odocoileus hemionus","Ovis canadensis","original")]
        self.assertAlmostEqual(ill["month_conditioned_daylength_slope"],-0.0846572503992727,places=13)
        self.assertAlmostEqual(wyo["month_conditioned_daylength_slope"],0.010266184385105494,places=13)
        self.assertLess(ill["month_block_bootstrap95"][1],0)
        self.assertGreater(wyo["month_block_bootstrap95"][0],0)
        self.assertGreater(ill["month_block_bootstrap95"][1],-.002)
        self.assertLess(wyo["month_block_bootstrap95"][0],.001)

    def test_zero_spanning_month_intervals_are_not_imputed_zero_effects(self):
        rows=[x for x in self.r["source_case_results"]
              if x["month_conditioned_source_nominal_sign"]=="INCLUDES_ZERO"]
        self.assertEqual(len(rows),4)
        self.assertTrue(all(x["month_conditioned_daylength_slope"] != 0 for x in rows))
        self.assertTrue(all(a<=0<=b for a,b in
                            (x["month_block_bootstrap95"] for x in rows)))

    def test_source_bio_external_negative_and_ERA5_version_holds_unchanged(self):
        d=self.r
        self.assertTrue(d["old_12_ERA5_source_provider_raw_SHA_HOLDS_UNCHANGED"])
        self.assertEqual(d["old_external_48_source_negative_radiation_direction_nonreplication_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")
        self.assertEqual(len(d["original_source_first_input_git_blobs"]),3)

    def test_source_result_first_SHA_mismatch_halts_interpretation(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for src in PIN:
                p=root/src
                p.parent.mkdir(parents=True,exist_ok=True)
                p.write_bytes((ROOT/src).read_bytes())
            target=root/list(PIN)[2]
            d=json.loads(target.read_text())
            d["original_studies"]=3
            target.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"first original GPS source"):
                source_interpretation(root)

    def test_original_private_tag_dates_coords_not_exported(self):
        raw=json.dumps(self.r)
        for forbidden in ('"tag_id"','"individual_id"','"animal_id"',
                          '"timestamp"','"location-lat"','"location-long"',
                          '"event_year"','"event_month"','"event_date"',
                          '"latitude"','"longitude"','"raw_gps"'):
            self.assertNotIn(forbidden,raw)


if __name__=="__main__":
    unittest.main()
