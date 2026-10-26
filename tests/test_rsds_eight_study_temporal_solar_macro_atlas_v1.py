"""Frozen eight study source macro ecological denominator and sign support."""
import json
import tempfile
import unittest
from pathlib import Path

from scripts.build_rsds_eight_study_temporal_solar_macro_atlas_v1 import (
    ROOT, SOURCES, NEGATIVE, macro, observed
)


class OriginalEightStudyMacroEcologyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.o = macro()

    def test_eight_original_studies_twelve_cohorts_one_egret_HOLD(self):
        d=self.o
        self.assertEqual(d["original_study_denominator"],8)
        self.assertEqual(d["source_taxon_rows"],12)
        self.assertEqual(d["original_studies_with_admissible_q"],7)
        self.assertEqual(d["source_taxa_with_q"],11)
        self.assertEqual(d["source_taxa_with_solar_raw"],10)
        self.assertEqual(d["source_taxa_with_paired_original_tags"],9)
        self.assertEqual(d["source_tag_context_pair_occurrences_NOT_global_biomes"],288)
        egret=next(x for x in d["all_sources_including_egret_original_metadata_HOLD"]
                   if x["taxon"]=="Ardea alba")
        self.assertEqual(egret["source_status"],"ORIGINAL_METADATA_HOLD")
        self.assertEqual(egret["q"]["status"],"UNAVAILABLE")

    def test_original_night_day_q_signs_not_all_one_direction(self):
        d=self.o["observed_nominal_source_animal_bootstrap_sign_counts"]
        self.assertEqual(d["night_day_q"],{
            "NEGATIVE_NOMINAL_SOURCE_TAG_INTERVAL":7,
            "POSITIVE_NOMINAL_SOURCE_TAG_INTERVAL":4,
            "UNAVAILABLE":1
        })
        self.assertEqual(d["q_association_with_daylength_per_hour"],{
            "POSITIVE_NOMINAL_SOURCE_TAG_INTERVAL":5,
            "NEGATIVE_NOMINAL_SOURCE_TAG_INTERVAL":1,
            "INCLUDES_ZERO":5,
            "UNAVAILABLE":1
        })

    def test_original_heading_raw_separate_from_joint_uncertainty(self):
        d=self.o["observed_nominal_source_animal_bootstrap_sign_counts"]
        self.assertEqual(d["mean_sun_relative_heading_cos"],{
            "POSITIVE_NOMINAL_SOURCE_TAG_INTERVAL":2,
            "NEGATIVE_NOMINAL_SOURCE_TAG_INTERVAL":2,
            "INCLUDES_ZERO":6,
            "UNAVAILABLE":2
        })
        self.assertEqual(d["joint_heading_sun_covariance"],{
            "POINT_ONLY_NO_ORIGINAL_CI":6,
            "INCLUDES_ZERO":4,
            "UNAVAILABLE":2
        })

    def test_same_source_taxa_rows_are_not_eight_independent_ecosystems(self):
        d=self.o["all_sources_including_egret_original_metadata_HOLD"]
        self.assertEqual(len({x["source_study_doi"] for x in d}),8)
        self.assertTrue(self.o["source_equal_unit_not_12_independent_ecosystems"])
        self.assertEqual(len({x["taxon"] for x in d}),10)
        self.assertEqual(self.o["old_48_radiation_study_external_negative_direction_test_UNCHANGED"],
                         NEGATIVE)

    def test_point_only_old_joint_covariance_must_not_get_false_CI(self):
        val=observed(.01,None,15,point_only=True)
        self.assertEqual(val["status"],"POINT_ONLY_NO_ORIGINAL_CI")
        with self.assertRaises(ValueError):
            observed(.01,[.002,.02],15,point_only=True)
        self.assertEqual(observed(None,None,0)["status"],"UNAVAILABLE")

    def test_changed_previous_original_source_result_fails_sha_integrity(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for path in SOURCES:
                target=root/path
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes((ROOT/path).read_bytes())
            target=root/list(SOURCES)[1]
            payload=json.loads(target.read_text(encoding="utf-8"))
            payload["frozen_source_original_study_count"]=7
            target.write_text(json.dumps(payload),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"changed unexpectedly"):
                macro(root)

    def test_no_gps_raw_individual_identifiers_coordinates_in_macro(self):
        x=json.dumps(self.o)
        for prohibited in ('"individual_id"', '"tag_ids"', '"location-lat"',
                           '"location-long"', '"UTC_timestamp"','"raw_gps_fixes"'):
            self.assertNotIn(prohibited,x)


if __name__=="__main__":
    unittest.main()
