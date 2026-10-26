"""Synthetic source/HOLD tests for six newly fixed Movebank original study route."""
from __future__ import annotations
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from scripts.run_rsds_new_six_original_gps_ecology_v0 import source_analysis


SOURCE = {"id":"fake-package","doi":"10.5441/001/1.test-source",
          "study_id":"fake-original-study","region":"synthetic","guild":"test",
          "taxon":"Synthetic species","reported_animals":7,"reported_fixes":1000}


def admissible():
    return {"package_uuid":SOURCE["id"],"source_metadata_status":"ORIGINAL_FILE_METADATA_ADMISSIBLE",
            "original_file":{"bitstream_uuid":"original-gps","filename":"orig.csv",
                             "bytes":1000,"checksum":{"checkSumAlgorithm":"MD5","value":"a"*32}}}


class OriginalSixGPSKnownTruth(unittest.TestCase):
    def test_metadata_hold_does_not_download_original_gps(self):
        p={"source_metadata_status":"SOURCE_METADATA_HOLD",
           "metadata_hold_reason":"NO_UNAMBIGUOUS_SINGLE_PRIMARY_ORIGINAL_FILE"}
        with TemporaryDirectory() as d:
            with patch("scripts.run_rsds_new_six_original_gps_ecology_v0.original_files",
                       side_effect=AssertionError("Animal source opened without file qualification")):
                r=source_analysis(SOURCE,p,Path(d))
        self.assertEqual(r["source_status"],"ORIGINAL_METADATA_HOLD")
        self.assertIsNone(r["all_outcome_axes"]["sunward_raw"])
        self.assertIsNone(r["all_outcome_axes"]["same_tag_two_contexts"])
        self.assertNotIn('"animal_id"',json.dumps(r))

    def test_changed_original_md5_held_and_does_not_inspect_events(self):
        with TemporaryDirectory() as d:
            altered={"bitstream_uuid":"original-gps","bytes":1000,
                     "checksum":{"checkSumAlgorithm":"MD5","value":"b"*32}}
            with patch("scripts.run_rsds_new_six_original_gps_ecology_v0.original_files",
                       return_value=[altered]),patch(
                "scripts.run_rsds_new_six_original_gps_ecology_v0.load_source_taxa",
                side_effect=AssertionError("Must not read source when MD5 fails")
            ):
                r=source_analysis(SOURCE,admissible(),Path(d))
        self.assertEqual(r["source_status"],"ORIGINAL_GPS_SOURCE_HOLD")
        self.assertEqual(r["source_hold_error_type"],"StrictSourceHold")
        self.assertIsNone(r["all_outcome_axes"]["day_night_q"])

    def test_good_source_preserves_q_heading_and_same_animal_axes(self):
        q={"taxon":SOURCE["taxon"],"total_original_animals":7,"admitted_animals_diel":6,
           "admitted_animals_photoperiod_slope":5,
           "original_adjacent_hour_steps":1000,"qualified_daypart_step_observations":600,
           "day_night_speed_contrast":{"n_individuals":6,"mean":.05,
                                       "animal_bootstrap_95_percentile":[.01,.09]},
           "daylength_association":{"n_individuals":5,"mean":-.01,
                                    "animal_bootstrap_95_percentile":[-.04,.02]}}
        solar={"source_animals":7,"supported_individuals":5,
               "daylight_directional_qualified_steps":444,
               "status":"descriptive_source_estimated",
               "results":{"raw_toward_sun_mean_cos":{"n_individuals":5,"mean":.08,"bootstrap95":[.01,.14]},
                          "compass_independence_excess_cos":{"n_individuals":5,"mean":.002,
                             "bootstrap95":[-.01,.03]}}}
        geometry={"source_status":"PAIRED_TAG_UNAVAILABLE","paired_context_effect":None,
                  "original_q_tags_in_both_selected_contexts":0,
                  "selection_is_coverage_only":True}
        with TemporaryDirectory() as d:
            p=Path(d)
            with patch("scripts.run_rsds_new_six_original_gps_ecology_v0.original_files",
                       return_value=[{"bitstream_uuid":"original-gps","bytes":1000,
                       "checksum":{"checkSumAlgorithm":"MD5","value":"a"*32}}]),patch(
                 "scripts.run_rsds_new_six_original_gps_ecology_v0.checksum_and_download",
                 return_value={"bytes":1000,"checksum_type":"md5",
                               "provider_checksum_verified":True}),patch(
                 "scripts.run_rsds_new_six_original_gps_ecology_v0.load_source_taxa",
                 return_value=({SOURCE["taxon"]:{"private-tag":{}}},
                               {SOURCE["taxon"]:{"source_individuals":7}},{"original_total_rows":1000})),patch(
                 "scripts.run_rsds_new_six_original_gps_ecology_v0.one_taxon",
                 return_value=(q,{})),patch(
                 "scripts.run_rsds_new_six_original_gps_ecology_v0.solar_taxon",
                 return_value=solar),patch(
                 "scripts.run_rsds_new_six_original_gps_ecology_v0.paired_geometry_taxon",
                 return_value=geometry):
                r=source_analysis(SOURCE,admissible(),p)
        self.assertEqual(r["source_status"],"ORIGINAL_SOURCE_GPS_ECOLOGY_EVALUATED")
        self.assertEqual(r["all_outcome_axes"]["day_night_q"]["mean"],.05)
        self.assertEqual(r["all_outcome_axes"]["sunward_joint"]["mean"],.002)
        self.assertIsNone(r["all_outcome_axes"]["same_tag_two_contexts"]["paired_context_effect"])
        self.assertNotIn("private-tag",json.dumps(r))
        self.assertEqual(r["old_negative_external_shortwave_replication_UNCHANGED"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")


if __name__=="__main__":
    unittest.main()
