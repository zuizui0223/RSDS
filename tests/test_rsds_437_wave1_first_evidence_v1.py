"""Immutable evidence crosswalk test for first 437-source selection + 16 GPS wave."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from scripts.rsds_437_systematic_source_frame_v1 import HERE,git_blob

X=HERE/"research/RSDS_437_WAVE1_METHOD_EVIDENCE_CROSSWALK_V1.json"

def load_crosswalk():
    d=json.loads(X.read_text(encoding="utf-8"))
    if d["version"]!="RSDS_437_WAVE1_EVIDENCE_CROSSWALK_V1":
        raise ValueError("Unfrozen first original event evidence crosswalk")
    data={}
    for key,info in d["source_blobs"].items():
        raw=(HERE/info["path"]).read_bytes()
        if git_blob(raw)!=info["git_blob"]:
            raise ValueError("First study evidence blob modified: "+key)
        if info["path"].endswith(".json"):
            data[key]=json.loads(raw)
    return d,data

class FirstSourceEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cross,cls.data=load_crosswalk()

    def test_all_source_blobs_identity_matched(self):
        self.assertEqual(len(self.cross["source_blobs"]),9)
        self.assertEqual(self.data["original_437_catalog"]["package_count"],437)
        self.assertEqual(len(self.data["original_437_catalog"]["packages"]),437)
        self.assertEqual(self.data["wave_selected_sources"]["animal_event_records_accessed"],0)
        self.assertEqual(self.data["wave_original_file_metadata"]["original_GPS_event_records_opened"],0)
        self.assertEqual(self.cross["ordered_execution_receipts"]["native_GPS_all_16_source_first_successful_workflow"],43690283579)

    def test_first_16_source_record_exact_original_numerical_denominators(self):
        d=self.data["wave_first_native_GPS_effects"]
        n=self.cross["exact_denominators"]
        self.assertEqual(d["original_catalog_package_count"],n["catalog_packages"])
        self.assertEqual(d["original_GPS_GNSS_candidates"],n["gps_sensor_packages"])
        self.assertEqual(d["metadata_eligible_finite_frame"],n["stage1_eligible_metadata"])
        self.assertEqual(d["original_selection_wave_slots"],n["wave_total"])
        self.assertEqual(d["original_source_metadata_file_HOLD_slots"],n["original_file_metadata_HOLD"])
        self.assertEqual(len(d["all_first_wave_source_study_summaries"]),16)
        self.assertEqual(len(d["all_first_wave_registered_native_group_effects"]),27)
        self.assertEqual(
            d["source_status_counts"],
            {"ORIGINAL_SOURCE_GPS_HOLD":2,"ORIGINAL_SOURCE_GPS_Q_ESTIMABLE":6,
             "ORIGINAL_SOURCE_GPS_Q_UNAVAILABLE_NOT_ZERO":5,
             "PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD":3}
        )
        self.assertEqual(
            d["stratified_source_complete_effect_sensitivity"]["q"]["first_wave_source_complete_response"],4)
        self.assertEqual(
            d["stratified_source_complete_effect_sensitivity"]["daylength"]["first_wave_source_complete_response"],3)
        self.assertEqual(sum(x["native_covariance_status"]=="SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED"
                             for x in d["all_first_wave_source_study_summaries"]),0)

    def test_first_16_study_source_averages_and_nonresponse_sensitivity_not_global(self):
        d=self.data["wave_first_native_GPS_effects"]
        v=self.cross["descriptive_observational_results"]
        studies=d["all_first_wave_source_study_summaries"]
        supported=[s["source_q_equal_native_if_complete"] for s in studies
                   if s["source_q_equal_native_if_complete"] is not None]
        self.assertEqual(len(supported),4)
        for x,y in zip(sorted(supported),sorted(v["source_full_q_examples"].values())):
            self.assertAlmostEqual(x,y,places=12)
        sample=d["stratified_source_complete_effect_sensitivity"]["q"]
        self.assertAlmostEqual(
            sample["response_only_stratum_poststratified_point_DESCRIPTIVE_NOT_DESIGN_UNBIASED"],
            v["respondent_only_source_response_weighted_q"],places=12)
        self.assertEqual(
            sample["hypothetical_fill_missing_selected_sources_with_observed_range"],
            v["hypothetical_missing_q_range_under_observed_complete_source_range"])
        self.assertTrue(sample["no_qualified_world_population_inference"])
        self.assertTrue(d["not_a_world_population_meta_analysis"])

    def test_frozen_prior_nulls_and_no_raw_private_data(self):
        d=self.data["wave_first_native_GPS_effects"]
        for s in d["all_first_wave_registered_native_group_effects"]:
            for k in ("q","daylength"):
                if s[k] is not None:
                    self.assertGreaterEqual(s[k]["n_individuals"],5)
                    self.assertEqual(len(s[k]["animal_bootstrap_95_percentile"]),2)
            self.assertFalse(any(x in s for x in (
                "animal_id","tag_id","event_date","timestamp","latitude","longitude")))
        self.assertEqual(
            self.cross["historical_external_48_signed_radiation_nonreplication"],
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")
        self.assertTrue(self.cross["historical_12_original_ERA5_raw_HOLD"])
        self.assertTrue(self.data["wave_selected_sources"]["old_ERA5_raw_SHA_HOLD"])

    def test_modified_original_first_result_fails_closed(self):
        previous=self.cross["source_blobs"]["wave_first_native_GPS_effects"]["git_blob"]
        self.cross["source_blobs"]["wave_first_native_GPS_effects"]["git_blob"]="0000000000000000000000000000000000000000"
        try:
            # check the same data using a patched crosswalk file in place
            text=json.dumps(self.cross)
            from tempfile import TemporaryDirectory
            with TemporaryDirectory() as temp:
                candidate=Path(temp)/"crosswalk.json"
                candidate.write_text(text,encoding="utf-8")
                with patch(__name__+".X",candidate):
                    with self.assertRaisesRegex(ValueError,"First study evidence blob modified"):
                        load_crosswalk()
        finally:
            self.cross["source_blobs"]["wave_first_native_GPS_effects"]["git_blob"]=previous

if __name__=="__main__":
    unittest.main()
