"""First evidence crosswalk for 437 -> 120 -> 35 original source meta."""
from __future__ import annotations
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from scripts.rsds_437_systematic_source_frame_v1 import HERE,git_blob

CROSS=HERE/"research/RSDS_WAVE2_35_SOURCE_META_EVIDENCE_CROSSWALK_V1.json"

def load():
    m=json.loads(CROSS.read_text(encoding="utf-8"))
    if m["version"]!="RSDS_WAVE2_35_SOURCE_ORIGINAL_META_METHOD_EVIDENCE_CROSSWALK_V1":
        raise ValueError("Unknown original source first evidence receipt")
    evidence={}
    for key,ref in m["source_blobs"].items():
        raw=(HERE/ref["path"]).read_bytes()
        if git_blob(raw)!=ref["git_blob_sha"]:
            raise ValueError("Frozen first source evidence altered: "+key)
        if ref["path"].endswith(".json"):
            evidence[key]=json.loads(raw)
    return m,evidence

class OriginalWave2Evidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest,cls.source=load()

    def test_all_first_blobs_replay_and_access_phases(self):
        self.assertEqual(len(self.manifest["source_blobs"]),8)
        self.assertEqual([x["run_id"] for x in self.manifest["original_ci_order"]],
                         [43691550739,43691736136,43692435152,
                          43693018390,43693286399,43693427256])
        self.assertEqual(
            self.source["wave2_19_SHA_source_selection"]["original_wave2_GPS_events_opened"],0)
        self.assertEqual(
            self.source["wave2_19_original_file_first_metadata"]["all_new_wave2_animal_GPS_events_opened"],0)
        self.assertFalse(self.source["wave2_19_original_file_first_metadata"][
            "human_GPS_or_position_data_retrieved"])

    def test_exact_35_source_and_original_registered_native_rows(self):
        d=self.source["original_first_35_source_meta"]
        s=self.manifest["exact_denominators"]
        self.assertEqual(d["original_source_study_units"],s["total_source_studies"])
        self.assertEqual(d["original_metadata_candidate_frame"],120)
        self.assertEqual(d["combined_original_registered_native_rows"],50)
        self.assertEqual(d["wave2_source_study_units"],19)
        self.assertEqual(d["wave2_status_counts"],{
            "HUMAN_LOCATION_PRIVACY_HOLD":1,
            "SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE":7,
            "SOURCE_FILE_METADATA_HOLD":3,
            "SOURCE_GPS_HOLD":1,
            "SOURCE_NO_NATIVE_Q_ESTIMABLE":7
        })
        self.assertEqual(len(d["combined_source_study_summary_rows"]),35)
        self.assertEqual(d["finite_120_metadata_frame_source_measurability"]["q"][
            "by_source_stratum"],{"single_native":8,"multi_native":2})

    def test_source_level_q_not_world_claim(self):
        m=self.source["study_source_sign_fragility_first"]
        a=self.manifest["results"]
        x=m["response_only_poststratified_q"]
        self.assertAlmostEqual(x["wave1_original_16_source_response_only"],
                               a["first_wave_respondent_only_stratum_adjusted_q"],places=12)
        self.assertAlmostEqual(x["wave1_plus_wave2_35_source_response_only"],
                               a["combined_35_respondent_only_stratum_adjusted_q"],places=12)
        self.assertEqual(x["leave_one_source_out_range"],
                         a["single_DOI_leave_out_respondent_mean_range"])
        self.assertEqual(x["leave_two_source_out_range"],
                         a["two_DOI_leave_out_respondent_mean_range"])
        self.assertEqual(x["hypothetical_missing_effect_sensitivity_for_35_selected"],
                         a["source_q_missing_effect_conditional_observed_range"])
        self.assertTrue(x["not_population_estimand_not_confirmatory"])
        self.assertEqual(m["source_effect_signs"],{
            "respondent_full_source_q_negative":3,
            "respondent_full_source_q_positive":7,
            "respondent_full_source_q_zero":0})

    def test_privacy_holds_and_old_signed_irradiance_failure(self):
        m=self.manifest
        self.assertTrue(m["historic_ERA5_original_raw_SHA_12_HOLD"])
        self.assertEqual(m["historic_external_48_signed_shortwave"],
                         "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")
        a=self.source["original_first_35_source_meta"]
        self.assertTrue(a["human_location_source_no_access"])
        self.assertTrue(a["previous_48_source_sunlight_signed_external_not_replicated"])
        self.assertTrue(a["previous_12_ERA5_raw_SHA_HOLD"])
        self.assertTrue(a["no_original_GPS_tracks_IDs_exact_dates_coordinates_or_bootstrap_draws_exported"])

    def test_changed_original_source_Git_blob_fails_closed(self):
        old=self.manifest["source_blobs"]["original_first_35_source_meta"]["git_blob_sha"]
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            candidate=Path(d)/"wrong-crosswalk.json"
            changed=json.loads(CROSS.read_text())
            changed["source_blobs"]["original_first_35_source_meta"]["git_blob_sha"]="0"*40
            candidate.write_text(json.dumps(changed),encoding="utf-8")
            with patch(__name__+".CROSS",candidate):
                with self.assertRaisesRegex(ValueError,"Frozen first source evidence altered"):
                    load()
        self.assertEqual(self.manifest["source_blobs"]["original_first_35_source_meta"]["git_blob_sha"],old)

if __name__=="__main__":
    unittest.main()
