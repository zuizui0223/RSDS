"""Immutable first scientific evidence proof for 35 original source GPS gates."""
from __future__ import annotations
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.rsds_437_systematic_source_frame_v1 import HERE,git_blob

CROSSWALK=HERE/"research/RSDS_35_ORIGINAL_SOURCE_GPS_OBSERVATION_EVIDENCE_CROSSWALK_V1.json"


def read():
    receipt=json.loads(CROSSWALK.read_text(encoding="utf-8"))
    if receipt["version"]!="RSDS_35_ORIGINAL_SOURCE_GPS_OBSERVATION_EVIDENCE_CROSSWALK_V1":
        raise ValueError("Original source first scientific receipt changed")
    evidence={}
    for name,info in receipt["original_sources"].items():
        p=HERE/info["path"]
        raw=p.read_bytes()
        if git_blob(raw)!=info["git_blob_sha"]:
            raise ValueError("Original 35-source GPS first scientific evidence blob altered: "+name)
        if p.suffix==".json":
            evidence[name]=json.loads(raw)
    return receipt,evidence


class Original35Evidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt,cls.records=read()

    def test_every_source_or_document_has_pinned_original_git_blob(self):
        self.assertEqual(len(self.receipt["original_sources"]),6)
        self.assertEqual([a for a in self.receipt["first_CI"] if "run" in a].__len__(),4)
        self.assertEqual(self.records["original_35_first_q"]["original_source_study_units"],35)
        self.assertEqual(self.records["original_GPS_funnel_first_result"]["original_source_funnel_replayed"],25)
        self.assertEqual(self.records["original_q_measurability_first_rank_result"]["GPS_parsed_original_sources"],25)

    def test_stagewise_original_GPS_admission_and_native_denominators(self):
        x=self.records["original_GPS_funnel_first_result"]
        y=self.receipt["verified"]
        self.assertEqual(x["original_public_Movebank_source_frame"],[437,259,120,35])
        self.assertEqual(len(x["individual_original_source_diagnostic_rows"]),35)
        self.assertEqual(len(x["original_native_aggregate_metric_rows"]),34)
        self.assertEqual(x["source_original_registered_native_rows"],50)
        z=x["original_per_fix_gate_census_NESTED_but_NOT_independent_samples"]
        self.assertEqual([z[k] for k in (
            "original_adjacent_source_fix_intervals",
            "accepted_original_0p5_to_2h_time_intervals",
            "accepted_same_solar_date_and_daypart_intervals",
            "accepted_finite_speed_max60kmh_intervals")],
            [1502093,1218362,1068433,1066632])
        self.assertEqual(y["parsed_source_studies"],25)
        self.assertEqual(y["first_original_q_complete_source_studies"],10)
        self.assertEqual(y["original_native_q_any_source_studies"],13)

    def test_source_equal_observation_means_and_posthoc_rank_AUC_not_global(self):
        a=self.records["original_GPS_funnel_first_result"][
          "source_equal_observation_gate_contrasts_by_archived_q_availability"]
        b=self.records["original_q_measurability_first_rank_result"][
          "original_measurement_scores"]
        self.assertEqual(a["complete_native_q"]["parsed_source_study_count"],10)
        self.assertEqual(a["not_complete_native_q"]["parsed_source_study_count"],15)
        self.assertAlmostEqual(
            a["complete_native_q"]["source_equal_median_total_accepted_over_adjacent"],
            .708111810119529,places=12)
        self.assertAlmostEqual(
            a["not_complete_native_q"]["source_equal_median_total_accepted_over_adjacent"],
            .7077159415602918,places=12)
        self.assertAlmostEqual(
            b["full_original_segment_acceptance_fraction"]["all_parsed_sources"][
              "AUC_HIGHSCORE_PREDICTS_MEASURABLE"],
            .5428571428571428,places=12)
        self.assertAlmostEqual(
            b["time_0p5_to_2h_fraction"]["all_parsed_sources"][
              "AUC_HIGHSCORE_PREDICTS_MEASURABLE"],
            .7071428571428572,places=12)
        self.assertTrue(self.records["original_q_measurability_first_rank_result"][
           "original_48_sunlight_external_negative_nonreplication_UNCHANGED"])
        self.assertTrue(self.receipt["old_ERA5_source_raw_SHA_12_HOLD"])

    def test_human_privacy_and_historical_sun_null_never_reclassified(self):
        x=self.records["original_GPS_funnel_first_result"]
        self.assertFalse(x["original_human_GPS_source_opened"])
        self.assertEqual(x["original_source_diagnostic_status_counts"],{
            "REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL":25,
            "ARCHIVED_ORIGINAL_SOURCE_HOLD_UNCHANGED":9,
            "HUMAN_LOCATION_PRIVACY_HOLD_NO_ACCESS":1
        })
        self.assertEqual(self.receipt["old_48_external_solar_negative"],
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE")

    def test_changing_first_frozen_Git_blob_must_stop(self):
        copied=json.loads(CROSSWALK.read_text(encoding="utf-8"))
        copied["original_sources"]["original_GPS_funnel_first_result"][
            "git_blob_sha"]="0"*40
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"falsified.json"
            path.write_text(json.dumps(copied),encoding="utf-8")
            with patch(__name__+".CROSSWALK",path):
                with self.assertRaisesRegex(ValueError,"scientific evidence blob altered"):
                    read()

if __name__=="__main__":
    unittest.main()
