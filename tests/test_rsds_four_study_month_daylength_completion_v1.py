"""Unit-level four-field-source month ecology inventory and HOLD protection."""
from __future__ import annotations

import json
import unittest
from pathlib import Path

from scripts.summarize_rsds_four_original_shared_months_with_utah_v1 import (
    FOUR, FOUR_SHA, blob, combine_four, NEG, SOURCES
)


def fake_utah(status="UTAH_SOURCE_MONTH_CENTERED_ORIGINAL_EVALUATED"):
    def one(clock,mean,ci,sign):
        return {
            "source_GPS_UTC_clock_schedule":clock,
            "source_native_first":"Odocoileus hemionus",
            "source_native_second":"Ovis canadensis",
            "all_months_prior_source_slope_replay":"ARCHIVED_UNADJUSTED_SLOPE_REPLAYED",
            "all_months_original_shared_daylength_slope_per_hour":-.001,
            "month_centered_source_daylength_response":{
                "source_year_month_blocks":10,
                "eligible_source_month_blocks":9,
                "eligible_source_calendar_days":181,
                "within_month_daylength_slope":mean,
                "month_cluster_bootstrap95":ci,
                "nominal_interval_sign":sign,
            }
        }
    return {
        "version":"RSDS_UTAH_TWO_NATIVE_COOOBSERVED_MONTH_CENTERED_DAYLENGTH_V1_FIRST",
        "original_study_doi":"10.5441/001/1.711",
        "all_original_native_study_taxa":[
            "Antilocapra americana","Cervus elaphus",
            "Odocoileus hemionus","Ovis canadensis"
        ],
        "fixed_two_target_taxa":["Odocoileus hemionus","Ovis canadensis"],
        "status":status,
        "provider_original_GPS_bitstream_MD5_verified":True,
        "all_four_native_original_q_and_daylength_replayed_before_month_test":True,
        "original_two_taxa_triple_source_GPS_replayed_before_month_test":True,
        "strict_shared_source_date_support":{
            "source_shared_original_calendar_days":199,
            "original_source_year_month_blocks":10
        },
        "month_conditioned_original_pair_results":[
            one("original",.005,[-.03,.04],"INCLUDES_ZERO"),
            one("even",.02,[.01,.03],"POSITIVE"),
            one("odd",-.02,[-.03,-.01],"NEGATIVE")
        ],
        "earlier_four_study_Utah_within_month_not_tested_UNCHANGED":True,
        "old_external_48_source_negative_radiation_replication_UNCHANGED":NEG,
        "old_ERA5_raw_JSON_12_source_HOLDS_UNCHANGED":True
    }


class FourOriginalFieldStudyWithinMonthSourceAtlas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=FOUR.read_bytes()
        cls.first=json.loads(raw)
        cls.sha=blob(raw)
        if cls.sha!=FOUR_SHA:
            raise ValueError("Original first four-source result changed")

    def test_four_original_studies_eight_pairs_one_original_utah_new_result(self):
        z=combine_four(self.first,fake_utah(),self.sha,"x"*40)
        self.assertEqual(z["original_field_studies"],4)
        self.assertEqual(z["original_study_native_pair_rows"],8)
        self.assertEqual(z["old_original_three_studies_native_pair_rows"],7)
        self.assertEqual(z["new_Utah_month_source_pair_rows"],1)
        self.assertEqual(set(z["original_study_equal_sign_support"]),set(SOURCES))
        self.assertEqual(z["source_month_nominal_bootstrap_sign_counts"],{
            "POSITIVE":1,"NEGATIVE":2,"INCLUDES_ZERO":5
        })
        ut=next(x for x in z["all_four_studies_native_pair_month_results"]
                if x["original_source_field_study"]=="Utah")
        self.assertEqual(ut["original_native_pair"],["Odocoileus hemionus","Ovis canadensis"])
        self.assertEqual(ut["source_original_year_month_blocks"],10)
        self.assertEqual(ut["month_slope"],.005)

    def test_historic_utah_unestimated_state_not_reclassified_as_biological_zero(self):
        ut=fake_utah("ORIGINAL_SOURCE_OR_MONTH_SUPPORT_HOLD")
        ut["month_conditioned_original_pair_results"]=[]
        ut["strict_shared_source_date_support"]=None
        z=combine_four(self.first,ut,self.sha,"r"*40)
        entry=next(x for x in z["all_four_studies_native_pair_month_results"]
                   if x["original_source_field_study"]=="Utah")
        self.assertEqual(entry["nominal_source_month_sign"],"UNAVAILABLE")
        self.assertIsNone(entry["month_slope"])
        self.assertIsNone(entry["month_block_bootstrap95"])
        self.assertEqual(z["source_month_nominal_bootstrap_sign_counts"],{
            "POSITIVE":1,"NEGATIVE":2,"INCLUDES_ZERO":4,"UNAVAILABLE":1
        })

    def test_unverified_original_provider_or_daylight_slope_fails_closed(self):
        x=fake_utah()
        x["provider_original_GPS_bitstream_MD5_verified"]=False
        with self.assertRaises(ValueError):
            combine_four(self.first,x,self.sha,"x"*40)
        x=fake_utah()
        x["month_conditioned_original_pair_results"][0][
            "all_months_prior_source_slope_replay"]="ORIGINAL_NOT_REPLAYED"
        with self.assertRaises(ValueError):
            combine_four(self.first,x,self.sha,"x"*40)

    def test_prior_four_original_study_and_sunlight_negative_cannot_change(self):
        x=fake_utah()
        prior=dict(self.first)
        prior["old_48_source_negative_radiation_direction_external_replication_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            combine_four(prior,x,self.sha,"x"*40)
        with self.assertRaises(ValueError):
            combine_four(self.first,x,"0"*40,"x"*40)
        x["old_ERA5_raw_JSON_12_source_HOLDS_UNCHANGED"]=False
        with self.assertRaises(ValueError):
            combine_four(self.first,x,self.sha,"x"*40)

    def test_month_ci_sign_and_world_ecosystem_inference_protection(self):
        x=fake_utah()
        x["month_conditioned_original_pair_results"][0][
            "month_centered_source_daylength_response"]["month_cluster_bootstrap95"]=[.01,.04]
        with self.assertRaises(ValueError):
            combine_four(self.first,x,self.sha,"x"*40)
        z=combine_four(self.first,fake_utah(),self.sha,"x"*40)
        self.assertTrue(any("not eight independently" in q.lower()
                            for q in z["science_ceiling"]))
        self.assertNotIn('"tag_id"',json.dumps(z))
        self.assertNotIn('"latitude"',json.dumps(z))


if __name__=="__main__":
    unittest.main()
