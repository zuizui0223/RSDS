"""First actual original Illinois three-native shared-calendar source result gates."""
from __future__ import annotations
import hashlib
import json
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/"research/RSDS_ILLINOIS_THREE_NATIVE_SHARED_DATE_V1_FIRST_RESULT_20261029.json"
ORIGINAL_SHA="924156e8f98d0bb9fda2f9d23aabc03060af2a31"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def first():
    raw=FILE.read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
    if sha!=ORIGINAL_SHA:
        raise ValueError("Original Illinois native predator/prey source results were changed")
    return json.loads(raw.decode("utf-8"))


class OriginalIllinoisTripleNativeFirstResult(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=first()

    def test_original_MD5_prior_native_source_and_source_calendar_support(self):
        d=self.d
        self.assertEqual(d["status"],"THREE_NATIVE_SHARED_SOURCE_DATE_EVALUATED")
        self.assertEqual(d["study_DOI"],"10.5441/001/1.649")
        self.assertEqual(d["native_taxa"],["Canis latrans","Lynx rufus","Odocoileus virginianus"])
        self.assertTrue(d["MD5_verified"])
        self.assertTrue(d["first_all_native_q_daylength_and_solar_source_replayed"])
        cov=d["source_all_three_native_support"]
        self.assertEqual(cov["shared_source_calendar_dates_all_three_native_taxa"],169)
        self.assertEqual(cov["shared_source_year_month_blocks"],9)
        self.assertEqual(cov["native_original_tag_count_on_all_three_shared_source_dates"],
                         {"Canis latrans":8,"Lynx rufus":5,"Odocoileus virginianus":12})

    def test_three_source_native_time_q_order_significant_nominal_pairs(self):
        d=self.d
        native=d["source_all_native_q_levels_on_same_original_calendar"]
        expected={
            "Canis latrans":0.268521642425637,
            "Lynx rufus":0.09339268528468585,
            "Odocoileus virginianus":0.0024497395467369303
        }
        for sp,point in expected.items():
            self.assertAlmostEqual(native[sp]["q_original"]["source_animal_mean"],point,places=13)
        pairs=d["three_source_native_pair_q_and_daylength"]
        self.assertEqual(len(pairs),3)
        for key,rec in pairs.items():
            q=rec["q"]
            self.assertEqual(q["status"],"POSITIVE")
            self.assertTrue(q["original_tag_bootstrap95"][0]>0)
            self.assertTrue(q["first_minus_second_mean"]>0)
        self.assertAlmostEqual(pairs["Canis latrans MINUS Lynx rufus"]["q"]["first_minus_second_mean"],
                               .1751289571409511,places=13)

    def test_common_date_photoperiod_q_slopes_differ_by_native_species(self):
        d=self.d
        native=d["source_all_native_q_levels_on_same_original_calendar"]
        self.assertEqual(native["Canis latrans"]["daylength_original"]["nominal_sign"],"POSITIVE")
        self.assertEqual(native["Lynx rufus"]["daylength_original"]["nominal_sign"],"NEGATIVE")
        self.assertEqual(native["Odocoileus virginianus"]["daylength_original"]["nominal_sign"],"INCLUDES_ZERO")
        pairs=d["three_source_native_pair_q_and_daylength"]
        self.assertEqual(pairs["Canis latrans MINUS Lynx rufus"]["q_daylength"]["status"],"POSITIVE")
        self.assertEqual(pairs["Canis latrans MINUS Odocoileus virginianus"]["q_daylength"]["status"],"POSITIVE")
        self.assertEqual(pairs["Lynx rufus MINUS Odocoileus virginianus"]["q_daylength"]["status"],"NEGATIVE")
        for pair in pairs.values():
            self.assertEqual(pair["daily_mean_q_difference_vs_standard_daylength"]["source_calendar_days"],169)
            self.assertEqual(pair["daily_mean_q_difference_vs_standard_daylength"]["source_year_month_blocks"],9)

    def test_old_Washington_temporal_hold_and_external_shortwave_results_are_unchanged(self):
        d=self.d
        self.assertTrue(d["previous_washington_temporal_no_eligible_tags_HOLD_UNCHANGED"])
        self.assertTrue(d["old_12_ERA5_raw_SHA_HOLDS_UNCHANGED"])
        self.assertEqual(d["prior_48_negative_shortwave_direction_external_UNCHANGED"],NEGATIVE)

    def test_no_original_animal_or_exact_date_telemetry_republished(self):
        s=json.dumps(self.d)
        for forbidden in ('"animal_id"','"tag_id"','"individual_id"','"original_gps_fix"',
                          '"GPS_lat"','"GPS_lon"','"location-lat"','"location-long"',
                          '"timestamp"','"calendar_date"','"animal_year"'):
            self.assertNotIn(forbidden,s)


if __name__=="__main__":
    unittest.main()
