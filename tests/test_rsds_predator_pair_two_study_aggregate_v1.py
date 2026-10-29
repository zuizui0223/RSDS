"""Original 2-study/5-native-taxon predator temporal q source-complete known truth."""
from __future__ import annotations
import copy
import unittest

from scripts.aggregate_rsds_predator_pair_two_study_temporal_v1 import aggregate
from scripts.run_rsds_predator_pair_shared_date_temporal_v1 import frozen,NEGATIVE


def fake_sources():
    c,_=frozen()
    out={}
    for s in c["source_studies"]:
        out[s["region"]]={
            "version":"RSDS_COYOTE_BOBCAT_SHARED_DATE_TEMPORAL_ORIGINAL_STUDY_V1",
            "source_region":s["region"],
            "original_study_DOI":s["doi"],
            "native_taxa":s["native_taxa"],
            "status":"ORIGINAL_SOURCE_OR_SHARED_CALENDAR_HOLD",
            "provider_MD5_verified":False,
            "first_solar_source_native_hour_fix_counts_replayed":False,
            "original_native_source_full_temporal_q_and_daylength":[
                {"taxon":t,"q":None,"q_daylength":None} for t in s["native_taxa"]],
            "shared_original_calendar_comparison":None,
            "shared_triple_GPS_clock_comparison":None,
            "old_12_ERA5_RAW_SHA_HOLDS_UNCHANGED":True,
            "old_48_source_sun_direction_external_nonreplication_UNCHANGED":NEGATIVE
        }
    return out


class TwoCoyoteBobcatOriginalStudyTemporalAtlas(unittest.TestCase):
    def test_all_original_pinned_studies_and_five_native_taxa_with_ineligible_effects(self):
        d=aggregate(fake_sources())
        self.assertEqual(d["original_study_DOI_count"],2)
        self.assertEqual(d["native_study_taxon_rows"],5)
        self.assertEqual(set(d["source_original_outcome_states"]),{"Washington","Illinois"})
        self.assertEqual(d["two_study_native_pair_repeatability"]["ORIGINAL_ONLY"]["original"]
                         ["two_source_same_sign_classification"]["q"],
                         "INCONCLUSIVE_OR_SOURCE_HOLD")
        deer=d["all_original_source_study_native_taxa"][1]["source_native_taxa"]
        self.assertIn("Odocoileus virginianus",deer)
        self.assertIsNone(deer["Odocoileus virginianus"]["q"]["mean"])

    def test_missing_original_study_or_deer_is_not_allowed(self):
        a=fake_sources()
        del a["Washington"]
        with self.assertRaises(ValueError):
            aggregate(a)
        a=fake_sources()
        a["Illinois"]["original_native_source_full_temporal_q_and_daylength"].pop()
        with self.assertRaises(ValueError):
            aggregate(a)

    def test_two_original_study_nominal_native_pair_q_signs_without_new_globe_claim(self):
        r=fake_sources()
        def contrast(q,day="INCLUDES_ZERO"):
            return {"status":"BOTH_NATIVE_TAXA_SOURCE_COOCCURRENT_ESTIMABLE",
                    "source_calendar_dates_with_both_taxa":80,
                    "source_month_year_blocks":5,
                    "source_native_taxon_tag_support_on_shared_calendar":{
                        "Canis latrans":12,"Lynx rufus":8},
                    "native_pair":{
                        "coyote_minus_bobcat_q":{
                            "original":{"coyote_minus_bobcat_mean":.05 if q=="POSITIVE" else -.05,
                                        "nominal_sign":q}},
                        "coyote_minus_bobcat_daylength":{
                            "original":{"coyote_minus_bobcat_mean":.005,"nominal_sign":day}}
                    },
                    "date_level_pair_difference_vs_daylength":None}
        for x in r.values():
            x["status"]="ORIGINAL_SOURCE_SHARED_CALENDAR_EVALUATED"
            x["provider_MD5_verified"]=True
            x["first_solar_source_native_hour_fix_counts_replayed"]=True
            x["shared_original_calendar_comparison"]=contrast("POSITIVE")
        d=aggregate(r)
        self.assertEqual(d["two_study_native_pair_repeatability"]["ORIGINAL_ONLY"]["original"]
                         ["two_source_same_sign_classification"]["q"],
                         "SAME_NOMINAL_SIGN_TWO_SOURCE_STUDIES")
        r["Illinois"]["shared_original_calendar_comparison"]=contrast("NEGATIVE")
        d=aggregate(r)
        self.assertEqual(d["two_study_native_pair_repeatability"]["ORIGINAL_ONLY"]["original"]
                         ["two_source_same_sign_classification"]["q"],
                         "OPPOSITE_NOMINAL_SIGNS_TWO_STUDIES")

    def test_historical_external_radiation_failure_and_ERA5_hold_immutable(self):
        a=fake_sources()
        a["Illinois"]["old_48_source_sun_direction_external_nonreplication_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            aggregate(a)
        a=fake_sources()
        a["Washington"]["old_12_ERA5_RAW_SHA_HOLDS_UNCHANGED"]=False
        with self.assertRaises(ValueError):
            aggregate(a)


if __name__=="__main__":
    unittest.main()
