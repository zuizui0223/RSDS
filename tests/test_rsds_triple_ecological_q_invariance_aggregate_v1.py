"""Source study level tests for comparable q/night-day and q~daylength ecology."""
from __future__ import annotations
import unittest
from copy import deepcopy

from scripts.aggregate_rsds_triple_ecological_q_invariance_v1 import consistency
from scripts.run_rsds_triple_matched_ecological_q_v1 import (
    Q_KEYS,SLOPE_KEYS,freeze,NEGATIVE
)


def fake():
    _, first=freeze()
    src={}
    for old in first["all_original_source_studies_and_native_tag_triple_q"]:
        idx=old["original_study_index"]
        taxa=[]
        for x in old["taxa"]:
            n=(x.get("triple_original_animal_day_support") or {}).get(
                "source_original_tags_with_at_least_ten_triple_days",0)
            n_slope=(x.get("triple_original_animal_day_support") or {}).get(
                "source_original_tags_with_daylength_slope_triple_support",0)
            groups={}
            for name,count in [(k,n) for k in Q_KEYS]+[(k,n_slope) for k in SLOPE_KEYS]:
                groups[name]={
                    "status":"POSITIVE" if count>=5 else "UNAVAILABLE",
                    "n_source_tags":count,
                    "group_original_tag_mean":.1 if count>=5 else None,
                    "nominal_source_tag_bootstrap95":[.05,.15] if count>=5 else None
                }
            taxa.append({
                "taxon":x["taxon"],
                "source_status":"THREE_SCHEDULE_ECOLOGICAL_Q_ESTIMABLE" if n>=5 else "INSUFFICIENT_ORIGINAL_TRIPLE_TAG_DATE_SUPPORT",
                "original_three_schedule_source_support":{
                    "original_tags_with_at_least_ten_triple_qualified_days":n,
                    "matched_original_tagged_animal_calendar_days":(x.get("triple_original_animal_day_support") or {}).get("triply_qualified_original_animal_days",0),
                    "original_tags_with_common_photoperiod_slope_support":n_slope
                } if x.get("triple_original_animal_day_support") is not None else None,
                "q_original_tag_means_on_exact_same_dates":({k:groups[k] for k in Q_KEYS}
                    if x.get("triple_original_animal_day_support") is not None else None),
                "q_daylength_per_hour_on_exact_same_dates":({k:groups[k] for k in SLOPE_KEYS}
                    if x.get("triple_original_animal_day_support") is not None else None),
                "source_q_sign_recurrence":"THREE_SCHEDULES_NOMINAL_POSITIVE" if n>=5 else "UNAVAILABLE",
                "source_daylength_slope_sign_recurrence":"THREE_SCHEDULES_NOMINAL_POSITIVE" if n_slope>=5 else "UNAVAILABLE",
            })
        src[idx]={
            "version":"RSDS_8_STUDY_TRIPLE_DATE_ECOLOGICAL_Q_AND_DAYLENGTH_V1_SOURCE",
            "index":idx,"original_study_doi":old["original_study_doi"],
            "native_taxa":[x["taxon"] for x in old["taxa"]],
            "status":"PREDECLARED_FILE_METADATA_HOLD" if idx==6 else "SOURCE_TRIPLE_ECOLOGICAL_Q_COMPLETED",
            "original_provider_bitstream_MD5_confirmed":idx!=6,
            "original_source_q_and_daylength_replayed_before_new_q_levels":idx!=6,
            "previous_12_ERA5_raw_version_HOLDS_UNCHANGED":True,
            "previous_48_source_signed_shortwave_external_NONREPLICATED_UNCHANGED":NEGATIVE,
            "native_taxon_results":taxa
        }
    return src


class SourceDateMatchedTimeNicheStudyEcology(unittest.TestCase):
    def test_original_study_and_native_taxa_denominator_with_one_egret_hold(self):
        d=consistency(fake())
        self.assertEqual(d["original_source_studies"],8)
        self.assertEqual(d["original_study_taxon_rows"],12)
        self.assertEqual(d["study_status_counts"],{
            "SOURCE_TRIPLE_ECOLOGICAL_Q_COMPLETED":7,
            "PREDECLARED_FILE_METADATA_HOLD":1
        })
        self.assertEqual(d["native_taxa_q_three_schedule_nominal_sign_invariance"][
            "THREE_SCHEDULES_NOMINAL_POSITIVE"],8)

    def test_4_utah_native_taxa_are_only_one_original_published_field_study(self):
        r=consistency(fake())
        ut=r["original_study_DOI_equal_source_q_levels"]["10.5441/001/1.711"]
        self.assertEqual(ut["native_original_taxa_with_5_all_schedule_q_tags"],3)
        self.assertAlmostEqual(ut["source_study_equal_native_taxon_means_on_triple_dates"][
            "q_original"],.1)
        self.assertEqual(len(r["original_study_DOI_equal_source_q_levels"]),8)

    def test_missing_source_native_taxon_is_never_silently_dropped(self):
        d=fake()
        d[0]["native_taxon_results"].pop()
        with self.assertRaises(ValueError):
            consistency(d)
        d=fake()
        d.pop(6)
        with self.assertRaises(ValueError):
            consistency(d)

    def test_source_valid_effect_requires_original_source_MD5_and_q_replay(self):
        d=fake()
        d[1]["original_provider_bitstream_MD5_confirmed"]=False
        with self.assertRaises(ValueError):
            consistency(d)
        d=fake()
        d[1]["original_source_q_and_daylength_replayed_before_new_q_levels"]=False
        with self.assertRaises(ValueError):
            consistency(d)

    def test_under_five_tagged_original_animals_must_not_be_numeric_zero(self):
        d=fake()
        x=next(t for t in d[0]["native_taxon_results"] if t["taxon"]=="Antilocapra americana")
        self.assertEqual(x["q_original_tag_means_on_exact_same_dates"]["q_original"]["status"],"UNAVAILABLE")
        x["q_original_tag_means_on_exact_same_dates"]["q_original"]["group_original_tag_mean"]=0
        with self.assertRaises(ValueError):
            consistency(d)

    def test_original_egret_HOLD_and_old_shortwave_nonreplication_cannot_change(self):
        d=fake()
        d[6]["status"]="SOURCE_TRIPLE_ECOLOGICAL_Q_COMPLETED"
        with self.assertRaises(ValueError):
            consistency(d)
        d=fake()
        d[4]["previous_48_source_signed_shortwave_external_NONREPLICATED_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            consistency(d)

    def test_a_single_source_taxon_sign_flip_is_detected_as_not_invariant(self):
        d=fake()
        c=d[1]["native_taxon_results"][0]
        c["q_original_tag_means_on_exact_same_dates"]["q_odd"].update(
            status="NEGATIVE",group_original_tag_mean=-.1,
            nominal_source_tag_bootstrap95=[-.15,-.05])
        c["source_q_sign_recurrence"]="OPPOSITE_NOMINAL_SOURCE_TAG_SIGN"
        r=consistency(d)
        self.assertEqual(r["native_taxa_q_three_schedule_nominal_sign_invariance"]["OPPOSITE_NOMINAL_SOURCE_TAG_SIGN"],1)


if __name__=="__main__":
    unittest.main()
