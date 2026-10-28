"""Eight original field studies/twelve native taxa remain in strict triple-date aggregate."""
import unittest
from copy import deepcopy

from scripts.aggregate_rsds_three_phase_exact_original_gps_v2 import aggregate
from scripts.run_rsds_gps_three_phase_exact_dates_v2 import AXES,frozen,NEGATIVE


def placeholder():
    _, first=frozen()
    result={}
    for old in first["all_original_source_study_and_native_taxa"]:
        i=old["fixed_original_study_index"]
        taxa=[r["taxon"] for r in old["source_native_taxon_results"]]
        state="PREDECLARED_FILE_METADATA_HOLD" if i==6 else "ORIGINAL_SOURCE_OR_TRIPLE_MATCH_HOLD"
        result[i]={
            "version":"RSDS_EXACT_TRIPLE_MATCHED_ORIGINAL_GPS_STUDY_V2",
            "original_study_index":i,
            "original_study_doi":old["original_study_doi"],
            "original_native_taxa":taxa,
            "source_status":state,
            "original_source_provider_checksum_verified":False,
            "native_q_daylength_baseline_replayed_before_new_phase_comparison":False,
            "historical_12_ERA5_JSON_SOURCE_HOLDS_UNCHANGED":True,
            "historical_signed_shortwave_external_nonreplication_UNCHANGED":NEGATIVE,
            "native_taxon_results":[{
                "taxon":t,"source_status":"SOURCE_HOLD",
                "triple_source_effects":None,
            } for t in taxa]
        }
    return result


class TripleSourceEightOriginalStudyAggregator(unittest.TestCase):
    def test_eight_original_field_studies_and_nulls_not_zero(self):
        result=aggregate(placeholder())
        self.assertEqual(result["fixed_original_source_studies"],8)
        self.assertEqual(result["native_source_study_taxon_rows"],12)
        self.assertEqual(result["source_study_states"]["PREDECLARED_FILE_METADATA_HOLD"],1)
        self.assertEqual(len(result["all_original_source_studies_and_native_tag_triple_q"]),8)
        self.assertEqual(result["strict_three_phase_axis_nominal_tag_bootstrap_sign_counts"][
            "q_even_minus_odd"],{})
        self.assertTrue(result["old_12_ERA5_source_RAW_SHA_HOLDS_UNCHANGED"])

    def test_positive_triple_effect_rejected_if_original_md5_not_replayed(self):
        d=placeholder()
        x=d[1]
        x["source_status"]="TRIPLE_SOURCE_SUPPORT_EVALUATED"
        with self.assertRaises(ValueError):
            aggregate(d)

    def test_missing_fixed_study_taxon_or_old_signed_radiation_verdict_rejected(self):
        a=placeholder()
        del a[3]
        with self.assertRaises(ValueError):
            aggregate(a)
        a=placeholder()
        a[2]["historical_signed_shortwave_external_nonreplication_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            aggregate(a)
        a=placeholder()
        a[0]["native_taxon_results"].pop()
        with self.assertRaises(ValueError):
            aggregate(a)

    def test_triple_estimate_needs_five_tagged_animals_and_CI(self):
        p=placeholder()
        x=p[1]
        x["source_status"]="TRIPLE_SOURCE_SUPPORT_EVALUATED"
        x["original_source_provider_checksum_verified"]=True
        x["native_q_daylength_baseline_replayed_before_new_phase_comparison"]=True
        for t in x["native_taxon_results"]:
            t["source_status"]="TRIPLE_ORIGINAL_TAG_UTC_PARITY_ESTIMABLE"
            t["original_same_tag_date_triple_source_support"]={
                "source_original_tags_with_at_least_ten_triple_days":6,
                "source_original_tags_with_daylength_slope_triple_support":6,
                "triply_qualified_tag_years_with_at_least_ten_identical_dates":6,
                "triply_qualified_original_animal_days":180,
            }
            t["triple_source_effects"]={axis:{
                "source_original_same_tags_and_triple_dates":6,
                "group_source_mean_delta":0.02,
                "descriptive_source_tag_bootstrap95":[0.01,0.03],
                "status":"NOMINAL_POSITIVE"
            } for axis in AXES}
        result=aggregate(p)
        self.assertEqual(result["strict_three_phase_axis_nominal_tag_bootstrap_sign_counts"][
            "q_even_minus_odd"],{"NOMINAL_POSITIVE":2})
        p[1]["native_taxon_results"][0]["triple_source_effects"][
            "q_even_minus_odd"]["source_original_same_tags_and_triple_dates"]=4
        with self.assertRaises(ValueError):
            aggregate(p)

    def test_original_egret_historical_hold_cannot_be_reclassified(self):
        p=placeholder()
        p[6]["source_status"]="TRIPLE_SOURCE_SUPPORT_EVALUATED"
        with self.assertRaises(ValueError):
            aggregate(p)


if __name__=="__main__":
    unittest.main()
