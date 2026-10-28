"""Source-denominator and typed missingness tests for eight original GPS cadence studies."""
from __future__ import annotations
import unittest
from copy import deepcopy

from scripts.aggregate_rsds_eight_study_gps_cadence_v1 import aggregate
from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    verify_frame,PAIRS,NEGATIVE
)


def placeholder():
    old,new,reg,meta=verify_frame()
    a={
        0:("10.5441/001/1.711",[t["taxon"] for t in old["sources"][0]["taxon_results"]]),
        1:("10.5441/001/1.703",[t["taxon"] for t in old["sources"][1]["taxon_results"]])
    }
    for k,z in enumerate(reg,2):
        a[k]=(z["doi"],[z["taxon"]])
    rows={}
    for i,(doi,taxa) in a.items():
        rows[i]={
            "version":"RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_V1_SOURCE",
            "source_index":i,"original_study_doi":doi,
            "source_native_taxa":taxa,
            "source_status":"PREDECLARED_FILE_METADATA_HOLD",
            "original_GPS_provider_MD5_verified":False,
            "baseline_1h_q_and_daylength_source_effects_replayed":False,
            "native_taxon_results":[{"taxon":t,"source_status":"SOURCE_HOLD",
                "source_original_even_odd_2h_phase_sensitivity":None} for t in taxa],
            "unchanged_old_ERA5_12_RAW_SHA_HOLDS":True,
            "unchanged_48_source_signed_radiation_external_NONREPLICATION":NEGATIVE,
        }
    return rows


def observed_axes():
    def effect():
        return {
          "n_original_tagged_animals_paired_between_cadences":6,
          "status":"SOURCE_TAGGED_ANIMAL_DESCRIPTIVE",
          "paired_tag_mean_2h_minus_same_day_1h":-0.1,
          "original_tag_resample95":[-0.2,-0.05],
          "n_positive_tag_deltas":1,"n_negative_tag_deltas":5,
          "nominal_interval_sign":"NEGATIVE"
        }
    return {
        k:{
          "matched_q_velocity_proxy_difference":effect(),
          "matched_q_daylength_slope_difference_per_hour":effect(),
        } for k in PAIRS
    }


class EightStudyCadenceAggregateTests(unittest.TestCase):
    def test_all_native_studies_HOLD_without_zero_imputation(self):
        a=aggregate(placeholder())
        self.assertEqual(a["original_field_studies"],8)
        self.assertEqual(a["original_source_taxon_rows"],12)
        self.assertEqual(len(a["all_original_source_study_and_native_taxa"]),8)
        self.assertEqual(a["status_by_original_source_study"],{"PREDECLARED_FILE_METADATA_HOLD":8})
        self.assertEqual(a["original_q_paired_cadence_support_counts"]["EVEN_UTC_HOURS"]["matched_q"],{})
        self.assertTrue(all(x is None for rec in
                            a["one_original_DOI_one_study_phase_delta_mean"].values()
                            for x in rec.values()))

    def test_one_supported_native_taxon_requires_MD5_and_baseline_replay(self):
        x=placeholder()
        x[0]["source_status"]="SOURCE_CADENCE_SENSITIVITY_COMPLETED"
        x[0]["original_GPS_provider_MD5_verified"]=True
        x[0]["baseline_1h_q_and_daylength_source_effects_replayed"]=True
        for entry in x[0]["native_taxon_results"]:
            entry["source_original_even_odd_2h_phase_sensitivity"]=observed_axes()
            entry["source_original_full_hour_baseline_support"]={"tags_with_some_q_observable":7}
        a=aggregate(x)
        self.assertEqual(a["original_q_paired_cadence_support_counts"]["EVEN_UTC_HOURS"]["matched_q"],
                         {"NEGATIVE":4})
        self.assertEqual(a["original_q_paired_cadence_support_counts"]["ODD_UTC_HOURS"]["matched_daylength"],
                         {"NEGATIVE":4})
        x[0]["original_GPS_provider_MD5_verified"]=False
        with self.assertRaises(ValueError):
            aggregate(x)

    def test_unsupported_q_pair_with_fake_numeric_mean_fails(self):
        x=placeholder()
        x[2]["source_status"]="SOURCE_CADENCE_SENSITIVITY_COMPLETED"
        x[2]["original_GPS_provider_MD5_verified"]=True
        x[2]["baseline_1h_q_and_daylength_source_effects_replayed"]=True
        axes=observed_axes()
        axes["EVEN_UTC_HOURS"]["matched_q_velocity_proxy_difference"].update({
            "status":"INSUFFICIENT_MATCHED_TAGS",
            "n_original_tagged_animals_paired_between_cadences":3,
            "paired_tag_mean_2h_minus_same_day_1h":1.0
        })
        x[2]["native_taxon_results"][0].update(
            source_original_even_odd_2h_phase_sensitivity=axes,
            source_original_full_hour_baseline_support={"tags_with_some_q_observable":7})
        with self.assertRaises(ValueError):
            aggregate(x)

    def test_all_8_original_study_IDs_and_previous_48_source_failure_immutable(self):
        x=placeholder()
        x[7]["original_study_doi"]="10.5441/001/1.999"
        with self.assertRaises(ValueError):
            aggregate(x)
        x=placeholder()
        x[4]["unchanged_48_source_signed_radiation_external_NONREPLICATION"]="REPLICATED"
        with self.assertRaises(ValueError):
            aggregate(x)
        x=placeholder()
        del x[3]
        with self.assertRaises(ValueError):
            aggregate(x)


if __name__=="__main__":
    unittest.main()
