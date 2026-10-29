"""Source study and scientific status preservation for photoperiod within-month contrasts."""
from __future__ import annotations
import json
import unittest
from copy import deepcopy

from scripts.aggregate_rsds_within_month_shared_calendar_v1 import aggregate
from scripts.run_rsds_shared_calendar_month_photoperiod_v1 import frozen


def fake_source(region):
    _,ref,first=frozen(region)
    a={
        "version":"RSDS_TEST_ONLY","region":region,
        "source_DOI":ref["doi"],
        "native_taxa":ref["taxa"],
        "status":"ORIGINAL_IL_THREE_NATIVE_MONTH_CONTRAST_COMPLETE" if region=="Illinois" else
                 "ORIGINAL_WY_TWO_NATIVE_MONTH_CONTRAST_COMPLETE",
        "provider_original_file_MD5_replayed":True,
        "prior_native_q_and_daylength_replayed_before_month_diagnosis":True,
        "same_original_calendar_days":ref["shared_native_calendar_days"],
        "source_year_month_blocks":ref["year_month_blocks"],
        "source_native_pair_month_contrasts":[],
        "historic_12_ERA5_source_HOLDS_UNCHANGED":True,
        "historic_48_external_sunlight_signed_nonreplication_UNCHANGED":
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
    }
    for pair in ref["pairs"]:
        for schedule in ref["schedules"]:
            old=first["three_source_native_pair_q_and_daylength"][f"{pair[0]} MINUS {pair[1]}"]["daily_mean_q_difference_vs_standard_daylength"] \
                if region=="Illinois" else first["calendar_date_paired_source_difference"]["slopes"][schedule]
            x=old["slope"] if region=="Illinois" else old["mule_q_minus_bighorn_q_vs_shared_reference_daylength"]
            a["source_native_pair_month_contrasts"].append({
                "native_first":pair[0],"native_second":pair[1],
                "source_GPS_clock_schedule":schedule,
                "same_study_original_calendar_days":ref["shared_native_calendar_days"],
                "previous_all_months_source_daylength_contrast_coefficient":x,
                "previous_unadjusted_source_nominal_sign":old["status"],
                "archived_first_unadjusted_original_coefficient_replayed":"ARCHIVED_UNADJUSTED_SLOPE_REPLAYED",
                "within_calendar_month_analysis":{
                    "source_calendar_dates":ref["shared_native_calendar_days"],
                    "source_year_month_blocks":ref["year_month_blocks"],
                    "eligible_source_month_blocks":6,
                    "eligible_source_calendar_days":70,
                    "within_month_daylength_centered_sum_squares_hours2":3.,
                    "status":"MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE",
                    "within_month_daylength_slope":.12,
                    "month_cluster_bootstrap95":[.05,.19],
                    "leave_one_month_block_out_range":[.1,.2],
                    "bootstrap_finite_replicates":2000,
                    "nominal_interval_sign":"POSITIVE"
                },
                "between_vs_within_month_status":"PREVIOUS_POSITIVE_ASSOCIATION_PERSISTS_WITHIN_MONTH_EXPLORATORY"
            })
    return a


class SourceStudyMonthResultAggregation(unittest.TestCase):
    def test_all_two_original_studies_six_fixed_native_source_pair_schedules(self):
        x=aggregate({"Illinois":fake_source("Illinois"),
                     "Wyoming":fake_source("Wyoming")})
        self.assertEqual(x["original_studies"],2)
        self.assertEqual(x["original_study_native_taxa"],{"Illinois":3,"Wyoming":2})
        self.assertEqual(x["original_source_pair_clock_comparisons_expected"],6)
        self.assertEqual(x["original_source_pair_clock_comparisons_reported"],6)
        self.assertEqual(x["contrasts_with_estimable_month_adjusted_daylength_slope"],6)
        self.assertEqual(x["month_conditional_nominal_bootstrap_sign_counts"],{"POSITIVE":6})

    def test_six_fixed_contrasts_cannot_become_effect_sign_selected_subset(self):
        s=fake_source("Illinois")
        s["source_native_pair_month_contrasts"].pop()
        with self.assertRaises(ValueError):
            aggregate({"Illinois":s,"Wyoming":fake_source("Wyoming")})

    def test_source_with_weak_month_support_keeps_unavailable_not_zero(self):
        s=fake_source("Illinois")
        x=s["source_native_pair_month_contrasts"][0]["within_calendar_month_analysis"]
        x.update(status="INSUFFICIENT_COMMON_SOURCE_MONTH_SUPPORT",
                 within_month_daylength_slope=None,month_cluster_bootstrap95=None,
                 nominal_interval_sign="UNAVAILABLE")
        d=aggregate({"Illinois":s,"Wyoming":fake_source("Wyoming")})
        self.assertEqual(d["contrasts_with_estimable_month_adjusted_daylength_slope"],5)
        self.assertEqual(d["month_conditional_nominal_bootstrap_sign_counts"]["UNAVAILABLE"],1)
        x["within_month_daylength_slope"]=0.0
        with self.assertRaises(ValueError):
            aggregate({"Illinois":s,"Wyoming":fake_source("Wyoming")})

    def test_prior_source_DOI_MD5_q_daylight_replay_required(self):
        s=fake_source("Wyoming")
        s["prior_native_q_and_daylength_replayed_before_month_diagnosis"]=False
        with self.assertRaises(ValueError):
            aggregate({"Illinois":fake_source("Illinois"),"Wyoming":s})

    def test_prior_external_negative_radiation_is_immutable(self):
        s=fake_source("Illinois")
        s["historic_48_external_sunlight_signed_nonreplication_UNCHANGED"]="REPLICATED"
        with self.assertRaises(ValueError):
            aggregate({"Illinois":s,"Wyoming":fake_source("Wyoming")})


if __name__=="__main__":
    unittest.main()
