"""Quantify observation-stage attrition and source-level meta sign fragility.

Reads ONLY frozen first source-level aggregate receipts (no GPS events).
Every scientific source-level HOLD, the human privacy exclusion and the
original 16+19 study identities are preserved. Leave-one/two source results
are POST-OUTCOME descriptive sensitivities, not valid refits or p-values.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
import math
from pathlib import Path
from statistics import mean

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned

FIRST=("research/RSDS_437_WAVE1_PLUS_WAVE2_35_ORIGINAL_GPS_FIRST_META_RESULT_20261115.json",
       "9b02774f7c31d29b2eec5039a80431a895dccb36")
EARLIER=("research/RSDS_437_WAVE1_16_SOURCE_GPS_FIRST_ALL_RESULTS_20261115.json",
         "57029f30b65991118325c7d0cec1d8b3a36f4427")
BASE_WEIGHTS={"single_native":109/120,"multi_native":11/120}
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"

def weighted_effect(rows):
    sampled={k:[float(r["source_q_mean_equal_registered_native"])
                for r in rows if r["source_stratum"]==k
                and r["source_q_mean_equal_registered_native"] is not None]
             for k in BASE_WEIGHTS}
    if not all(sampled.values()):
        return None
    return sum(BASE_WEIGHTS[k]*mean(sampled[k]) for k in BASE_WEIGHTS)

def run():
    current=pinned(*FIRST)
    earlier=pinned(*EARLIER)
    rows=current["combined_source_study_summary_rows"]
    if (current["original_source_study_units"]!=35 or len(rows)!=35
        or current["wave1_source_study_units"]!=16
        or current["wave2_source_study_units"]!=19
        or current["combined_original_registered_native_rows"]!=50
        or current["original_metadata_candidate_frame"]!=120
        or current["human_location_source_no_access"] is not True
        or current["previous_48_source_sunlight_signed_external_not_replicated"] is not True
        or current["previous_12_ERA5_raw_SHA_HOLD"] is not True
        or len(earlier["all_first_wave_source_study_summaries"])!=16
        or len({r["source_original_repository_study_id"] for r in rows})!=35):
        raise ValueError("Frozen 16+19 original GPS source/old biological receipts changed")
    earlier_q=earlier["stratified_source_complete_effect_sensitivity"]["q"]
    present=current["source_complete_effects_and_nonresponse_sensitivity"]["q"]
    qualified=[r for r in rows if r["source_q_mean_equal_registered_native"] is not None]
    if (len(qualified)!=10 or present["supported_complete_native_source_count"]!=10
        or len([r for r in qualified if r["source_stratum"]=="single_native"])!=8
        or len([r for r in qualified if r["source_stratum"]=="multi_native"])!=2):
        raise ValueError("Original supported source q strata or number changed")
    source_actual=weighted_effect(qualified)
    frozen=present["respondent_only_poststratified_mean_EXPLORATORY"]
    if (not math.isclose(source_actual,frozen,abs_tol=1e-12,rel_tol=0)
        or not math.isclose(earlier_q[
            "response_only_stratum_poststratified_point_DESCRIPTIVE_NOT_DESIGN_UNBIASED"
            ],-.0477935014187846,abs_tol=1e-12,rel_tol=0)):
        raise ValueError("Wave1/wave2 source weight calculation was altered")
    previous_response_only=earlier_q[
        "response_only_stratum_poststratified_point_DESCRIPTIVE_NOT_DESIGN_UNBIASED"]
    if not previous_response_only < 0 < source_actual:
        raise ValueError("Original source-weighted respondent-only q sign reversal no longer holds")
    leave_one=[]
    for r in qualified:
        subgroup=[s for s in qualified if s is not r]
        q=weighted_effect(subgroup)
        if q is None:raise ValueError("Missing entire response stratum on source removal")
        leave_one.append({
           "held_out_original_source_DOI":r["source_DOI_labels"][0],
           "held_out_wave":r["wave"],
           "held_out_stratum":r["source_stratum"],
           "held_out_source_q":r["source_q_mean_equal_registered_native"],
           "remaining_respondent_only_poststratified_q":q
        })
    leave_two=[]
    for a,b in combinations(qualified,2):
        others=[r for r in qualified if r is not a and r is not b]
        q=weighted_effect(others)
        if q is not None:
            leave_two.append({
               "held_out_source_DOIs":sorted([a["source_DOI_labels"][0],b["source_DOI_labels"][0]]),
               "remaining_respondent_only_poststratified_q":q
            })
    by_status=Counter((x["wave"],x["original_source_status"]) for x in rows)
    pre_hold=sum(v for (wave,s),v in by_status.items()
                 if s in ("PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD",
                          "SOURCE_FILE_METADATA_HOLD","HUMAN_LOCATION_PRIVACY_HOLD"))
    source_hold=sum(v for (wave,s),v in by_status.items()
                    if s in ("ORIGINAL_SOURCE_GPS_HOLD","SOURCE_GPS_HOLD"))
    native_none=sum(r["native_q_supported"]==0 for r in rows
                    if r["original_source_status"] not in
                    ("PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD",
                     "SOURCE_FILE_METADATA_HOLD","HUMAN_LOCATION_PRIVACY_HOLD",
                     "ORIGINAL_SOURCE_GPS_HOLD","SOURCE_GPS_HOLD"))
    any_q=sum(r["native_q_supported"]>0 for r in rows)
    full_q=sum(r["complete_source_native_q"] for r in rows)
    full_day=sum(r["complete_source_daylength"] for r in rows)
    multi_joint=sum(r["joint_native_month_covariance_status"] in
                    ("SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED",
                     "SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED")
                    for r in rows)
    if (pre_hold!=7 or source_hold!=3 or native_none!=12
        or any_q!=13 or full_q!=10 or full_day!=8 or multi_joint!=0):
        raise ValueError("Stage-specific original source q measurement attrition changed")
    loos=[x["remaining_respondent_only_poststratified_q"] for x in leave_one]
    twos=[x["remaining_respondent_only_poststratified_q"] for x in leave_two]
    if len(leave_two)!=44: # dropping both multi-native respondent sources leaves no multi stratum
        raise ValueError("Unexpected pairwise respondent-stratum coverage")
    return {
      "version":"RSDS_SOURCE_MEASURABILITY_AND_Q_META_SIGN_FRAGILITY_V1",
      "stage":"POST_ORIGINAL_GPS_FIRST_RESULT_METHOD_SENSITIVITY_NO_NEW_EVENTS",
      "original_437_to_120_to_35_source_frame":[437,259,120,35],
      "source_study_attrition":{
         "selected_original_source_studies":35,
         "provider_file_or_privacy_HOLD":pre_hold,
         "original_GPS_source_wide_HOLD":source_hold,
         "original_GPS_parsed_but_no_native_q":native_none,
         "source_has_at_least_one_native_q":any_q,
         "source_complete_all_registered_native_q":full_q,
         "source_complete_all_native_daylength":full_day,
         "multi_native_same_calendar_month_covariance_estimated":multi_joint,
         "partial_native_q_sources":any_q-full_q,
         "remaining_untouched_single_native_metadata_source_candidates":85
      },
      "source_effect_signs":{
         "respondent_full_source_q_positive":sum(r["source_q_mean_equal_registered_native"]>0
                                                   for r in qualified),
         "respondent_full_source_q_negative":sum(r["source_q_mean_equal_registered_native"]<0
                                                   for r in qualified),
         "respondent_full_source_q_zero":sum(r["source_q_mean_equal_registered_native"]==0
                                               for r in qualified)
      },
      "response_only_poststratified_q":{
         "wave1_original_16_source_response_only":earlier_q[
           "response_only_stratum_poststratified_point_DESCRIPTIVE_NOT_DESIGN_UNBIASED"],
         "wave1_plus_wave2_35_source_response_only":source_actual,
         "sign_reversal_after_wave2_sources":True,
         "source_n_with_full_q":len(qualified),
         "source_stratum_weights":{"single_native":109/120,"multi_native":11/120},
         "leave_one_source_out_range":[min(loos),max(loos)],
         "leave_one_source_out_records":sorted(leave_one,key=lambda x:x["remaining_respondent_only_poststratified_q"]),
         "leave_two_source_out_range":[min(twos),max(twos)],
         "most_negative_two_source_removal":min(leave_two,key=lambda x:x["remaining_respondent_only_poststratified_q"]),
         "leave_two_sensitivity_count":len(leave_two),
         "hypothetical_missing_effect_sensitivity_for_35_selected":present[
             "hypothetical_fill_missing_among_selected_35_with_observed_source_range"],
         "not_population_estimand_not_confirmatory":True
      },
      "finite_120_metadata_frame_source_measurability":current[
          "finite_120_metadata_frame_source_measurability"],
      "interpretation_limits":[
        "The 35 studies include 25 with incomplete source-wide q, and among 120 metadata-eligible packages 85 additional single-native sources remain unexamined.",
        "The original source frame includes one human-location privacy HOLD; it is not an appropriate wildlife target domain. Never open its GPS or substitute another source.",
        "Post-wave1 source quota selected after examining first outcomes. Deterministic SHA order is not a randomized population survey.",
        "Leave-one/two source exclusions are deliberate post-outcome robustness diagnostics, not extra independent biological studies or cherry-picked new sampling designs.",
        "This q is GPS log1p displacement speed night-minus-day, not true nocturnal activity or causal sunlight selection.",
        "Complete respondent source effects must not be pooled with missing effect sources as if q=0 or missing at random.",
        "Older 48-source shortwave external negative sign NOT REPLICATED and 12 ERA5 original raw-SHA HOLD remain unchanged."
      ],
      "old_48_source_external_solar_direction_negative_nonreplicated":True,
      "old_12_ERA5_raw_json_source_hold":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    print("RSDS_META_Q_SIGN_FRAGILITY",json.dumps({
        "stages":result["source_study_attrition"],
        "source_signs":result["source_effect_signs"],
        "q_response_only":result["response_only_poststratified_q"][
            "wave1_plus_wave2_35_source_response_only"],
        "leave_one":result["response_only_poststratified_q"]["leave_one_source_out_range"],
        "leave_two":result["response_only_poststratified_q"]["leave_two_source_out_range"],
        "largest_negative_two_source_case":result["response_only_poststratified_q"]["most_negative_two_source_removal"]
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
