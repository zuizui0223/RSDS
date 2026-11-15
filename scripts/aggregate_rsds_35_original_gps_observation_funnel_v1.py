"""Aggregate privacy-safe ORIGINAL source GPS observation funnel (35 studies).

Source groups, not GPS fixes, are the units of description. No outcome- or
source-dependent thresholds are adjusted. Original archived response labels,
GPS parser failures, privacy/source-file HOLDs remain unchanged. The source
time/solar/geometry rates are conditional sampling diagnostics and need not
reflect animal ecology or global time-niche prevalence.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from statistics import mean,median

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned

CONTRACT=HERE/"research/RSDS_35_SOURCE_GPS_OBSERVATION_FUNNEL_V1_CONTRACT.json"
EARLIER=("research/RSDS_437_WAVE1_PLUS_WAVE2_35_ORIGINAL_GPS_FIRST_META_RESULT_20261115.json",
         "9b02774f7c31d29b2eec5039a80431a895dccb36")
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
METRICS=(
    "original_time_0p5_to_2h_over_adjacent",
    "same_solar_date_daypart_over_time_eligible",
    "plausible_speed_over_solar_eligible",
    "total_accepted_over_adjacent",
    "UTC_even_hour_fraction",
)
COUNT_KEYS=(
    "original_adjacent_source_fix_intervals",
    "accepted_original_0p5_to_2h_time_intervals",
    "accepted_same_solar_date_and_daypart_intervals",
    "accepted_finite_speed_max60kmh_intervals",
)


def between_zero_one(value):
    return value is None or (isinstance(value,(float,int)) and
           math.isfinite(value) and 0.0<=value<=1.0)


def source_summary(profiles):
    """One aggregate per source. All native species are retained."""
    if not profiles:return None
    result={"n_registered_native_taxa":len(profiles)}
    for key in METRICS:
        present=[p["original_observation_gate_ratios"][key]
                 for p in profiles
                 if p["original_observation_gate_ratios"][key] is not None]
        result["source_native_equal_"+key]=mean(present) if present else None
        result["source_min_native_"+key]=min(present) if present else None
        result["source_native_metric_supported_"+key]=len(present)
    for k in COUNT_KEYS:
        result["source_total_"+k]=sum(p["source_original_support_counts"][k]
                                      for p in profiles)
    result["source_total_tags_q_eligible"]=sum(
       p["source_original_support_counts"]["accepted_native_tags_with_ge10_daypart_dates"]
       for p in profiles)
    result["source_total_original_tags"]=sum(
       p["source_original_support_counts"]["source_original_tagged_animals_in_hourly_parsed_track"]
       for p in profiles)
    result["source_registered_native_q_supported"]=sum(
        p["replayed_archived_native_q_support"]["q_estimable"] for p in profiles)
    return result


def summarize(rows):
    if len(rows)!=35:
        raise ValueError("All original 35 source groups required")
    previous=pinned(*EARLIER)
    old=previous["combined_source_study_summary_rows"]
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (len(old)!=35 or previous["combined_original_registered_native_rows"]!=50
        or c["original_35_source_study_ids"]!=35
        or c["original_GPS_source_wide_HOLD"]!=3
        or c["previous_48_signed_sunlight_external"]!=NEG
        or c["previous_ERA5_raw_SHA_12_HOLD"] is not True):
        raise ValueError("Immutable original GPS study frame/legacy null changed")
    frame=[];status=Counter();all_native=[]
    for i,(current,expected) in enumerate(zip(rows,old)):
        wave,index=("wave1",i) if i<16 else ("wave2",i-16)
        if (current["wave"]!=wave or current["index"]!=index
            or current["source_original_DOI"]!=expected["source_DOI_labels"][0]
            or current["source_repository_study_id"]!=expected["source_original_repository_study_id"]
            or len(current["original_registered_native_taxa"])!=expected["registered_native_count"]
            or current["prior_first_source_q_status"]!=expected["original_source_status"]
            or current["original_first_q_study_complete"]!=expected["complete_source_native_q"]
            or current["original_first_daylength_study_complete"]!=expected["complete_source_daylength"]
            or current["original_48_signed_sunlight_external_NONREPLICATION_UNCHANGED"]!=NEG
            or current["original_12_raw_ERA5_source_HOLD_UNCHANGED"] is not True
            or current["no_individual_ID_GPS_location_date_timestamp_or_q_values_exported"] is not True):
            raise ValueError("Original registered source or first q status altered")
        parsed=(expected["original_source_status"] not in (
            "PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD","SOURCE_FILE_METADATA_HOLD",
            "HUMAN_LOCATION_PRIVACY_HOLD","ORIGINAL_SOURCE_GPS_HOLD",
            "SOURCE_GPS_HOLD"))
        success=(current["source_status"]==
                 "REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL")
        if parsed!=success:
            raise ValueError("Previously parsed original GPS source was promoted to HOLD or vice versa")
        profiles=current["original_native_funnel"]
        if parsed:
            if (len(profiles)!=expected["registered_native_count"]
                or not current["original_MD5_replayed"]
                or current["original_raw_original_event_count"]!=
                   expected["original_source_event_row_count"]
                or [p["native_taxon"] for p in profiles]!=
                   current["original_registered_native_taxa"]
                or sum(p["replayed_archived_native_q_support"]["q_estimable"]
                       for p in profiles)!=expected["native_q_supported"]
                or sum(p["replayed_archived_native_q_support"]["daylength_estimable"]
                       for p in profiles)!=expected["native_daylength_supported"]):
                raise ValueError("Original GPS source-native measured support not replayed")
            for p in profiles:
                counts=p["source_original_support_counts"]
                vals=[counts[k] for k in COUNT_KEYS]
                if vals!=sorted(vals,reverse=True):
                    raise ValueError("Original native observation gate not nested")
                if any(not between_zero_one(p["original_observation_gate_ratios"][k])
                       for k in METRICS):
                    raise ValueError("Invalid original native cadence gate proportion")
                if (counts["UTC_even_hour_selected_fix_count"]
                    +counts["UTC_odd_hour_selected_fix_count"]!=
                    counts["original_selected_UTC_hour_fix_count"]):
                    raise ValueError("Original UTC even/odd monitoring phases do not partition")
                all_native.append(p)
        else:
            if (profiles or current["original_MD5_replayed"] or
                current["original_raw_original_event_count"] is not None):
                raise ValueError("Archived original source HOLD unexpectedly read GPS")
        status[current["source_status"]]+=1
        per=source_summary(profiles)
        frame.append({
           "source_DOI":current["source_original_DOI"],
           "source_study_ID":current["source_repository_study_id"],
           "wave":wave,"stratum":expected["source_stratum"],
           "previous_original_source_status":expected["original_source_status"],
           "native_registered":expected["registered_native_count"],
           "complete_native_q":expected["complete_source_native_q"],
           "complete_native_daylength":expected["complete_source_daylength"],
           "diagnostic_status":current["source_status"],
           "native_observation_funnel_summary":per
        })
    if (status!=Counter({
        "REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL":25,
        "ARCHIVED_ORIGINAL_SOURCE_HOLD_UNCHANGED":9,
        "HUMAN_LOCATION_PRIVACY_HOLD_NO_ACCESS":1,
    }) or len({x["source_study_ID"] for x in frame})!=35
        or sum(x["native_registered"] for x in frame)!=50):
        raise ValueError("Original full source stage count lost or substituted")
    strata={}
    for group in ("complete_native_q","not_complete_native_q"):
        eligible=[x for x in frame if x["native_observation_funnel_summary"] is not None
                  and bool(x["complete_native_q"])==(group=="complete_native_q")]
        sub={"parsed_source_study_count":len(eligible)}
        for key in METRICS:
            measurements=[x["native_observation_funnel_summary"][
                              "source_native_equal_"+key] for x in eligible]
            finite=[v for v in measurements if v is not None]
            sub["source_equal_median_"+key]=median(finite) if finite else None
            sub["source_equal_range_"+key]=[min(finite),max(finite)] if finite else None
            sub["source_studies_with_measurement_"+key]=len(finite)
        strata[group]=sub
    source_all=[x for x in frame if x["native_observation_funnel_summary"] is not None]
    totals={k:sum(x["native_observation_funnel_summary"]["source_total_"+k]
                  for x in source_all) for k in COUNT_KEYS}
    out={
        "version":"RSDS_ORIGINAL_35_SOURCE_GPS_OBSERVATION_GATE_FUNNEL_V1",
        "stage":"POST_FIRST_WAVE_OUTCOME_METHOD_DIAGNOSTIC_NO_NEW_BIOLOGICAL_ENDPOINT",
        "original_public_Movebank_source_frame":[437,259,120,35],
        "source_original_registered_native_rows":50,
        "first_source_study_q_complete":10,
        "first_source_study_daylength_complete":8,
        "original_source_funnel_replayed":25,
        "original_preaccess_file_or_privacy_HOLDS":7,
        "original_sourcewide_GPS_HOLDS":3,
        "original_human_GPS_source_opened":False,
        "individual_original_source_diagnostic_rows":frame,
        "original_native_aggregate_metric_rows":all_native,
        "original_source_diagnostic_status_counts":dict(sorted(status.items())),
        "source_equal_observation_gate_contrasts_by_archived_q_availability":strata,
        "original_per_fix_gate_census_NESTED_but_NOT_independent_samples":totals,
        "scientific_claim_ceiling":[
           "This is retrospective original GPS timestamp/time/daypart/movement-speed data usability, not an experimental manipulation of sampling cadence or direct causal daylight response.",
           "All original q effects remain unchanged. Time gate, astronomy continuity and speed gates have distinct denominators.",
           "GPS fixes, segments and nested native taxa are NOT independent study units; no event-level p-values.",
           "Postoutcome grouping by q supported/unsupported is descriptive and does not establish that cadence caused missingness.",
           "One human source privacy hold, six provider size holds and three original source-wide GPS holds retain their original IDs/statuses.",
           "Previously selected 35 source DOIs are not a random global wildlife ecology source sample; two-wave quota was adjusted after Wave1.",
           "Original 48-source negative sunlight-direction external effect remains NOT REPLICATED; 12 old ERA5 raw JSON HOLDS remain."
        ],
        "original_48_solar_negative_external_NONREPLICATION_UNCHANGED":NEG,
        "old_12_ERA5_raw_source_versions_HOLD_UNCHANGED":True
    }
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    rows=[
        json.loads((a.directory/("wave1" if i<16 else "wave2")/
           ("source_%02d.json"%(i if i<16 else i-16))).read_text(encoding="utf-8"))
        for i in range(35)
    ]
    result=summarize(rows)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                     encoding="utf-8")
    print("RSDS_35_ORIGINAL_SOURCE_OBSERVATION_FUNNEL",json.dumps({
       "source_gates":result["original_source_diagnostic_status_counts"],
       "q_group_ratios":result["source_equal_observation_gate_contrasts_by_archived_q_availability"],
       "source_stage_totals":result["original_per_fix_gate_census_NESTED_but_NOT_independent_samples"],
       "not_causal_or_new_biological_result":True
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
