"""SHA-frozen original ecological interpretation of same-month shared GPS calendars.

This is a SECONDARY source-result audit; no animal GPS or raw month/year
labels are read. Protect previous Illinois and Wyoming unadjusted slopes,
study DOI/taxon roster, all six new month-conditional source contrasts,
and historical negative external shortwave/ERA5 raw JSON version HOLDs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PIN={
"research/RSDS_ILLINOIS_THREE_NATIVE_SHARED_DATE_V1_FIRST_RESULT_20261029.json":"924156e8f98d0bb9fda2f9d23aabc03060af2a31",
"research/RSDS_WYOMING_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json":"9b04b18ee73227b4e7a1715ebfb97f610af66952",
"research/RSDS_SHARED_MONTH_ADJUSTED_PHOTOPERIOD_TWO_STUDY_V1_FIRST_RESULT_20261029.json":"ce14cf34bcb1237671a8ad4c5f325bd91eb35356"
}
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
OUT=ROOT/"research/RSDS_SHARED_MONTH_PHOTOPERIOD_VALIDATED_ECOLOGY_20261029.json"


def locked(root=ROOT):
    result={}
    for p,expected in PIN.items():
        raw=(root/p).read_bytes()
        got=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
        if got!=expected:
            raise ValueError("A first original GPS source ecological result changed: "+p)
        result[p]=json.loads(raw.decode("utf-8"))
    return list(result.values())


def source_interpretation(root=ROOT):
    ill,wy,first=locked(root)
    if (ill["study_DOI"]!="10.5441/001/1.649"
            or wy["study_Doi"]!="10.5441/001/1.703"
            or first["original_studies"]!=2
            or first["original_source_pair_clock_comparisons_expected"]!=6
            or first["original_source_pair_clock_comparisons_reported"]!=6
            or first["original_study_native_taxa"]!={"Illinois":3,"Wyoming":2}
            or not first["old_12_ERA5_raw_JSON_HOLDS_UNCHANGED"]
            or first["old_48_signed_radiation_external_nonreplication_UNCHANGED"]!=NEGATIVE):
        raise ValueError("Previous ecological source study/negative verdict altered")
    rows=[]
    signed=Counter()
    pair_keys=set()
    old_repro=0
    for study in first["source_study_results"]:
        reg=study["original_study_region"]
        if reg not in ("Illinois","Wyoming") or not study["MD5_and_native_source_effects_replayed"]:
            raise ValueError("New source month inference without original source MD5/q/daylight replay")
        count=0
        for x in study["native_pair_schedule_results"]:
            a,b,schedule=x["native_first"],x["native_second"],x["schedule"]
            k=(reg,a,b,schedule)
            if k in pair_keys:raise ValueError("Duplicated source-native calendar pair")
            pair_keys.add(k)
            if reg=="Illinois":
                old=ill["three_source_native_pair_q_and_daylength"][f"{a} MINUS {b}"][
                    "daily_mean_q_difference_vs_standard_daylength"]
                old_value=old["slope"]
                old_state=old["status"]
            else:
                old=wy["calendar_date_paired_source_difference"]["slopes"][schedule]
                old_value=old["mule_q_minus_bighorn_q_vs_shared_reference_daylength"]
                old_state=old["status"]
            if (not math.isclose(old_value,x["original_unadjusted_daylength_slope"],rel_tol=0,abs_tol=1e-12)
                    or x["original_unadjusted_source_nominal_sign"]!=old_state):
                raise ValueError("Previously observed original native q~daylength source coefficient was modified")
            old_repro+=1
            fresh=x["within_year_month_demeaned_q_daylength"]
            state=fresh["nominal_interval_sign"]
            if state not in ("NEGATIVE","POSITIVE","INCLUDES_ZERO","UNAVAILABLE"):
                raise ValueError("Unknown original month-block source result")
            if state!="UNAVAILABLE":
                if (fresh["status"]!="MONTH_CONDITIONAL_DESCRIPTIVE_ESTIMABLE"
                        or fresh["eligible_source_month_blocks"]<5
                        or fresh["eligible_source_calendar_days"]<40
                        or fresh["within_month_daylength_centered_sum_squares_hours2"]<.25
                        or fresh["bootstrap_finite_replicates"]<1900):
                    raise ValueError("Within-month slope identified without adequate year-month source exposure")
                ci=fresh["month_cluster_bootstrap95"]
                y=fresh["within_month_daylength_slope"]
                if not isinstance(ci,list) or len(ci)!=2 or not all(math.isfinite(v) for v in (y,*ci)):
                    raise ValueError("Month-block bootstrap or estimand undefined despite supported source")
                signed[state]+=1
            else:
                if fresh["within_month_daylength_slope"] is not None or fresh["month_cluster_bootstrap95"] is not None:
                    raise ValueError("Unavailable monthly animal source contrast converted to biological zero")
                signed["UNAVAILABLE"]+=1
            count+=1
            rows.append({
                "original_study_doi":study["original_study_DOI"],
                "region":reg,
                "taxon_first":a,"taxon_second":b,"GPS_timestamp_schedule":schedule,
                "same_original_source_calendar_days":x["source_original_calendar_days"],
                "original_source_month_blocks":study["source_original_year_month_blocks"],
                "month_conditioned_year_month_blocks":fresh["eligible_source_month_blocks"],
                "month_conditioned_source_calendar_days":fresh["eligible_source_calendar_days"],
                "unadjusted_daylength_slope":old_value,
                "unadjusted_source_nominal_status":old_state,
                "month_conditioned_daylength_slope":fresh["within_month_daylength_slope"],
                "month_block_bootstrap95":fresh["month_cluster_bootstrap95"],
                "month_block_leave_one_out_range":fresh["leave_one_month_block_out_range"],
                "month_conditioned_source_nominal_sign":state,
                "interpretation":x["before_vs_after_calendar_month_control"]
            })
        if count!=(3 if reg=="Illinois" else 3):
            raise ValueError("Original three native-pairs or three GPS clock schedules dropped")
    if len(rows)!=6 or old_repro!=6:
        raise ValueError("Source first seasonal coefficients not fully reproduced")
    if signed!={"INCLUDES_ZERO":4,"NEGATIVE":1,"POSITIVE":1}:
        raise ValueError("The first within-calendar-month source experimental signs changed")
    expected_persistent={
        ("Illinois","Lynx rufus","Odocoileus virginianus","original"):"NEGATIVE",
        ("Wyoming","Odocoileus hemionus","Ovis canadensis","original"):"POSITIVE",
    }
    if { (x["region"],x["taxon_first"],x["taxon_second"],x["GPS_timestamp_schedule"]):x["month_conditioned_source_nominal_sign"]
         for x in rows if x["month_conditioned_source_nominal_sign"] in ("POSITIVE","NEGATIVE")
       } != expected_persistent:
        raise ValueError("The two original source conditional positive/negative native taxa comparisons changed")
    return {
        "version":"RSDS_COOBSERVED_SOURCE_WITHIN_MONTH_PHOTOPERIOD_VALIDATED_V1",
        "date":"2026-10-29",
        "stage":"SOURCE_FIRST_TWO_STUDY_MONTH_CONFOUNDED_SEASONAL_ECOLOGY_SENSITIVITY",
        "original_source_field_studies":2,
        "source_native_taxon_and_schedule_comparisons":6,
        "all_original_unadjusted_native_daylength_slopes_replayed":6,
        "month_conditioned_comparisons_nominal_bootstrap_interval_signs":dict(signed),
        "nominal_source_intervals_not_spanning_zero":len(expected_persistent),
        "source_case_results":rows,
        "old_12_ERA5_source_provider_raw_SHA_HOLDS_UNCHANGED":True,
        "old_external_48_source_negative_radiation_direction_nonreplication_UNCHANGED":NEGATIVE,
        "original_source_first_input_git_blobs":PIN,
        "scientific_claim_ceiling":[
            "Within-source calendar-month daylength association persists in only two of six predeclared study×pair×schedule source conditional comparisons; four other estimates are not numerically zero but have descriptive 95% month-block intervals spanning zero.",
            "Both two zero-excluding month-block bootstrap intervals border zero, and only eight eligible month-year blocks from each original study support these exploratory results. No multiplicity correction, causal effect or world generality follows.",
            "Some adjusted point slopes are numerically larger than old all-month slopes; a zero-spanning CI does NOT show between-month seasons entirely explain the former associations.",
            "A within-calendar-month daylength trend is still almost a deterministic day-of-month trend; weather, food, breeding, tag composition, source GPS duty-cycle and spatial route remain confounded.",
            "Illinois three taxa and Wyoming two taxa are each from exactly one original Movebank source study, not five independent new ecological populations.",
            "Historical first native source effects, old raw ERA5 JSON 12 holds and external 48 study negative signed radiation replication unchanged.",
            "No original GPS tracking IDs, coordinates, date/month/year labels or source-private daily q data are published."
        ],
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,default=OUT)
    args=cli.parse_args()
    d=source_interpretation()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_MONTH_PHOTOPERIOD_ORIGINAL_ECOLOGY_REPLAY",json.dumps({
        "studies":d["original_source_field_studies"],
        "pair_schedules":d["source_native_taxon_and_schedule_comparisons"],
        "old_slopes_replayed":d["all_original_unadjusted_native_daylength_slopes_replayed"],
        "new_signs":d["month_conditioned_comparisons_nominal_bootstrap_interval_signs"],
        "old_negative":d["old_external_48_source_negative_radiation_direction_nonreplication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
