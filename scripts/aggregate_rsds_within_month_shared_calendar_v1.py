"""Source-complete two-study original shared-calendar month-adjusted astronomy results.

Retain all Illinois three native species pair and Wyoming mule deer/bighorn
three clock schedules, including any observation support/seasonal HOLDS.
This is exploratory within source-year-month centering, not a causal
daylength response or population-general macro effect.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"research/RSDS_SHARED_MONTH_ADJUSTED_PHOTOPERIOD_SOURCE_V1_CONTRACT.json"
RESULT=ROOT/"research/RSDS_SHARED_MONTH_ADJUSTED_PHOTOPERIOD_TWO_STUDY_V1_FIRST_RESULT_20261029.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def aggregate(studies):
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    fixed={x["region"]:x for x in spec["fixed_original_sources"]}
    if (spec["version"]!="RSDS_WITHIN_MONTH_SHARED_CALENDAR_PHOTOPERIOD_V1"
            or set(studies)!=set(fixed) or set(fixed)!={"Illinois","Wyoming"}):
        raise ValueError("Two original fixed source field studies absent")
    records=[]
    status=Counter()
    native=Counter()
    states=Counter()
    supported=0
    for region in ("Illinois","Wyoming"):
        z=studies[region]
        ref=fixed[region]
        if (z.get("region")!=region or z.get("source_DOI")!=ref["doi"]
                or z.get("native_taxa")!=ref["taxa"]
                or not z.get("historic_12_ERA5_source_HOLDS_UNCHANGED")
                or z.get("historic_48_external_sunlight_signed_nonreplication_UNCHANGED")!=NEGATIVE):
            raise ValueError("Source DOI/taxon/historical negative radiation source changed")
        rstatus=z["status"]
        status[rstatus]+=1
        if rstatus in ("ORIGINAL_IL_THREE_NATIVE_MONTH_CONTRAST_COMPLETE",
                       "ORIGINAL_WY_TWO_NATIVE_MONTH_CONTRAST_COMPLETE"):
            if not z.get("provider_original_file_MD5_replayed") or not z.get(
                "prior_native_q_and_daylength_replayed_before_month_diagnosis"):
                raise ValueError("New within-month source result before MD5/q/daylength replay")
            if z.get("same_original_calendar_days")!=ref["shared_native_calendar_days"] or z.get(
                "source_year_month_blocks")!=ref["year_month_blocks"]:
                raise ValueError("Previous same-source calendar/month support changed")
        pairs=z.get("source_native_pair_month_contrasts",[])
        expected={(a,b,sch) for a,b in ref["pairs"] for sch in ref["schedules"]}
        if len(pairs)!=len(expected) or {
            (a.get("native_first"),a.get("native_second"),a.get("source_GPS_clock_schedule"))
            for a in pairs
        }!=expected:
            if rstatus in ("ORIGINAL_IL_THREE_NATIVE_MONTH_CONTRAST_COMPLETE",
                           "ORIGINAL_WY_TWO_NATIVE_MONTH_CONTRAST_COMPLETE"):
                raise ValueError("One of fixed six native study/schedule pairs missing")
        full=[]
        for p in pairs:
            now=p["within_calendar_month_analysis"]
            if (p["same_study_original_calendar_days"]!=ref["shared_native_calendar_days"]
                    or p["archived_first_unadjusted_original_coefficient_replayed"]!="ARCHIVED_UNADJUSTED_SLOPE_REPLAYED"):
                raise ValueError("Previous source datewise seasonal contrast was not numerically replayed")
            status_new=now["nominal_interval_sign"]
            if status_new not in ("UNAVAILABLE","POSITIVE","NEGATIVE","INCLUDES_ZERO"):
                raise ValueError("Unexpected within-month source status")
            if status_new=="UNAVAILABLE":
                if now["within_month_daylength_slope"] is not None or now["month_cluster_bootstrap95"] is not None:
                    raise ValueError("Missing source/month support converted to fake biological q")
            else:
                n=now["eligible_source_month_blocks"]
                k=now["eligible_source_calendar_days"]
                sxx=now["within_month_daylength_centered_sum_squares_hours2"]
                slope=now["within_month_daylength_slope"]
                ci=now["month_cluster_bootstrap95"]
                if (n<5 or k<40 or sxx<.25 or ci is None or len(ci)!=2 or
                        not all(math.isfinite(v) for v in [slope,*ci])):
                    raise ValueError("Within-month inference without source block/daylight variation support")
                sign=("POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO")
                if sign!=status_new:
                    raise ValueError("Month bootstrap sign classification altered")
                supported+=1
            states[status_new]+=1
            full.append({
                "native_first":p["native_first"],
                "native_second":p["native_second"],
                "schedule":p["source_GPS_clock_schedule"],
                "source_original_calendar_days":p["same_study_original_calendar_days"],
                "original_unadjusted_daylength_slope":p["previous_all_months_source_daylength_contrast_coefficient"],
                "original_unadjusted_source_nominal_sign":p["previous_unadjusted_source_nominal_sign"],
                "within_year_month_demeaned_q_daylength":now,
                "before_vs_after_calendar_month_control":p["between_vs_within_month_status"]
            })
        records.append({
            "original_study_region":region,
            "original_study_DOI":z["source_DOI"],
            "original_native_taxa":z["native_taxa"],
            "status":rstatus,
            "MD5_and_native_source_effects_replayed":z.get("provider_original_file_MD5_replayed") and z.get(
                "prior_native_q_and_daylength_replayed_before_month_diagnosis"),
            "native_pair_schedule_results":full,
            "source_original_shared_calendar_days":z.get("same_original_calendar_days"),
            "source_original_year_month_blocks":z.get("source_year_month_blocks"),
            "hold_error_type":z.get("hold_error_type"),
            "hold_stage":z.get("hold_stage")
        })
        native[region]=len(ref["taxa"])
    if len(records)!=2:
        raise ValueError("Wrong number original field studies")
    return {
        "version":"RSDS_ILLINOIS_WYOMING_WITHIN_MONTH_PHOTOPERIOD_TWO_SOURCE_V1",
        "date":"2026-10-29",
        "stage":"POST_COOBSERVED_NATIVE_SOURCE_PHOTOPERIOD_RESULTS_WITHIN_MONTH_SEASONAL_CONFOUND_CHECK",
        "original_studies":2,
        "original_study_native_taxa":dict(native),
        "original_source_pair_clock_comparisons_expected":6,
        "original_source_pair_clock_comparisons_reported":sum(len(x["native_pair_schedule_results"]) for x in records),
        "contrasts_with_estimable_month_adjusted_daylength_slope":supported,
        "month_conditional_nominal_bootstrap_sign_counts":dict(states),
        "original_study_source_states":dict(status),
        "source_study_results":records,
        "prior_first_original_study_results_untouched":True,
        "old_12_ERA5_raw_JSON_HOLDS_UNCHANGED":True,
        "old_48_signed_radiation_external_nonreplication_UNCHANGED":NEGATIVE,
        "interpretation_limits":[
            "The 6 pair×schedule rows represent ONLY TWO original field-study DOI units, with correlated native species sampled within each field-study, NOT six independent world populations.",
            "The same original source calendar was fixed before examining these WITHIN-MONTH slope values, but the study/pair selection is post prior-sign and is exploratory not blinded replication.",
            "Year-month intercept adjustment removes across-month differences, not individual DATE trends within month. Astronomical daylength is deterministic in calendar day and cannot be isolated causally from other seasonal cues without independent conditions.",
            "Only nine year-month source blocks Illinois and eight Wyoming support original estimates. Month block bootstrap 95% is descriptive, unadjusted, and not a population confidence interval.",
            "Same shared source date does not mean the same actual GPS event time, collar acquisition schedule, habitat/weather or predator-prey interaction across different tagged original animals.",
            "NIGHT-DAY GPS q is a sampled chord-displacement contrast, NOT true activity duration, direct sun cue use, interspecific competition or predation.",
            "Missing or weak within-month astronomical support yields UNAVAILABLE, not zero. Do not select favorable date subsets post hoc.",
            "Previous 48 source irradiance direction external negative replication and original ERA5 response raw-SHA twelve historical holds unchanged.",
            "No private animal IDs, original GPS locations, precise source date/month-year labels, raw telemetry or individual q slopes are republished."
        ]
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",type=Path,required=True)
    parser.add_argument("--out",type=Path,default=RESULT)
    args=parser.parse_args()
    sources={}
    for k in ("illinois","wyoming"):
        p=args.inputs/(k+"_shared_month.json")
        if not p.is_file(): raise ValueError("Required frozen source output missing: "+k)
        sources[k.capitalize()]=json.loads(p.read_text(encoding="utf-8"))
    result=aggregate(sources)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_WITHIN_MONTH_SHARED_SOURCE_ECOLOGY",json.dumps({
        "field_studies":result["original_studies"],
        "pair_schedule_results":result["original_source_pair_clock_comparisons_reported"],
        "supported":result["contrasts_with_estimable_month_adjusted_daylength_slope"],
        "signs":result["month_conditional_nominal_bootstrap_sign_counts"],
        "old_external_negative_unchanged":result["old_48_signed_radiation_external_nonreplication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
