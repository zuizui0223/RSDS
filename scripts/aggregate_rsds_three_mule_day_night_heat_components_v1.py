"""Source-study-gated day/night thermal contribution across 3 original GPS DOIs."""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path

from scripts.rsds_three_mule_day_night_heat_components_v1 import (
    frozen,SOURCE_VERSION,KEYS,TOL
)
VERSION="RSDS_THREE_MULE_ORIGINAL_SOURCE_DAY_NIGHT_HEAT_COMPONENT_META_V1"

def verify_ci(d):
    if d["status"]=="ESTIMABLE":
        v=d["estimate"];ci=d["tag_bootstrap95"];n=d["original_tag_n"]
        if (not isinstance(n,int) or n<5 or
            not isinstance(v,(int,float)) or not math.isfinite(v) or
            not isinstance(ci,list) or len(ci)!=2 or
            not all(isinstance(x,(int,float)) and math.isfinite(x) for x in ci) or ci[0]>ci[1]):
            raise ValueError("Invalid original source tagged animal q component CI")
        expected="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
        if d["nominal_sign"]!=expected:
            raise ValueError("Original source component nominal CI sign changed")
        return
    if (d.get("estimate") is not None or d.get("tag_bootstrap95") is not None or
        d.get("nominal_sign")!="UNAVAILABLE"):
        raise ValueError("Missing original source DAY/NIGHT effect must not be replaced by zero")

def pattern(components):
    night=components["month_night_tmax"]["nominal_sign"]
    day=components["month_day_tmax"]["nominal_sign"]
    q=components["month_q_tmax"]["nominal_sign"]
    if q!="POSITIVE":
        return "Q_HEAT_POSITIVE_NOT_SUPPORTED_DO_NOT_ATTRIBUTE"
    if day=="NEGATIVE" and night=="POSITIVE":
        return "NEGATIVE_DAY_AND_POSITIVE_NIGHT_COEFFICIENTS"
    if day=="NEGATIVE" and night=="INCLUDES_ZERO":
        return "NEGATIVE_DAY_COEFFICIENT_NIGHT_UNRESOLVED"
    if day=="INCLUDES_ZERO" and night=="POSITIVE":
        return "POSITIVE_NIGHT_COEFFICIENT_DAY_UNRESOLVED"
    if day=="NEGATIVE" and night=="NEGATIVE":
        return "BOTH_COMPONENT_SLOPES_NEGATIVE_DAY_MORE_NEGATIVE"
    if day=="POSITIVE" and night=="POSITIVE":
        return "BOTH_COMPONENT_SLOPES_POSITIVE_NIGHT_MORE_POSITIVE"
    if day=="INCLUDES_ZERO" and night=="INCLUDES_ZERO":
        return "COMPONENT_INTERVALS_INDIVIDUALLY_UNRESOLVED_DESPITE_Q"
    return "MIXED_DIRECTION_OR_COMPONENT_UNCERTAINTY_NO_SIMPLE_MECHANISM"

def aggregate(rows):
    c,sources,old_records,prior=frozen()
    if len(rows)!=3:raise ValueError("All 3 original source study records required")
    report=[]
    for i,(r,s,p) in enumerate(zip(rows,sources,prior)):
        if (r.get("version")!=SOURCE_VERSION or r.get("source_index")!=i or
            r.get("original_study_DOI")!=s["doi"] or
            r.get("source_key")!=s["source_key"] or
            r.get("original_study_id")!=str(s["study_id"]) or
            r.get("original_registered_taxa")!=s["registered_taxa"] or
            r.get("old_48_irradiance_heading_direction_signed_test_NOT_REPLICATED") is not True or
            r.get("old_12_ERA5_raw_SHA_HOLD") is not True or
            r.get("no_original_GPS_tag_ids_positions_exact_dates_climate_grid_component_coefficients_exported") is not True):
            raise ValueError("Frozen public source ID/original native taxonomy or historical science changed")
        state=r["status"]
        if state not in ("SOURCE_HOLD","SOURCE_COMPONENT_UNAVAILABLE","SOURCE_COMPONENTS_ESTIMABLE"):
            raise ValueError("Unknown original source first component status")
        eff=r["source_component_summary"]
        if state=="SOURCE_HOLD":
            if eff is not None:
                raise ValueError("Partial original source HOLD leaked individual component")
            classification="SOURCE_HOLD_NOT_ZERO"
        else:
            if not (r["source_original_GPS_MD5_bytes_verified"] and
                    r["original_full_native_q_daylength_replayed"] and
                    r["original_source_daily_q_with_day_night_components_reconstructed"] and
                    r["same_20261116_daily_ERA5_physical_SHA_replayed"] and
                    r["original_six_weather_source_model_means_replayed"] and
                    r["source_weather_grid_HOLD_count"]==0 and eff is not None and
                    eff["all_three_endpoints_identical_tag_date_weather_support"] is True):
                raise ValueError("DAY/NIGHT source accepted without original q/weather full replay")
            if eff["eligible_original_tagged_animals"]!=p["new_month_same_tag_model"]["eligible_original_tagged_animals_joint_all_three"]:
                raise ValueError("Original month joint tag response subset changed")
            if state=="SOURCE_COMPONENTS_ESTIMABLE":
                if (eff["source_q_temperature_exact_mean_and_95_CI_replayed"] is not True or
                    eff["q_night_minus_day_algebra_exact_in_all_tag_years_and_bootstrap"] is not True or
                    set(eff["components"])!=set(KEYS)):
                    raise ValueError("Original algebra/bootstrap q native identities missing")
                for key in KEYS:
                    verify_ci(eff["components"][key])
                for prefix in ("raw","month"):
                    b=eff["components"]
                    if abs(b[prefix+"_night_tmax"]["estimate"]-
                           b[prefix+"_day_tmax"]["estimate"]-
                           b[prefix+"_q_tmax"]["estimate"])>TOL:
                        raise ValueError("Source q heat coefficient != NIGHT minus DAY")
                classification=pattern(eff["components"])
            else:
                if eff["eligible_original_tagged_animals"]>=5:
                    raise ValueError("Original q and component estimates suppressed despite tag support")
                classification="SOURCE_COMPONENT_INSUFFICIENT_TAGS_NOT_ZERO"
        vals=(eff or {}).get("components") or {}
        report.append({
            "source_index":i,
            "source_original_DOI":s["doi"],
            "source_key":s["source_key"],
            "source_native_taxon":"Odocoileus hemionus",
            "source_status":state,
            "source_original_q_daylength_model_replayed":r["original_full_native_q_daylength_replayed"],
            "source_original_q_DAY_NIGHT_algebra_reconstructed":r["original_source_daily_q_with_day_night_components_reconstructed"],
            "first_2026_11_16_ERA5_physical_fingerprint_replayed":r["same_20261116_daily_ERA5_physical_SHA_replayed"],
            "original_six_source_q_weather_model_means_replayed":r["original_six_weather_source_model_means_replayed"],
            "original_q_month_bootstrap_mean_CI_replayed":r["original_q_Tmax_month_first_bootstrap_replayed"],
            "source_original_individual_tag_n":eff["eligible_original_tagged_animals"] if eff else None,
            "source_original_month_replayed_q_Tmax":p["new_month_same_tag_model"]["month_partial_joint_estimands"]["month_Tmax_given_L_water"],
            "same_original_individual_tag_day_month_DAY_Tmax":vals.get("month_day_tmax"),
            "same_original_individual_tag_day_month_NIGHT_Tmax":vals.get("month_night_tmax"),
            "same_original_individual_tag_day_month_Q_Tmax":vals.get("month_q_tmax"),
            "same_original_individual_tag_day_month_component_regression_means_and_paired_95":eff,
            "source_original_temperature_components_pattern":classification,
            "original_source_weather_grid_count":r["source_weather_grid_count"],
            "original_source_weather_grid_HOLD_count":r["source_weather_grid_HOLD_count"],
            "source_HOLD_stage":r["hold_stage"],
            "source_HOLD_type":r["hold_type"],
        })
    usable=[r for r in report if r["source_status"]=="SOURCE_COMPONENTS_ESTIMABLE"]
    signs={
        "source_Q_Tmax_nominal_positive":sum(r["same_original_individual_tag_day_month_Q_Tmax"]["nominal_sign"]=="POSITIVE" for r in usable),
        "source_DAY_Tmax_nominal_negative":sum(r["same_original_individual_tag_day_month_DAY_Tmax"]["nominal_sign"]=="NEGATIVE" for r in usable),
        "source_DAY_Tmax_nominal_positive":sum(r["same_original_individual_tag_day_month_DAY_Tmax"]["nominal_sign"]=="POSITIVE" for r in usable),
        "source_NIGHT_Tmax_nominal_positive":sum(r["same_original_individual_tag_day_month_NIGHT_Tmax"]["nominal_sign"]=="POSITIVE" for r in usable),
        "source_NIGHT_Tmax_nominal_negative":sum(r["same_original_individual_tag_day_month_NIGHT_Tmax"]["nominal_sign"]=="NEGATIVE" for r in usable),
    }
    return {
        "version":VERSION,
        "stage":"POST_PRIOR_POSITIVE_Q_HEAT_SOURCE_NATIVE_DAY_NIGHT_EXACT_SAME_TAG_DAY_DECOMPOSITION",
        "registered_original_source_DOI_studies":3,
        "all_three_original_source_studies_including_SOURCE_HOLD":report,
        "original_source_status_counts":dict(sorted(Counter(x["source_status"] for x in report).items())),
        "same_source_original_tag_date_weather_algebra_eligible_studies":len(usable),
        "original_tag_source_component_nominal_CI_sign_counts":signs,
        "by_source_component_interpretation":{
            r["source_key"]:r["source_original_temperature_components_pattern"]
            for r in report
        },
        "previous_first_month_q_heat_source_blob_unchanged":
            c["original_old_three_source_month_climate_first_result_blob"],
        "not_a_prospective_or_world_random_wildlife_test":True,
        "not_causal_temperature_behavioral_heat_avoidance":True,
        "source_populations_independence_unverified":True,
        "prior_48_original_source_irradiance_heading_signed_result_NOT_REPLICATED":True,
        "prior_12_ERA5_old_raw_SHA_cases_HOLD":True,
        "original_GPS_day_night_tag_positions_dates_individual_effects_not_exported":True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--directory",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    arg=parser.parse_args()
    fs=list(arg.directory.glob("source_*.json"))
    if len(fs)!=3:raise ValueError("Exactly 3 original source study first result files required")
    by={}
    for f in fs:
        x=json.loads(f.read_text(encoding="utf-8"));i=x.get("source_index")
        if type(i)!=int or i not in (0,1,2) or i in by:
            raise ValueError("Duplicate/missing original source study identity")
        by[i]=x
    result=aggregate([by[i] for i in range(3)])
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_DAY_NIGHT_TMAX_THREE_SOURCE_META",json.dumps({
        "source_studies":3,
        "status":result["original_source_status_counts"],
        "sign_counts":result["original_tag_source_component_nominal_CI_sign_counts"],
        "classified":result["by_source_component_interpretation"],
        "original_q_exact_algebra_required":True
    },sort_keys=True),flush=True)


if __name__=="__main__":main()
