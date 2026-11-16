"""Three-study q~Tmax,L,water original within-tag-month source-only meta."""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path
from scripts.rsds_three_mule_within_month_joint_weather_v1 import (
    freeze,SOURCE_VERSION,RAW_KEYS,MONTH_KEYS,DELTA_KEYS
)
VERSION="RSDS_THREE_MULE_SOURCE_NATIVE_MONTH_JOINT_CLIMATE_FIRST_META_V1"

def validate_effect(x):
    if x["status"]!="ESTIMABLE":
        if x["estimate"] is not None or x["nominal_CI_sign"]!="UNAVAILABLE":
            raise ValueError("Unsupported original source effect assigned sign")
        return
    v=x["estimate"];ci=x["tag_bootstrap95"];n=x["original_tag_n"]
    if (not isinstance(n,int) or n<5 or not isinstance(v,(int,float)) or
        not math.isfinite(v) or not isinstance(ci,list) or len(ci)!=2 or
        not all(isinstance(a,(float,int)) and math.isfinite(a) for a in ci) or
        ci[0]>ci[1]):
        raise ValueError("Source native bootstrap original tag n/CI invalid")
    sign="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
    if x["nominal_CI_sign"]!=sign:raise ValueError("Source original tag bootstrap sign altered")


def aggregate(results):
    c,studies,prior=freeze()
    if len(results)!=3:raise ValueError("Three registered source studies mandatory")
    out=[]
    for i,(study,p,r) in enumerate(zip(studies,prior,results)):
        if (r.get("version")!=SOURCE_VERSION or
            r.get("source_index")!=i or
            r.get("source_key")!=study["source_key"] or
            r.get("source_original_DOI")!=study["doi"] or
            r.get("source_original_study_id")!=str(study["study_id"]) or
            r.get("original_source_registered_native_taxa")!=study["registered_taxa"] or
            r.get("original_48_solar_heading_signed_external_not_replicated") is not True or
            r.get("original_12_ERA5_old_raw_sha_hold") is not True or
            r.get("GPS_animal_ID_coordinates_dates_weather_grid_key_and_individual_slopes_not_exported") is not True):
            raise ValueError("Original source DOI/species/previous hold changed")
        status=r["status"]
        if status not in ("SOURCE_MONTH_JOINT_ESTIMABLE",
                          "SOURCE_MONTH_JOINT_UNAVAILABLE","SOURCE_HOLD"):
            raise ValueError("Unknown source first month climate status")
        fit=r["source_month_joint"]
        if status=="SOURCE_HOLD" and fit is not None:
            raise ValueError("Source HOLD must not carry partial effects")
        if status!="SOURCE_HOLD":
            if (not r["provider_original_MD5_verified"] or
                not r["original_native_q_and_daylength_mean_n_replayed_before_ERA5"] or
                not r["original_weather_version_exact_20261116_physical_sha_replayed"] or
                not r["prior_six_climate_model_mean_n_replayed"] or
                r["weather_grid_HOLD_count"]!=0 or
                r["weather_grid_count"]!=p["new_weather_source_grid_count"]):
                raise ValueError("Joined source GPS/weather physics insufficient for new inference")
            if fit is None or set(fit["month_partial_joint_estimands"])!=set(RAW_KEYS+MONTH_KEYS+DELTA_KEYS):
                raise ValueError("Source month model missing any prespecified component")
            names=RAW_KEYS+MONTH_KEYS+DELTA_KEYS
            for name in names:validate_effect(fit["month_partial_joint_estimands"][name])
            if status=="SOURCE_MONTH_JOINT_ESTIMABLE" and fit["eligible_original_tagged_animals_joint_all_three"]<5:
                raise ValueError("Original source month tagged sample insufficient")
            if status=="SOURCE_MONTH_JOINT_UNAVAILABLE" and fit["eligible_original_tagged_animals_joint_all_three"]>=5:
                raise ValueError("Source original month effect suppressed despite >=5 tags")
        allE=fit["month_partial_joint_estimands"] if fit else {}
        out.append({
            "source_original_DOI":study["doi"],
            "source_key":study["source_key"],
            "original_native":"Odocoileus hemionus",
            "source_status":status,
            "prior_first_Tmax_joint_source":p["source_weather_adjusted_endpoints"]["tmax_given_daylight_and_water"],
            "prior_first_daylength_joint_source":p["source_weather_adjusted_endpoints"]["daylight_given_tmax_and_water"],
            "prior_first_water_joint_source":p["source_weather_adjusted_endpoints"]["water_given_daylight_and_tmax"],
            "new_month_same_tag_model":fit,
            "source_month_Tmax_nominal_CI":allE.get("month_Tmax_given_L_water",{}).get("nominal_CI_sign","UNAVAILABLE"),
            "source_month_L_nominal_CI":allE.get("month_L_given_Tmax_water",{}).get("nominal_CI_sign","UNAVAILABLE"),
            "source_month_water_nominal_CI":allE.get("month_water_given_Tmax_L",{}).get("nominal_CI_sign","UNAVAILABLE"),
            "month_vs_raw_Tmax_paired_nominal_CI":allE.get("delta_month_minus_raw_Tmax",{}).get("nominal_CI_sign","UNAVAILABLE"),
            "exact_original_weather_physical_snapshot_replayed":r["original_weather_version_exact_20261116_physical_sha_replayed"],
            "source_q_and_old_climate_means_replayed":r["original_native_q_and_daylength_mean_n_replayed_before_ERA5"] and r["prior_six_climate_model_mean_n_replayed"],
            "old_weather_source_grid_count":p["new_weather_source_grid_count"],
            "current_weather_grid_count":r["weather_grid_count"],
            "current_weather_grid_HOLDS":r["weather_grid_HOLD_count"],
            "source_HOLD_stage":r["hold_stage"],
            "source_HOLD_type":r["hold_type"],
        })
    return {
        "version":VERSION,
        "stage":"POST_PRIOR_HEAT_POSITIVE_THREE_SOURCE_WITHIN_ORIGINAL_TAG_YEAR_MONTH_NEW_JOINT",
        "original_source_study_DOI_denominator":3,
        "all_original_three_native_source_studies_including_HOLD":out,
        "source_status_counts":dict(sorted(Counter(z["source_status"] for z in out).items())),
        "first_three_source_weather_results_git_blob_unchanged":c["original_first_three_source_result_blob"],
        "eligible_month_joint_source_studies":sum(z["source_status"]=="SOURCE_MONTH_JOINT_ESTIMABLE" for z in out),
        "nominal_within_month_Tmax_positive_source_count":sum(z["source_month_Tmax_nominal_CI"]=="POSITIVE" for z in out),
        "nominal_within_month_Tmax_negative_source_count":sum(z["source_month_Tmax_nominal_CI"]=="NEGATIVE" for z in out),
        "source_native_nonindependent_study_populations":True,
        "original_tag_month_bootstrap_nominal_not_population_random_effects_CI":True,
        "source_weather_all_grid_same_physical_sha_required_for_inference":True,
        "species_model_context_chosen_after_prior_outcomes":True,
        "no_inference_of_direct_temperature_thermal_avoidance_or_evolution":True,
        "old_solar_48_source_signed_direction_external_NOT_REPLICATED":True,
        "old_ERA5_12_original_raw_SHA_HOLD":True,
        "original_source_GPS_tracks_dates_tags_weather_grids_and_individual_estimates_not_exported":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    arg=p.parse_args()
    fs=list(arg.directory.glob("source_*.json"))
    if len(fs)!=3:raise ValueError("Exactly 3 independent original DOI source-level records")
    by={}
    for file in fs:
        d=json.loads(file.read_text(encoding="utf-8"));i=d.get("source_index")
        if type(i)!=int or i not in range(3) or i in by:
            raise ValueError("Original source result index missing/duplicated")
        by[i]=d
    result=aggregate([by[i] for i in range(3)])
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_THREE_MULE_NATIVE_MONTH_CLIMATE",json.dumps({
        "source_studies":3,
        "source_status":result["source_status_counts"],
        "within_month_positive_Tmax_sources":result["nominal_within_month_Tmax_positive_source_count"],
        "per_source":[{"source":x["source_key"],"heat":x["source_month_Tmax_nominal_CI"],
        "L":x["source_month_L_nominal_CI"],"water":x["source_month_water_nominal_CI"]}
        for x in result["all_original_three_native_source_studies_including_HOLD"]]
    },sort_keys=True),flush=True)


if __name__=="__main__":main()
