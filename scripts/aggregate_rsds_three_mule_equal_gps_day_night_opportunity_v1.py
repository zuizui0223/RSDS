"""Study-level exact source 3-step DAY/NIGHT equal-observation sensitivity/HOLD."""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path
from scripts.rsds_three_mule_equal_gps_day_night_opportunity_v1 import (
    frozen,SOURCE_VERSION,KEYS,TOL
)

VERSION="RSDS_THREE_MULE_SOURCE_EQUAL_GPS_INTERVAL_DAY_NIGHT_FIRST_META_V1"

def validate(x):
    if (x.get("status")!="ESTIMABLE" or
        type(x.get("source_tag_n")) is not int or x["source_tag_n"]<5 or
        not isinstance(x.get("mean"),(float,int)) or
        not math.isfinite(x["mean"]) or
        not isinstance(x.get("tag_bootstrap95"),list) or
        len(x["tag_bootstrap95"])!=2 or
        not all(isinstance(v,(int,float)) and math.isfinite(v)
                for v in x["tag_bootstrap95"]) or
        x["tag_bootstrap95"][0]>x["tag_bootstrap95"][1]):
        raise ValueError("Bad original source-tag q/3-step bootstrap first result")
    ci=x["tag_bootstrap95"]
    sign="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
    if x["nominal_sign"]!=sign:
        raise ValueError("Original source q/3step sign changed")

def aggregate(rows):
    c,source,old,prior=frozen()
    if len(rows)!=3:raise ValueError("All 3 original source DOIs required")
    table=[]
    for i,(r,s,p) in enumerate(zip(rows,source,prior)):
        if (r.get("version")!=SOURCE_VERSION or r.get("source_index")!=i
            or r.get("source_key")!=s["source_key"] or
            r.get("original_study_DOI")!=s["doi"] or
            r.get("original_study_ID")!=str(s["study_id"]) or
            r.get("source_native")!="Odocoileus hemionus" or
            r.get("original_registered_native_taxa")!=s["registered_taxa"] or
            r.get("source_48_irradiance_heading_direction_external_sign_NOT_REPLICATED") is not True or
            r.get("source_12_old_ERA5_raw_SHA_HOLD") is not True or
            r.get("no_original_tag_ID_GPS_exact_date_climate_cell_or_individual_speeds_exported") is not True):
            raise ValueError("Source original identity, source rights, old negative science changed")
        status=r["status"]
        if status not in ("SOURCE_HOLD","SOURCE_BALANCED_3_EACH_COMPONENT_ESTIMABLE"):
            raise ValueError("Invalid original 3step balance result status")
        result=r["result"]
        if status=="SOURCE_HOLD" and result is not None:
            raise ValueError("Partial original source 3step result must be HOLD")
        if status!="SOURCE_HOLD":
            if (not r["provider_source_MD5_verified"]
                or not r["source_native_full_q_and_daylength_replayed"]
                or not r["original_source_q_qualified_tag_solar_date_set_reconstructed"]
                or not r["exact_first_ERA5_daily_physical_SHA_verified"]
                or not r["previous_original_source_six_weather_model_means_replayed"]
                or not r["original_first_three_source_DAY_NIGHT_Q_tag_bootstraps_replayed"]
                or r["source_weather_grids_HOLD"]!=0
                or r["source_weather_grids"]!=old[i]["new_weather_source_grid_count"]
                or result is None):
                raise ValueError("Original GPS/source weather/first q model not replayed")
            if (result["original_q_weather_month_tagged_animals_identical_to_first"]!=
                p["source_original_individual_tag_n"] or
                not result["original_source_daily_q_and_first_heat_bootstrap_exact_replay"] or
                not result["equal_three_new_q_NIGHT_minus_DAY_algebra_exact"] or
                set(result["source_effects"])!=set(KEYS)):
                raise ValueError("Balanced source q estimate missing algebra/previous tag subset")
            for key in KEYS:validate(result["source_effects"][key])
            e=result["source_effects"]
            for prefix in ("old_month","balanced_month"):
                if abs(e[prefix+"_NIGHT_Tmax"]["mean"]-
                       e[prefix+"_DAY_Tmax"]["mean"]-
                       e[prefix+"_Q_Tmax"]["mean"])>TOL:
                    raise ValueError("Old or balanced source native q temperature != night-day")
            if not 0<=result["original_tag_solar_dates_with_unequal_count_fraction"]<=1:
                raise ValueError("Bad original source q date opportunity inequality fraction")
        effects=result["source_effects"] if result else {}
        table.append({
            "source_key":s["source_key"],
            "original_source_DOI":s["doi"],
            "original_native_taxon":s["taxon"],
            "source_status":status,
            "original_source_independent_tag_count":(
                result["original_q_weather_month_tagged_animals_identical_to_first"] if result else None),
            "original_source_tag_solar_days_in_model_not_independent":(
                result["original_tag_solar_dates_in_admitted_models_not_independent"] if result else None),
            "original_source_fraction_solar_tag_dates_unequal_day_night_step_count":(
                result["original_tag_solar_dates_with_unequal_count_fraction"] if result else None),
            "original_solar_q_speed_temperature":effects.get("old_month_Q_Tmax"),
            "balanced_three_solar_q_speed_temperature":effects.get("balanced_month_Q_Tmax"),
            "balanced_three_DAY_speed_temperature":effects.get("balanced_month_DAY_Tmax"),
            "balanced_three_NIGHT_speed_temperature":effects.get("balanced_month_NIGHT_Tmax"),
            "equal_three_minus_original_q_temperature":effects.get("balanced_minus_old_month_Q_Tmax"),
            "source_verified_original_MD5_and_physical_ERA5_q_source":status!="SOURCE_HOLD",
            "hold_stage":r["hold_stage"],
            "hold_type":r["hold_type"],
        })
    qualified=[x for x in table if x["source_status"]=="SOURCE_BALANCED_3_EACH_COMPONENT_ESTIMABLE"]
    return {
        "version":VERSION,
        "stage":"POST_PRIOR_THREE_SOURCE_GPS_DAY_NIGHT_THERMAL_FIRST_BALANCED_INTERVAL_SENSITIVITY",
        "original_source_study_DOI_count":3,
        "all_original_three_study_source_rows_including_HOLD":table,
        "source_status_counts":dict(sorted(Counter(z["source_status"] for z in table).items())),
        "original_source_studies_with_equal_step_model":len(qualified),
        "equal_three_q_Tmax_nominal_positive_sources":sum(
            x["balanced_three_solar_q_speed_temperature"]["nominal_sign"]=="POSITIVE" for x in qualified),
        "equal_three_DAY_Tmax_nominal_negative_sources":sum(
            x["balanced_three_DAY_speed_temperature"]["nominal_sign"]=="NEGATIVE" for x in qualified),
        "equal_three_NIGHT_Tmax_nominal_positive_sources":sum(
            x["balanced_three_NIGHT_speed_temperature"]["nominal_sign"]=="POSITIVE" for x in qualified),
        "equal_three_vs_full_q_change_source_95_sign":{
            x["source_key"]:x["equal_three_minus_original_q_temperature"]["nominal_sign"]
            for x in qualified},
        "prior_original_day_night_source_git_blob_preserved":
            c["prior_three_source_day_night_first_result_blob"],
        "three_original_source_DOIs_not_verified_independent_populations":True,
        "source_data_selected_after_prior_positive_heat":True,
        "observed_GPS_sampling_robustness_not_causal_temperature_effect":True,
        "source_previous_48_irradiance_heading_signed_external_result_NOT_REPLICATED":True,
        "source_previous_12_ERA5_raw_SHA_HOLD":True,
        "no_original_GPS_animal_ID_source_exact_date_grid_or_individual_result_exported":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    files=list(a.directory.glob("source_*.json"))
    if len(files)!=3:raise ValueError("Exactly three original source results must be supplied")
    indexed={}
    for f in files:
        x=json.loads(f.read_text(encoding="utf-8"))
        idx=x.get("source_index")
        if type(idx)!=int or idx not in (0,1,2) or idx in indexed:
            raise ValueError("Duplicate original source ID or missing source data")
        indexed[idx]=x
    if len(indexed)!=3:raise ValueError("All original source ID slots required")
    res=aggregate([indexed[i] for i in range(3)])
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_THREE_SOURCE_EQUAL_INTERVAL_Q",json.dumps({
        "source_status":res["source_status_counts"],
        "q_positive_after_balancing":res["equal_three_q_Tmax_nominal_positive_sources"],
        "DAY_negative":res["equal_three_DAY_Tmax_nominal_negative_sources"],
        "NIGHT_positive":res["equal_three_NIGHT_Tmax_nominal_positive_sources"],
        "within_source_q_change":res["equal_three_vs_full_q_change_source_95_sign"]
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
