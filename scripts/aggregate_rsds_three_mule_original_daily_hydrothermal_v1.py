"""Three original mule deer source-study joint daily climate q associations.

No weather result from a source with a provider hold may be promoted to a
full-source ecological claim; keep all 3 frozen source IDs and historical HOLDS.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path
from scripts.rsds_three_mule_original_daily_hydrothermal_v1 import frozen,SOURCE_VERSION

ENDPOINTS=(
"unadjusted_daylight_on_climate_days",
"daylight_given_tmax",
"daylight_given_antecedent_water",
"daylight_given_tmax_and_water",
"tmax_given_daylight_and_water",
"water_given_daylight_and_tmax",
)

def categorize(effect):
    if effect is None or effect.get("mean") is None:return "UNAVAILABLE"
    ci=effect.get("animal_bootstrap_95")
    if (not isinstance(ci,list) or len(ci)!=2 or
        not all(isinstance(x,(int,float)) and math.isfinite(x) for x in ci) or
        ci[0]>ci[1] or not isinstance(effect["mean"],(int,float)) or
        not math.isfinite(effect["mean"]) or not isinstance(effect.get("admitted_animals"),int) or
        effect["admitted_animals"]<5):
        raise ValueError("Wrong original within-source tag bootstrap support")
    return "POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"


def aggregate(results):
    c,sources=frozen()
    if len(results)!=3:raise ValueError("All 3 original mule source ID slots including HOLD required")
    seen=set();rows=[]
    for i,(s,r) in enumerate(zip(sources,results)):
        if (r.get("version")!=SOURCE_VERSION or r.get("source_index")!=i or
            r.get("source_key")!=s["source_key"] or
            r.get("original_study_DOI")!=s["doi"] or
            r.get("original_study_id")!=str(s["study_id"]) or
            r.get("original_native_taxon")!="Odocoileus hemionus" or
            r.get("registered_original_native_taxa")!=s["registered_taxa"] or
            r.get("original_48_sunlight_heading_external_NOT_REPLICATED") is not True or
            r.get("original_12_ERA5_raw_SHA_HOLD_unchanged") is not True or
            r.get("no_original_tag_ID_exact_dates_coordinates_GPS_weather_events_exported") is not True):
            raise ValueError("Source registered ID/taxon/privacy or old science changed")
        if s["doi"] in seen:raise ValueError("Duplicate source study DOI")
        seen.add(s["doi"])
        status=r["source_status"]
        if status not in ("ORIGINAL_SOURCE_ERA5_FULL_GRID_COVERAGE_EVALUATED",
                          "SOURCE_PARTIAL_ERA5_GRID_HOLD_DESCRIPTIVE_ONLY",
                          "SOURCE_HOLD"):
            raise ValueError("Unrecognized daily climate source status")
        if status!="SOURCE_HOLD" and not (r["original_source_file_MD5_verified"] and
                r["original_source_native_q_and_daylength_replayed_before_weather"]):
            raise ValueError("Original q/source must replay before weather")
        if status=="SOURCE_HOLD" and r["all_source_native_weather_adjusted_effects"] is not None:
            raise ValueError("Source-wide HOLD leaked source partial biology")
        if status=="ORIGINAL_SOURCE_ERA5_FULL_GRID_COVERAGE_EVALUATED" and (
            r["source_ERA5_grid_HOLD_count"]!=0 or
            r["all_source_native_weather_adjusted_effects"]["all_weather_grid_inputs_retrieved"] is not True):
            raise ValueError("Zero weather HOLD required for complete source")
        d=r["all_source_native_weather_adjusted_effects"]
        endpoint={}
        if d:
            for name in ENDPOINTS:
                if name not in d:raise ValueError("Missing prespecified historical Wyoming weather outcome")
                endpoint[name]={
                    "original_source_tagged_animals":d[name]["admitted_animals"],
                    "source_original_animal_years":d[name]["admitted_animal_years"],
                    "mean":d[name]["mean"],
                    "original_tag_bootstrap95":d[name]["animal_bootstrap_95"],
                    "status":d[name]["status"],
                    "nominal_sign":categorize(d[name])
                }
        rows.append({
            "original_source_key":s["source_key"],
            "original_study_DOI":s["doi"],
            "source_native_taxon":s["taxon"],
            "source_gps_status":status,
            "native_GPS_source_MD5_replayed":r["original_source_file_MD5_verified"],
            "original_native_q_first_result_replayed":r["original_source_native_q_and_daylength_replayed_before_weather"],
            "first_original_source_q":s["baseline"]["q"],
            "first_original_source_full_year_daylength":s["baseline"]["daylength"],
            "weather_model_years_new_version_not_historically_equivalent":True,
            "current_provider_daily_physical_fingerprint":r["new_weather_physical_SHA256_not_coordinates"],
            "new_weather_source_grid_count":r["source_ERA5_cell_count"],
            "new_weather_source_grid_HOLD_count":r["source_ERA5_grid_HOLD_count"],
            "original_qualified_tag_year_count_before_weather":r["original_q_eligible_tag_years"],
            "model_eligible_qualified_tag_years":d["animal_years_with_at_least_ten_climate_complete_source_days_and_solar_span"] if d else None,
            "weather_complete_qualified_animal_days":d["weather_complete_original_animal_days"] if d else None,
            "source_weather_adjusted_endpoints":endpoint,
            "source_hold_stage":r["source_hold_stage"],
            "source_hold_type":r["source_hold_type"],
        })
    eligible=[z for z in rows if z["source_gps_status"]=="ORIGINAL_SOURCE_ERA5_FULL_GRID_COVERAGE_EVALUATED"]
    primary={z["original_source_key"]:z["source_weather_adjusted_endpoints"].get("daylight_given_tmax_and_water",{}).get("nominal_sign")
             for z in eligible}
    contexts={z["original_source_key"]:z["source_weather_adjusted_endpoints"].get("tmax_given_daylight_and_water",{}).get("nominal_sign")
             for z in eligible}
    return {
        "version":"RSDS_THREE_MULE_ORIGINAL_DAILY_HYDROTHERMAL_SOURCE_FIRST_META_V1",
        "stage":"FIRST_POSTDISCOVERY_ORIGINAL_SOURCE_Q_WATER_TMAX_EXTERNAL_DESCRIPTIVE",
        "source_study_DOI_denominator":3,
        "one_target_species":"Odocoileus hemionus",
        "original_source_status_counts":dict(sorted(Counter(z["source_gps_status"] for z in rows).items())),
        "all_three_original_source_rows_including_HOLD":rows,
        "full_weather_coverage_source_studies":len(eligible),
        "primary_daylength_joint_adjusted_nominal_CI_signs_by_fully_covered_source":primary,
        "Tmax_given_daylength_water_nominal_CI_signs_by_fully_covered_source":contexts,
        "original_source_weather_version_new_not_equivalent_to_historic_ERA5_v0":True,
        "historical_12_ERA5_raw_SHA_source_HOLD_UNCHANGED":True,
        "external_48_source_sunlight_heading_direction_NOT_REPLICATED_UNCHANGED":True,
        "not_independent_sample_of_world_mule_deer_populations":True,
        "no_meteorological_causality_or_vegetation_migration_state_measured":True,
        "no_source_to_source_tagged_population_overlap_test":True,
        "postoutcome_species_and_weather_model_selection":True,
        "missing_weather_grids_held_not_imputed":True,
        "source_group_only_no_original_tag_GPS_weather_or_dates_exported":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    fs=list(args.directory.glob("source_*.json"))
    if len(fs)!=3:raise ValueError("Exactly three original source result files including HOLD")
    indexed={}
    for f in fs:
        d=json.loads(f.read_text(encoding="utf-8"))
        i=d.get("source_index")
        if not isinstance(i,int) or i not in range(3) or i in indexed:
            raise ValueError("Duplicate or missing original source index")
        indexed[i]=d
    if len(indexed)!=3:raise ValueError("Three source original denominator incomplete")
    z=aggregate([indexed[i] for i in range(3)])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(z,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_THREE_MULE_WEATHER_AGG",json.dumps({
        "study_source_IDs":3,
        "status":z["original_source_status_counts"],
        "daylength_given_Tmax_water":z["primary_daylength_joint_adjusted_nominal_CI_signs_by_fully_covered_source"],
        "Tmax_given_L_water":z["Tmax_given_daylength_water_nominal_CI_signs_by_fully_covered_source"],
        "not_a_world_meta":True
    },sort_keys=True),flush=True)
if __name__=="__main__":main()
