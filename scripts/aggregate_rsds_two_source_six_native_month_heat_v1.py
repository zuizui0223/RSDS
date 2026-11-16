"""Two DOI × all original registered six native ungulate source study ecology."""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path
from scripts.rsds_two_source_all_native_ungulate_month_heat_v1 import (
    frozen,SOURCE_VERSION,TOL
)
VERSION="RSDS_TWO_STUDY_SIX_UG_NATIVE_WITHIN_MONTH_HEAT_FIRST_SOURCE_META_V1"

def validated(x):
    if x is None:return "UNAVAILABLE"
    if x.get("status") not in ("ESTIMABLE","INSUFFICIENT_ORIGINAL_TAGS"):
        raise ValueError("Nonfrozen original source native tagged-animal bootstrap status")
    if x["status"]!="ESTIMABLE":
        if x.get("estimate") is not None or x.get("tag_bootstrap95") is not None:
            raise ValueError("Unavailable original native source q must not become a zero response")
        return "UNAVAILABLE"
    ci=x["tag_bootstrap95"]
    if (not isinstance(x["original_tag_n"],int) or x["original_tag_n"]<5 or
        not isinstance(ci,list) or len(ci)!=2 or
        not all(isinstance(z,(int,float)) and math.isfinite(z) for z in ci) or
        not ci[0]<=ci[1] or not isinstance(x.get("estimate"),(int,float)) or
        not math.isfinite(x["estimate"])):
        raise ValueError("Native original source q/DAY/NIGHT bootstrap invalid")
    sign="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
    if x["nominal_CI_sign"]!=sign:
        raise ValueError("Altered original tagged-animal nominal CI class")
    return sign

def aggregate(source):
    c,manifest,history=frozen()
    if len(source)!=2:raise ValueError("Both original source native multi-ungulate study IDs are mandatory")
    rows=[];source_states=[];source_shared=[]
    for i,(s,r) in enumerate(zip(manifest,source)):
        if (r.get("version")!=SOURCE_VERSION or
            r.get("source_index")!=i or
            r.get("original_study_DOI")!=s["doi"] or
            r.get("original_repository_study_id")!=s["id"] or
            r.get("registered_native_taxa")!=s["taxa"] or
            r.get("no_original_animal_ID_tracked_coordinates_exact_date_GPS_weather_grid_or_individual_effect_exported") is not True or
            r.get("previous_48_source_irradiance_heading_direction_NOT_REPLICATED") is not True or
            r.get("previous_12_source_ERA5_raw_SHA_HOLD") is not True or
            r.get("source_new_weather_is_not_previous_20261023_or_20261116_physical_replay") is not True):
            raise ValueError("Original source/native privacy, DOI or old ecological science changed")
        status=r["source_status"]
        if status not in ("SOURCE_ALL_NATIVE_ERA5_MONTH_MODELS_EVALUATED","SOURCE_HOLD"):
            raise ValueError("Unknown original full source/registered native weather status")
        state={"source":s["doi"],"native_count":len(s["taxa"]),"status":status,
               "original_source_weather_grids":r["source_ERA5_grid_count"],
               "weather_HOLDs":r["source_ERA5_grid_HOLDs"],
               "new_source_ERA5_fingerprint":r["source_ERA5_daily_physical_SHA256"]}
        source_states.append(state)
        if status=="SOURCE_HOLD":
            if r["source_native_group_results"]:
                raise ValueError("Source original provider HOLD leaked partial registered-native outcome")
            rows.extend({
                "source_original_DOI":s["doi"],
                "native_taxon":sp,
                "status":"SOURCE_HOLD",
                "original_tag_n":None,
                "q_Tmax_source_CI_sign":"UNAVAILABLE",
                "day_Tmax_source_CI_sign":"UNAVAILABLE",
                "night_Tmax_source_CI_sign":"UNAVAILABLE",
                "original_q_Tmax_month":None,
                "original_day_Tmax_month":None,
                "original_night_Tmax_month":None,
            } for sp in s["taxa"])
            continue
        if (r["source_original_GPS_MD5_verified"] is not True or
            r["all_registered_native_q_and_solar_group_first_results_replayed"] is not True or
            r["original_q_NIGHT_minus_DAY_per_tag_date_reconstructed"] is not True or
            r["new_contemporaneous_source_native_ERA5_full_grid_coverage"] is not True or
            r["source_ERA5_grid_HOLDs"]!=0 or
            not isinstance(r["source_ERA5_daily_physical_SHA256"],str) or
            len(r["source_ERA5_daily_physical_SHA256"])!=64 or
            len(r["source_native_group_results"])!=len(s["taxa"]) or
            [z["original_native_taxon"] for z in r["source_native_group_results"]]!=s["taxa"]):
            raise ValueError("Original source native GPS/provider complete weather not verified")
        source_shared.append({
            "original_source_DOI":s["doi"],
            "mule_and_bighorn_original_shared_grid_yearmonth_count":
                r["original_source_mule_and_bighorn_shared_source_grid_yearmonth_count"],
            "shared_grid_yearmonth_both_taxa_at_least_three_original_tagged_animals":
                r["original_source_mule_and_bighorn_shared_grid_yearmonth_with_3_tags_each"],
        })
        for taxon,z in zip(s["taxa"],r["source_native_group_results"]):
            if z["source_native_status"] not in ("SOURCE_NATIVE_MONTH_COMPONENTS_ESTIMABLE","SOURCE_NATIVE_MONTH_INSUFFICIENT_TAGS"):
                raise ValueError("Unknown species native within-original-month q result")
            prior=history[s["doi"]][taxon]
            if (z["first_full_original_q_tagged_animals"]!=prior["q"]["source_tags"] or
                z["first_original_daylength_tagged_animals"]!=prior["daylength"]["source_tags"] or
                abs(z["first_full_original_q"]-prior["q"]["mean"])>TOL or
                abs(z["first_original_daylength_slope"]-prior["daylength"]["mean"])>TOL):
                raise ValueError("Previously pinned species-specific first source q/daylength replay not preserved")
            q=z["month_native_q_Tmax_given_daylength_antecedent_water"]
            day=z["month_native_DAY_Tmax_given_daylength_antecedent_water"]
            night=z["month_native_NIGHT_Tmax_given_daylength_antecedent_water"]
            qsign=validated(q);dsign=validated(day);nsign=validated(night)
            n=z["month_joint_original_tag_n"]
            if z["source_native_status"]=="SOURCE_NATIVE_MONTH_COMPONENTS_ESTIMABLE":
                if (day is None or night is None or
                    q["original_tag_n"]!=n or day["original_tag_n"]!=n or night["original_tag_n"]!=n or
                    abs((night["estimate"]-day["estimate"])-q["estimate"])>TOL or
                    z["original_native_q_NIGHT_minus_DAY_model_and_bootstrap_algebra_verified"] is not True or
                    z["first_source_native_q_Tmax_new_ERA5_tag_bootstrap_replayed"] is not True):
                    raise ValueError("DAY NIGHT q original source algebra/animal support not verified")
            elif qsign!="UNAVAILABLE" or day is not None or night is not None:
                raise ValueError("Unsupported per-source native tag month q changed into measured effect")
            rows.append({
                "source_original_DOI":s["doi"],
                "native_taxon":taxon,
                "status":z["source_native_status"],
                "original_tag_n":n,
                "original_full_source_q_baseline":z["first_full_original_q"],
                "original_full_source_q_tag_count":z["first_full_original_q_tagged_animals"],
                "original_full_source_daylength_slope_baseline":z["first_original_daylength_slope"],
                "original_q_Tmax_month":q,
                "original_day_Tmax_month":day,
                "original_night_Tmax_month":night,
                "q_Tmax_source_CI_sign":qsign,
                "day_Tmax_source_CI_sign":dsign,
                "night_Tmax_source_CI_sign":nsign,
            })
    if len(rows)!=6 or [r["native_taxon"] for r in rows]!=[
        "Antilocapra americana","Cervus elaphus","Odocoileus hemionus",
        "Ovis canadensis","Odocoileus hemionus","Ovis canadensis"]:
        raise ValueError("All originally registered six source-native ungulates required")
    complete=[x for x in rows if x["status"]=="SOURCE_NATIVE_MONTH_COMPONENTS_ESTIMABLE"]
    return {
        "version":VERSION,
        "stage":"POSTDISCOVERY_TWO_ORIGINAL_SOURCES_SIX_NATIVE_TAXA_SAME_TAG_YEAR_MONTH_DAY_NIGHT_TMAX",
        "original_source_studies":2,
        "original_native_taxon_study_rows":6,
        "original_source_statuses":source_states,
        "all_original_six_source_native_rows_including_HOLD":rows,
        "source_original_mule_vs_bighorn_shared_grid_yearmonth_opportunity_without_coordinates":source_shared,
        "source_native_q_Tmax_sign_counts":dict(sorted(Counter(x["q_Tmax_source_CI_sign"] for x in rows).items())),
        "source_native_day_Tmax_sign_counts":dict(sorted(Counter(x["day_Tmax_source_CI_sign"] for x in rows).items())),
        "source_native_night_Tmax_sign_counts":dict(sorted(Counter(x["night_Tmax_source_CI_sign"] for x in rows).items())),
        "original_source_complete_tag_month_native_rows":len(complete),
        "new_20261117_ERA5_snapshot_not_old_20261023_or_20261116_raw_physical_equivalence":True,
        "original_sources_and_native_taxa_selected_after_old_q_outcomes":True,
        "study_DOI_not_equal_taxon_replication_or_same_cell_month_experiment":True,
        "two_studies_not_six_independent_ecosystems":True,
        "no_formal_species_by_temperature_interaction_test":True,
        "prior_48_solar_irradiance_heading_direction_not_replicated":True,
        "prior_12_ERA5_original_raw_SHA_HOLD":True,
        "no_raw_GPS_weather_animal_ID_individual_effects_or_precise_dates_exported":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    fs=list(a.directory.glob("source_*.json"))
    if len(fs)!=2:raise ValueError("All two original study source first outputs including HOLD mandatory")
    d={}
    for f in fs:
        row=json.loads(f.read_text(encoding="utf-8"));i=row.get("source_index")
        if type(i)!=int or i not in (0,1) or i in d:
            raise ValueError("Source original DOI index duplicated/missing")
        d[i]=row
    if set(d)!={0,1}:raise ValueError("Two source DOI slots required")
    res=aggregate([d[i] for i in (0,1)])
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("RSDS_ORIGINAL_SIX_NATIVE_MONTH_TMAX",json.dumps({
        "source_studies":2,"all_native":6,
        "source":res["original_source_statuses"],
        "q":res["source_native_q_Tmax_sign_counts"],
        "day":res["source_native_day_Tmax_sign_counts"],
        "night":res["source_native_night_Tmax_sign_counts"],
        "within_source_mule_bighorn_shared_gridmonth":res["source_original_mule_vs_bighorn_shared_grid_yearmonth_opportunity_without_coordinates"],
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
