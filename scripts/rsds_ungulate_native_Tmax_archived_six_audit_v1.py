"""Validate ALL archived 2026-10-23 native ungulate q~Tmax source receipts.

No new GPS, new ERA5 or animal identifiers. Source-level original tagged-animal
nominal CIs are within SOURCE, NOT source-independent or species-interaction
meta CIs. Read only fixed original Utah 4-source-native contract/previous run
receipt and Wyoming 2-native source result. Never infer true activity time.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned

PATH=HERE/"research/RSDS_UTAH_WYOMING_SIX_UNGULATE_NATIVE_HEAT_Q_REFERENCE_20261116.json"
VERSION="RSDS_UTAH_WYOMING_SIX_SOURCE_NATIVE_UNGULATE_Q_TMAX_ARCHIVED_FIRST_V1"
UTAH_TAXA=("Antilocapra americana","Cervus elaphus","Odocoileus hemionus","Ovis canadensis")
WY_TAXA=("Odocoileus hemionus","Ovis canadensis")
SIGN={"POSITIVE","NEGATIVE","INCLUDES_ZERO","UNAVAILABLE"}

def audit():
    a=json.loads(PATH.read_text(encoding="utf-8"))
    if (a.get("version")!=VERSION
        or a.get("source_count")!=2 or a.get("source_native_rows")!=6
        or a.get("source_original_weather_ERA5_version", "").find("NOT the newly")<0
        or a.get("old_48_irradiance_solar_heading_signed_external_test")!="NOT_REPLICATED"
        or a.get("old_12_ERA5_raw_SHA_source_cases")!="HOLD"
        or a.get("raw_animal_GPS_tags_coordinates_exact_dates_weather_grid_not_exported") is not True):
        raise ValueError("Frozen six original source-native science/ERA5 caveat changed")
    u=a["source_provenance"]["utah"];w=a["source_provenance"]["wyoming"]
    if (u["original_run_id"]!=39878716871
        or u["original_job_id"]!=120170619463
        or u["original_signed_artifact_id"]!=12469413417
        or u["original_artifact_digest"]!="sha256:ab29be110408b8f028b2f1fc925741f6c7935c95c10a265e6fa9587a146493f0"
        or u["original_daily_ERA5_grid_count"]!=28 or u["original_source_weather_grid_hold_count"]!=0
        or w["first_original_source_run_id"]!=39878174016
        or w["original_daily_ERA5_grid_count"]!=6 or w["original_source_weather_grid_hold_count"]!=0):
        raise ValueError("Archived first original GPS/ERA5 source artifact provenance altered")
    spec=pinned(u["original_source_contract_path"],u["contract_git_blob"])
    original_wy=pinned(w["original_source_result_path"],w["result_git_blob"])
    if ([x["taxon"] for x in spec["source_reference_results"]]!=list(UTAH_TAXA)
        or [x["taxon"] for x in original_wy["source_group_results"]]!=list(WY_TAXA)
        or original_wy["status"]["previous_GPS_source_replay"] is not True
        or original_wy["status"]["weather_holds"]!=0
        or a["source_provenance"]["utah"]["original_q_and_daylength_results_replayed"] is not True):
        raise ValueError("Original public source-native cohort changed")
    rows=a["original_tmax_q_native_rows"]
    if (len(rows)!=6
        or [x["source"] for x in rows]!=["utah_2019_2020"]*4+["wyoming_2019_2020"]*2
        or [x["native_taxon"] for x in rows]!=list(UTAH_TAXA+WY_TAXA)):
        raise ValueError("One of six original taxon/source native rows is missing")
    taxonmap=defaultdict(list)
    for i,x in enumerate(rows):
        prov=u if i<4 else w
        if (x["study_doi"]!=prov["doi"]
            or not isinstance(x["tag_n"],int) or x["tag_n"]<5
            or not isinstance(x["animal_years"],int) or x["animal_years"]<x["tag_n"]
            or not isinstance(x["tmax_q_per_degree_C"],(float,int))
            or not math.isfinite(x["tmax_q_per_degree_C"])
            or not isinstance(x["original_tag_bootstrap_95"],list)
            or len(x["original_tag_bootstrap_95"])!=2
            or any(not isinstance(v,(int,float)) or not math.isfinite(v)
                   for v in x["original_tag_bootstrap_95"])
            or x["original_tag_bootstrap_95"][0]>x["original_tag_bootstrap_95"][1]
            or x["nominal_sign"] not in SIGN
            or not isinstance(x["original_source_daylevel_ERA5_complete_q_days"],int)
            or x["original_source_daylevel_ERA5_complete_q_days"]<x["animal_years"]):
            raise ValueError("Original source tagged-animal bootstrap or calendar support bad")
        lo,hi=x["original_tag_bootstrap_95"]
        sign="POSITIVE" if lo>0 else "NEGATIVE" if hi<0 else "INCLUDES_ZERO"
        if sign!=x["nominal_sign"]:
            raise ValueError("Original nominal 95 source sign reclassified")
        if i>=4:
            p=original_wy["source_group_results"][i-4]
            effect=p["tmax_given_daylight_and_water"]
            if (x["tag_n"]!=effect["admitted_animals"]
                or x["animal_years"]!=effect["admitted_animal_years"]
                or x["tmax_q_per_degree_C"]!=effect["mean"]
                or x["original_tag_bootstrap_95"]!=effect["animal_bootstrap_95"]):
                raise ValueError("Historical Wyoming original source weather result changed")
        else:
            if x["native_taxon"]!=spec["source_reference_results"][i]["taxon"]:
                raise ValueError("Historic Utah native source first q scope changed")
        taxonmap[x["native_taxon"]].append(x)
    counts=dict(sorted(Counter(x["nominal_sign"] for x in rows).items()))
    if counts!={"INCLUDES_ZERO":2,"NEGATIVE":1,"POSITIVE":3}:
        raise ValueError("Original source-native q~Tmax mixed 6-sign pattern changed")
    return {
        "version":"RSDS_SIX_ARCHIVED_UNGULATE_NATIVE_SOURCE_TMAX_AUDIT_V1",
        "distinct_original_field_study_DOIs":2,
        "original_source_native_taxon_study_rows":6,
        "source_native_q_Tmax_nominal95_CI_sign_counts":counts,
        "same_taxon_twice_tested_species":["Odocoileus hemionus","Ovis canadensis"],
        "within_Utah_positive_species":["Cervus elaphus","Odocoileus hemionus"],
        "within_Utah_negative_species":["Ovis canadensis"],
        "within_Utah_undetermined_species":["Antilocapra americana"],
        "same_source_does_not_establish_same_grid_month_animal_cohort":True,
        "source_bootstrap_not_interaction_CI_or_world_meta_test":True,
        "historical_ERA5_physical_version_not_equal_to_new_three_mule_20261116":True,
        "original_48_source_sunlight_heading_external_NOT_REPLICATED":True,
        "original_12_ERA5_raw_SHA_HOLD":True,
        "original_GPS_and_weather_event_rows_reopened":0,
        "source_archived_four_native_original_utah_signed_artifact_digest":u["original_artifact_digest"],
        "original_DIL_gps_q_is_not_true_activity_duration":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    result=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_SIX_UNGULATE_NATIVE_HEAT",json.dumps(result,sort_keys=True),flush=True)

if __name__=="__main__":main()
