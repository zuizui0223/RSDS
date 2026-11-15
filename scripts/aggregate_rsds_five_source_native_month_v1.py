"""Aggregate only all five SHA-frozen source-native postoutcome month checks.

Do not let a source MD5 HOLD or insufficient month coverage vanish. A single
original study DOI is one source. Original source-tag bootstraps are not
cross-source or causal confidence intervals.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import json
import math
from pathlib import Path
from scripts.rsds_five_source_native_tag_month_conditioned_v1 import frozen,SOURCE_VERSION

VERSION="RSDS_FIVE_ORIGINAL_SOURCE_NATIVE_TAG_YEAR_MONTH_META_V1"


def aggregate(rows:list):
    c=frozen()
    if len(rows)!=5:raise ValueError("All five original source ID outcomes required")
    unique=set()
    output=[]
    for i,(s,r) in enumerate(zip(c["source_manifest"],rows)):
        if (r.get("version")!=SOURCE_VERSION or r.get("source_index")!=i or
            r.get("original_study_DOI")!=s["doi"] or
            r.get("original_study_id")!=str(s["study_id"]) or
            r.get("original_native_target")!=s["taxon"] or
            r.get("all_original_registered_native_taxa")!=s["registered_taxa"] or
            r.get("old_solar_direction_external_replication_still_negative") is not True or
            r.get("old_original_ERA5_raw_SHA_hold") is not True or
            r.get("no_raw_animal_GPS_tracks_ID_exact_coordinates_dates_or_tag_slopes_exported") is not True):
            raise ValueError("Native source metadata or previous negative science changed")
        if r["original_study_DOI"] in unique:raise ValueError("Duplicated original source DOI")
        unique.add(r["original_study_DOI"])
        status=r["status"]
        if status not in ("SOURCE_HOLD","MONTH_NATIVE_ESTIMABLE",
                          "MONTH_NATIVE_UNAVAILABLE_INSUFFICIENT_SOURCE_SUPPORT"):
            raise ValueError("Unknown source month status")
        if status!="SOURCE_HOLD" and not (r["baseline_full_year_first_results_replayed"] and r["original_file_MD5_bytes_verified"]):
            raise ValueError("Month analysis used source lacking MD5 or full-source replay")
        if status=="SOURCE_HOLD" and r["month_conditioned_native"] is not None:
            raise ValueError("Source-wide HOLD leaked partial month result")
        new=r["month_conditioned_native"]
        local=(new or {}).get("local_original_solar_daylength") or {}
        fixed=(new or {}).get("fixed_source_reference_latitude_daylength") or {}
        if status=="MONTH_NATIVE_ESTIMABLE":
            if (local.get("status")!="ESTIMABLE" or fixed.get("status")!="ESTIMABLE" or
                local.get("tagged_animals_support",0)<5 or
                local.get("tagged_animals_support")!=fixed.get("tagged_animals_support")):
                raise ValueError("Month native source sign/animal support mismatched")
        if status=="MONTH_NATIVE_UNAVAILABLE_INSUFFICIENT_SOURCE_SUPPORT" and local.get("status")!="INSUFFICIENT_NATIVE_TAG_SUPPORT":
            raise ValueError("Source month HOLD not honest")
        for effect in (local,fixed):
            if effect.get("status")=="ESTIMABLE":
                ci=effect.get("original_tag_bootstrap95")
                v=effect.get("mean_slope_per_hour")
                if (not isinstance(ci,list) or len(ci)!=2 or not math.isfinite(v) or
                    not all(math.isfinite(x) for x in ci) or ci[0]>ci[1]):
                    raise ValueError("Ill-formed month conditional source bootstrap")
                sign="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
                if effect.get("nominal_interval_sign")!=sign:
                    raise ValueError("Changed nominal source bootstrap sign")
            elif effect.get("mean_slope_per_hour") is not None:
                raise ValueError("Unsupported month slope substituted with a number")
        oldq=s["baseline"]["q"]
        olddl=s["baseline"]["daylength"]
        new_ci=local.get("nominal_interval_sign","UNAVAILABLE")
        old_ci=olddl["nominal_source_CI_sign"]
        ref_ci=fixed.get("nominal_interval_sign","UNAVAILABLE")
        output.append({
            "source_index":i,
            "source_key":s["source_key"],
            "original_DOI":s["doi"],
            "source_native_species":s["taxon"],
            "original_registered_source_species":s["registered_taxa"],
            "source_gps_status":status,
            "previous_full_source_native_q":oldq,
            "previous_full_source_native_daylength":olddl,
            "source_first_result_baseline_replayed":r["baseline_full_year_first_results_replayed"],
            "source_original_MD5_verified":r["original_file_MD5_bytes_verified"],
            "month_native":new,
            "within_tag_year_month_local_vs_old_full_slope":
                ("OPPOSITE_NOMINAL_SOURCE_SIGN" if old_ci in ("POSITIVE","NEGATIVE") and
                 new_ci in ("POSITIVE","NEGATIVE") and old_ci!=new_ci else
                 "SAME_NOMINAL_SOURCE_SIGN" if old_ci in ("POSITIVE","NEGATIVE") and
                 new_ci==old_ci else
                 "UNRESOLVED_OR_UNAVAILABLE"),
            "fixed_lat_month_robustness_vs_local":
                ("BOTH_NOMINAL_MATCHING" if new_ci in ("POSITIVE","NEGATIVE") and ref_ci==new_ci else
                 "NOMINAL_DIFFERENT" if new_ci in ("POSITIVE","NEGATIVE") and ref_ci in ("POSITIVE","NEGATIVE") and ref_ci!=new_ci else
                 "UNRESOLVED_OR_UNAVAILABLE"),
            "source_hold_stage":r["original_source_hold_stage"],
            "source_hold_type":r["original_source_hold_type"],
            "previous_original_GPS_exposure_state":s["source_file_prior_access"]
        })
    by=defaultdict(list)
    for r in output:by[r["source_native_species"]].append(r)
    species_synth={}
    for sp,records in sorted(by.items()):
        local=[(x["source_key"],x["month_native"]["local_original_solar_daylength"]["nominal_interval_sign"])
               for x in records if x["source_gps_status"]=="MONTH_NATIVE_ESTIMABLE"]
        strong=[(a,b) for a,b in local if b in ("POSITIVE","NEGATIVE")]
        species_synth[sp]={
            "all_source_studies":len(records),
            "months_estimated_source_studies":len(local),
            "month_native_source_signs":{k:v for k,v in local},
            "month_native_nominal_positive_source_count":sum(v=="POSITIVE" for k,v in local),
            "month_native_nominal_negative_source_count":sum(v=="NEGATIVE" for k,v in local),
            "month_native_nominal_includes_zero_source_count":sum(v=="INCLUDES_ZERO" for k,v in local),
            "both_month_source_CI_signs_supported_and_opposite":any(v=="POSITIVE" for k,v in strong) and any(v=="NEGATIVE" for k,v in strong),
            "distinct_tagged_populations_confirmed":False,
            "cross_source_interaction_difference_bootstrap_estimated":False,
        }
    return {
        "version":VERSION,
        "stage":"POSTHOC_FIVE_SOURCE_TAG_YEAR_MONTH_CONDITIONED_OBSERVATIONAL_SENSITIVITY",
        "five_frozen_original_field_source_DOI":5,
        "two_previously_targeted_species":2,
        "source_status_count":dict(sorted(Counter(x["source_gps_status"] for x in output).items())),
        "all_five_frozen_original_source_native_analyses_including_HOLD":output,
        "per_species_original_source_native_month_sign_descriptions":species_synth,
        "source_tag_year_month_demeaning_not_a_season_treatment":True,
        "two_species_and_month_tests_selected_post_previous_outcomes":True,
        "not_prospective_global_macroecology":True,
        "between_source_population_and_tag_disjointness_unverified":True,
        "no_source_to_source_causal_daylength_comparison":True,
        "old_48_solar_direction_external_nonreplication_UNCHANGED":True,
        "old_12_ERA5_raw_SHA_HOLD_UNCHANGED":True,
        "original_GPS_positions_dates_IDs_and_individual_month_slopes_not_exported":True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--directory",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    arg=parser.parse_args()
    files=list(arg.directory.glob("source_*.json"))
    if len(files)!=5:raise ValueError("All five original source results mandatory")
    by={}
    for f in files:
        obj=json.loads(f.read_text(encoding="utf-8"))
        i=obj.get("source_index")
        if not isinstance(i,int) or not 0<=i<5 or i in by:raise ValueError("Bad original source index")
        by[i]=obj
    if set(by)!=set(range(5)):raise ValueError("Five original source IDs not complete")
    result=aggregate([by[i] for i in range(5)])
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_FIVE_SOURCE_NATIVE_MONTH_AGG",json.dumps({
        "source_status":result["source_status_count"],
        "species":result["per_species_original_source_native_month_sign_descriptions"],
        "no_original_GPS_export":True,
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
