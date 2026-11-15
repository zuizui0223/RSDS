"""Aggregate first wave1 original GPS source-effects without pseudo-replication.

The deterministic SHA-priority sample was balanced 8 single/8 multi from
metadata strata of 109/11 respectively. The *finite metadata frame*
inclusion probabilities are 8/109 and 8/11 only under uniform hashing.
Do NOT treat 16 equal-weight field DOI as a representative sample of 120
ecologically eligible populations or claim design inference with outcome HOLDs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from statistics import mean

from scripts.rsds_437_wave1_original_file_metadata_v1 import freeze

N={"single_native":109,"multi_native":11}
W={k:v/120 for k,v in N.items()}

def meta_snapshot(first):
    if len(first)!=16:
        raise ValueError("All sixteen fixed source slots must be supplied")
    old,selection=freeze()
    if sum(N.values())!=old["meta_stage1_candidate_count"]:
        raise ValueError("Source strata population denominator changed")
    statuses={}
    all_rows=[]
    all_study_ids=set()
    summaries=[]
    for i,(src,result) in enumerate(zip(selection,first)):
        if (result.get("index")!=i
            or result.get("original_ordinal")!=src["ordinal"]
            or result.get("source_package_uuid")!=src["package_uuid"]
            or result.get("original_repository_study_id")!=src["original_repository_study_id"]
            or result.get("source_DOI_label")!=src["source_DOI_label"]
            or result.get("registered_native_taxa")!=src["scientific_names"]
            or result.get("source_selection_stratum")!=src["selection_stratum"]
            or len(result.get("native_effects") or [])!=len(src["scientific_names"])
            or result.get("old_48_signed_solar_shortwave_external_NONREPLICATED_UNCHANGED")!=
              "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
            or result.get("old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED") is not True
            or result.get("no_private_GPS_or_individual_data_in_output") is not True):
            raise ValueError("Source selection or old first scientific endpoint reclassified")
        if src["original_repository_study_id"] in all_study_ids:
            raise ValueError("Repeat original source study ID among independent units")
        all_study_ids.add(src["original_repository_study_id"])
        expected=src["scientific_names"]
        if [x["native_taxon"] for x in result["native_effects"]]!=expected:
            raise ValueError("One of the registered original native taxa was selectively removed")
        allowed=("PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD",
                 "ORIGINAL_SOURCE_GPS_HOLD",
                 "ORIGINAL_SOURCE_GPS_Q_UNAVAILABLE_NOT_ZERO",
                 "ORIGINAL_SOURCE_GPS_Q_ESTIMABLE")
        if result["source_status"] not in allowed:
            raise ValueError("Unknown source original first outcome/HOLD status")
        if ((result["source_status"]=="PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD") !=
            (i in (0,1,13))):
            raise ValueError("Frozen 3 original source metadata HOLDS were changed")
        if (result["source_status"].startswith("ORIGINAL_SOURCE_GPS_Q")
            and not result["original_provider_MD5_validated"]):
            raise ValueError("Original native outcomes without MD5 verified")
        if (result["source_status"]=="ORIGINAL_SOURCE_GPS_HOLD"
            and any(x["q"] is not None for x in result["native_effects"])):
            raise ValueError("Source-wide first outcome failure must void partial native effects")
        if (result["source_status"]=="ORIGINAL_SOURCE_GPS_Q_UNAVAILABLE_NOT_ZERO"
            and any(x["q"] is not None for x in result["native_effects"])):
            raise ValueError("Unsupported source q cannot contain positive effect")
        q=[];day=[]
        for v in result["native_effects"]:
            for field in ("q","daylength"):
                f=v[field]
                if f is None:
                    continue
                if (f["n_individuals"]<5
                    or not math.isfinite(f["mean"])
                    or len(f["animal_bootstrap_95_percentile"])!=2
                    or not all(math.isfinite(x) for x in f["animal_bootstrap_95_percentile"])
                    or f["animal_bootstrap_95_percentile"][0]>
                       f["animal_bootstrap_95_percentile"][1]):
                    raise ValueError("Invalid native original tagged animal group effect")
            if v["q"] is not None:
                q.append(v["q"]["mean"])
            if v["daylength"] is not None:
                day.append(v["daylength"]["mean"])
        # A complete source study must include ALL previously registered native taxa.
        # Partial source effect is recorded but cannot serve as the full study response.
        complete_q=len(q)==len(expected)
        complete_day=len(day)==len(expected)
        if (result["source_status"]=="ORIGINAL_SOURCE_GPS_Q_ESTIMABLE" and len(q)==0):
            raise ValueError("Declared source q estimable with zero taxa")
        status=result["source_status"]
        statuses[status]=statuses.get(status,0)+1
        calendar=result["native_common_calendar_covariance"]
        if len(expected)>1 and calendar["status"]=="SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED":
            K=len(expected);V=calendar["joint_covariance_modes"]
            if (calendar["original_native_taxa"]!=expected or
                len(calendar["vector_order"])!=2*K or
                any(len(v["original_source_native_covariance_matrix"])!=2*K
                    for v in V.values())):
                raise ValueError("Joint native source covariance dimension/taxa drift")
            for v in V.values():
                if v["replicates"]!=1000:
                    raise ValueError("Original joint source covariance draw count changed")
                M=v["original_source_native_covariance_matrix"]
                if any(len(row)!=2*K for row in M):
                    raise ValueError("Partial source covariance block")
                if any(M[x][x]<0 or not math.isfinite(M[x][x]) for x in range(2*K)):
                    raise ValueError("Invalid covariance diagonal")
                if any(abs(M[x][y]-M[y][x])>1e-10
                    for x in range(2*K) for y in range(2*K)):
                    raise ValueError("Asymmetric source covariance")
        info={
            "original_study_DOI":src["source_DOI_label"],
            "original_repository_study_id":src["original_repository_study_id"],
            "stratum":src["selection_stratum"],
            "source_status":status,
            "total_registered_native":len(expected),
            "native_q_supported":len(q),
            "native_daylength_supported":len(day),
            "registered_source_complete_q":complete_q,
            "registered_source_complete_daylength":complete_day,
            "source_q_equal_native_if_complete":mean(q) if complete_q else None,
            "source_daylength_equal_native_if_complete":mean(day) if complete_day else None,
            "native_covariance_status":calendar["status"],
            "source_MD5_verified":result["original_provider_MD5_validated"],
            "original_source_animal_GPS_rows_opened":result["observed_original_source_event_rows_if_opened"],
            "source_hold_stage":result["source_hold_stage"],
            "source_hold_reason":result["source_hold_reason"]
        }
        summaries.append(info)
        for v in result["native_effects"]:
            all_rows.append({
              "source_study_DOI":src["source_DOI_label"],
              "native_taxon":v["native_taxon"],
              "source_status":status,
              "q":v["q"],
              "daylength":v["daylength"]
            })
    if len(summaries)!=16 or sum(v["total_registered_native"] for v in summaries)!=sum(
            len(v["scientific_names"]) for v in selection):
        raise ValueError("Missing native source rows in first wave summary")

    complete={}
    for field in ("q","daylength"):
        source_key="source_q_equal_native_if_complete" if field=="q" else "source_daylength_equal_native_if_complete"
        selected={h:[s for s in summaries if s["stratum"]==h] for h in N}
        available={h:[s[source_key] for s in selected[h] if s[source_key] is not None]
                   for h in N}
        x=[y for val in available.values() for y in val]
        if len(x)>0:
            low,high=min(x),max(x)
            # Sensitivity conditional on all 16 already selected wave slots;
            # not a bound on 120 packages absent an all-frame outcome range.
            min_completed=sum(W[h]*
                (sum(available[h])+(8-len(available[h]))*low)/8 for h in N)
            max_completed=sum(W[h]*
                (sum(available[h])+(8-len(available[h]))*high)/8 for h in N)
        else:
            low=high=min_completed=max_completed=None
        all_strata_have_response=all(len(available[h])>0 for h in N)
        complete[field]={
          "frame_single_native":109,"frame_multi_native":11,
          "selected_single_native":8,"selected_multi_native":8,
          "first_wave_source_complete_response":len(x),
          "first_wave_source_unavailable_or_incomplete":16-len(x),
          "source_complete_response_counts":{h:len(available[h]) for h in N},
          "first_wave_stratum_equal_respondent_mean":
            {h:mean(available[h]) if available[h] else None for h in N},
          "response_only_stratum_poststratified_point_DESCRIPTIVE_NOT_DESIGN_UNBIASED":(
              sum(W[h]*mean(available[h]) for h in N)
              if all_strata_have_response else None),
          "observed_full_source_effect_range":([low,high] if x else None),
          "hypothetical_fill_missing_selected_sources_with_observed_range":
              ([min_completed,max_completed] if x else None),
          "nonresponding_and_unsampled_sources_do_NOT_become_zero":True,
          "no_qualified_world_population_inference":True,
        }
    return {
      "version":"RSDS_437_WAVE1_16_SOURCE_ORIGINAL_GPS_AGGREGATE_V1",
      "original_catalog_package_count":437,
      "original_GPS_GNSS_candidates":259,
      "metadata_eligible_finite_frame":120,
      "original_selection_wave_slots":16,
      "original_source_metadata_file_HOLD_slots":3,
      "source_status_counts":dict(sorted(statuses.items())),
      "all_first_wave_source_study_summaries":summaries,
      "all_first_wave_registered_native_group_effects":all_rows,
      "stratified_source_complete_effect_sensitivity":complete,
      "first_wave_method_provenance":"Pre-original-GPS SHA-ranked metadata selection; source identifiers and historical 48/ERA5 no-success preserved.",
      "not_a_world_population_meta_analysis":True,
      "no_original_animal_ID_dates_coordinates_or_raw_GPS":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    source=[json.loads((args.directory/f"wave1_source_{i:02d}.json").read_text())
            for i in range(16)]
    result=meta_snapshot(source)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                        encoding="utf-8")
    print("RSDS_WAVE1_ORIGINAL_GPS_SOURCE_META",json.dumps({
      "sources":result["source_status_counts"],
      "q":result["stratified_source_complete_effect_sensitivity"]["q"],
      "daylength":result["stratified_source_complete_effect_sensitivity"]["daylength"],
      "registered_taxa":len(result["all_first_wave_registered_native_group_effects"]),
      "not_world_meta":True
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
