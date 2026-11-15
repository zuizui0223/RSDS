"""Study-first, stage-resolved 16 + 19 original source meta-admission audit.

Adds the second source-ID SHA-ranked wave to the unchanged first 16-source
result. The estimand is source *measurability* under a fixed observation
protocol, not a global mean of animal activity, and nonresponse is never
zero-imputed. All registered native taxa within a source must pass the
original independent-tag rules to admit a source-complete q.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path
from statistics import mean,stdev

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned
from scripts.rsds_437_wave2_original_file_metadata_v1 import freeze as frozen_wave2_sources

# Aggregate-only validation reads pinned source manifests, NEVER GPS parsers.
# This module intentionally has no numpy/source-event dependency.
NEW_CONTRACT=("research/RSDS_437_WAVE2_ORIGINAL_GPS_FIRST_EFFECTS_V1_CONTRACT.json",
              "89e9b39aa107b76aa280f3dd04674713e7dbb9b1")

def freeze_original():
    c=pinned(*NEW_CONTRACT)
    if (c["version"]!="RSDS_437_WAVE2_19_ORIGINAL_GPS_FIRST_EFFECTS_CONTRACT_V1"
        or c["wave2_source_slots"]!=19
        or c["wave2_file_metadata_admissible"]!=15
        or c["wave2_file_size_HOLD"]!=3
        or c["wave2_human_privacy_HOLD"]!=1
        or c["historical_48_signed_irradiance_external"]!=
             "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
        or c["historical_12_ERA5_original_raw_SHA_HOLD"] is not True):
        raise ValueError("Original source-wave2 study meta contract changed")
    roster=frozen_wave2_sources()
    d=pinned(c["source_provider_file_metadata_path"],
             c["source_provider_file_metadata_git_blob"])
    if (len(roster)!=19 or len(d["source_rows"])!=19
        or d["source_status_counts"]!={
            "HUMAN_LOCATION_PRIVACY_HOLD":1,
            "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":15,
            "SOURCE_METADATA_HOLD":3}
        or d["all_new_wave2_animal_GPS_events_opened"]!=0):
        raise ValueError("Original source identity/metadata HOLD changed")
    for a,b in zip(roster,d["source_rows"]):
        if (a["wave2_index"]!=b["wave2_index"] or
            a["package_uuid"]!=b["package_uuid"] or
            a["original_repository_study_id"]!=b["original_repository_study_id"] or
            a["scientific_name_labels"]!=b["registered_native_taxa"]):
            raise ValueError("Wave2 original source taxonomy/ID pre-event changed")
    return c,d["source_rows"]

N={"single_native":109,"multi_native":11}
OLD_COUNT={"single_native":8,"multi_native":8}
NEW_COUNT={"single_native":16,"multi_native":3}
COMBINED={"single_native":24,"multi_native":11}
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"

def validate_effect(f):
    if f is None: return None
    ci=f.get("animal_bootstrap_95_percentile")
    mu=f.get("mean")
    n=f.get("n_individuals")
    if (not isinstance(n,int) or n<5 or
        not isinstance(mu,(int,float)) or not math.isfinite(mu) or
        not isinstance(ci,list) or len(ci)!=2 or
        not all(isinstance(x,(float,int)) and math.isfinite(x) for x in ci) or
        ci[0]>ci[1]):
        raise ValueError("Invalid native original group effect/uncertainty")
    return float(mu)

def covariance_status(record,original_taxa):
    z=record.get("native_joint_calendar_covariance") or {}
    status=z.get("status")
    if status=="SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED":
        taxa=z["original_native_taxa"]
        k=len(taxa)
        if (k<2 or taxa!=original_taxa or
            z["source_dates_tags_positions_and_draws_never_exported"] is not True or
            z["neither_global_meta_analysis_nor_activity_duration"] is not True):
            raise ValueError("Missing native identity/privacy in real covariance")
        for key in ("month_only","month_plus_positive_tag_weights"):
            part=z["joint_covariance_modes"][key]
            V=part["original_source_native_covariance_matrix"]
            if part["replicates"]!=1000 or len(V)!=2*k or any(len(row)!=2*k for row in V):
                raise ValueError("Source joint covariance dimensions changed")
            for i in range(2*k):
                if V[i][i]<0 or not math.isfinite(V[i][i]):
                    raise ValueError("Nonfinite source joint covariance")
                for j in range(2*k):
                    if abs(V[i][j]-V[j][i])>1e-10:
                        raise ValueError("Asymmetric source-native covariance")
    elif status not in (
        "NOT_EVALUATED","NOT_APPLICABLE_SINGLE_NATIVE_SOURCE",
        "ALL_NATIVE_COMMON_CALENDAR_COVARIANCE_HOLD","SOURCE_WIDE_HOLD"):
        raise ValueError("Unknown new original source joint covariance stage")
    return status

def summarize(first_results):
    if len(first_results)!=19:
        raise ValueError("All nineteen new source slots including HOLD required")
    contract,metadata=freeze_original()
    old=pinned(contract["wave1_original_first_result"],contract["wave1_first_git_blob"])
    w1=old["all_first_wave_source_study_summaries"]
    if len(w1)!=16 or old["original_source_metadata_file_HOLD_slots"]!=3:
        raise ValueError("Wave1 first original source evidence reclassified")
    source_ids={x["original_repository_study_id"] for x in w1}
    if len(source_ids)!=16:
        raise ValueError("Wave1 original independent study-ID duplicates")
    counts=Counter()
    by_wave=Counter()
    wave2=[]
    for index,(source,record) in enumerate(zip(metadata,first_results)):
        if (record.get("wave2_index")!=index or
            record.get("source_package_uuid")!=source["package_uuid"] or
            record.get("source_original_repository_study_id")!=source["original_repository_study_id"] or
            record.get("source_original_DOI_labels")!=source["original_DOI_labels"] or
            record.get("source_metadata_stratum")!=source["selection_stratum"] or
            record.get("original_registered_native_taxa")!=source["registered_native_taxa"] or
            len(record.get("native_group_results") or [])!=len(source["registered_native_taxa"]) or
            record.get("old_48_signed_solar_external_nonreplication_UNCHANGED")!=NEG or
            record.get("old_12_ERA5_source_raw_SHA_HOLD_UNCHANGED") is not True or
            record.get("original_raw_GPS_dates_tag_IDs_coordinates_exported") is not False or
            record.get("human_mobility_GPS_data_opened") is not False):
            raise ValueError("Wave2 source ID, native roster or previous biology changed")
        sid=source["original_repository_study_id"]
        if sid in source_ids:
            raise ValueError("Wave1/Wave2 source study-ID pseudoreplication")
        source_ids.add(sid)
        status=record["source_science_status"]
        allowed=(
            "SOURCE_GPS_HOLD","SOURCE_FILE_METADATA_HOLD",
            "HUMAN_LOCATION_PRIVACY_HOLD","SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE",
            "SOURCE_NO_NATIVE_Q_ESTIMABLE")
        if status not in allowed:
            raise ValueError("Unsupported source scientific status")
        if status=="HUMAN_LOCATION_PRIVACY_HOLD":
            if ("Homo sapiens" not in source["registered_native_taxa"] or
                record["original_file_MD5_bytes_verified"] or
                record["original_source_GPS_event_rows_opened"]!=0 or
                record["source_hold_stage"]!="HUMAN_SOURCE_PRIVACY_NO_ACCESS"):
                raise ValueError("Human GPS private source was opened")
        elif "Homo sapiens" in source["registered_native_taxa"]:
            raise ValueError("Homo sapiens source no-access gate removed")
        if (source["source_status"]=="SOURCE_METADATA_HOLD") != (
            status=="SOURCE_FILE_METADATA_HOLD"):
            raise ValueError("Earlier source file metadata HOLD reclassified")
        if (status in ("SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE","SOURCE_NO_NATIVE_Q_ESTIMABLE")
            and not record["original_file_MD5_bytes_verified"]):
            raise ValueError("Original q outcome admitted without MD5")
        native=record["native_group_results"]
        if [x["taxon"] for x in native]!=source["registered_native_taxa"]:
            raise ValueError("Original native species removed after outcomes")
        q=[];day=[]
        for x in native:
            qval=validate_effect(x["q"])
            dval=validate_effect(x["daylength"])
            if qval is not None: q.append(qval)
            if dval is not None: day.append(dval)
        if ((status=="SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE") != bool(q)):
            if status not in ("SOURCE_GPS_HOLD","SOURCE_FILE_METADATA_HOLD",
                              "HUMAN_LOCATION_PRIVACY_HOLD","SOURCE_NO_NATIVE_Q_ESTIMABLE"):
                raise ValueError("Declared Q success without supported native")
        if (status in ("SOURCE_GPS_HOLD","SOURCE_FILE_METADATA_HOLD",
                       "HUMAN_LOCATION_PRIVACY_HOLD","SOURCE_NO_NATIVE_Q_ESTIMABLE")
            and q):
            raise ValueError("Source HOLD or no-supported-native recoded positive")
        cov=covariance_status(record,source["registered_native_taxa"])
        if status in ("SOURCE_GPS_HOLD","SOURCE_FILE_METADATA_HOLD",
                      "HUMAN_LOCATION_PRIVACY_HOLD") and cov=="SOURCE_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE_ESTIMATED":
            raise ValueError("Joint covariance without admitted original GPS")
        k=len(source["registered_native_taxa"])
        complete_q=(len(q)==k)
        complete_day=(len(day)==k)
        counts[status]+=1
        by_wave["wave2_"+source["selection_stratum"]]+=1
        wave2.append({
            "wave":"wave2",
            "source_DOI_labels":source["original_DOI_labels"],
            "source_original_repository_study_id":sid,
            "source_stratum":source["selection_stratum"],
            "original_source_status":status,
            "registered_native_count":k,
            "native_q_supported":len(q),
            "native_daylength_supported":len(day),
            "complete_source_native_q":complete_q,
            "complete_source_daylength":complete_day,
            "source_q_mean_equal_registered_native":mean(q) if complete_q else None,
            "source_daylength_mean_equal_registered_native":mean(day) if complete_day else None,
            "joint_native_month_covariance_status":cov,
            "source_MD5_confirmed":record["original_file_MD5_bytes_verified"],
            "original_source_event_row_count":record["original_source_GPS_event_rows_opened"],
            "source_hold_stage":record["source_hold_stage"],
            "source_hold_error_type":record["source_hold_reason"],
        })

    prior_rows=[]
    for x in w1:
        q=x["source_q_equal_native_if_complete"]
        d=x["source_daylength_equal_native_if_complete"]
        prior_rows.append({
            "wave":"wave1",
            "source_DOI_labels":[x["original_study_DOI"]],
            "source_original_repository_study_id":x["original_repository_study_id"],
            "source_stratum":x["stratum"],
            "original_source_status":x["source_status"],
            "registered_native_count":x["total_registered_native"],
            "native_q_supported":x["native_q_supported"],
            "native_daylength_supported":x["native_daylength_supported"],
            "complete_source_native_q":x["registered_source_complete_q"],
            "complete_source_daylength":x["registered_source_complete_daylength"],
            "source_q_mean_equal_registered_native":q,
            "source_daylength_mean_equal_registered_native":d,
            "joint_native_month_covariance_status":x["native_covariance_status"],
            "source_MD5_confirmed":x["source_MD5_verified"],
            "original_source_event_row_count":x["original_source_animal_GPS_rows_opened"],
            "source_hold_stage":x["source_hold_stage"],
            "source_hold_error_type":x["source_hold_reason"],
        })
    combined=prior_rows+wave2
    if len(combined)!=35 or len(source_ids)!=35:
        raise ValueError("Original thirty-five source ID roster not intact")
    if Counter(z["source_stratum"] for z in combined)!=COMBINED:
        raise ValueError("Wave1+2 original source strata quotas changed")
    source_status=Counter((z["wave"],z["original_source_status"]) for z in combined)
    first_status=Counter(z["source_status"] for z in w1)
    if first_status!=Counter(old["source_status_counts"]):
        raise ValueError("Original Wave1 source status first result altered")

    estimands={}
    for field in ("q","daylength"):
        name="source_q_mean_equal_registered_native" if field=="q" else "source_daylength_mean_equal_registered_native"
        data={k:[r for r in combined if r["source_stratum"]==k] for k in N}
        avail={k:[float(r[name]) for r in data[k] if r[name] is not None] for k in N}
        all_values=[y for k in N for y in avail[k]]
        complete_count=sum(len(v) for v in avail.values())
        if all_values:
            lo=min(all_values);hi=max(all_values)
            # Explicitly impute only among 35 already selected sources.
            conditional=[
                sum((N[k]/120) *
                    ((sum(avail[k])+(COMBINED[k]-len(avail[k]))*fill)/COMBINED[k])
                    for k in N)
                for fill in (lo,hi)
            ]
        else:
            lo=hi=None;conditional=None
        respondents_exist=all(avail[k] for k in N)
        estimands[field]={
            "original_source_frame_single_native":109,
            "original_source_frame_multi_native":11,
            "sampled_original_study_source_count":35,
            "sampled_original_single_native":24,
            "sampled_original_multi_native":11,
            "supported_complete_native_source_count":complete_count,
            "complete_original_source_by_stratum":{k:len(v) for k,v in avail.items()},
            "missing_complete_native_source_count":35-complete_count,
            "respondent_only_stratum_means":{
                k:mean(v) if v else None for k,v in avail.items()
            },
            "respondent_only_poststratified_mean_EXPLORATORY":(
                sum((N[k]/120)*mean(avail[k]) for k in N)
                if respondents_exist else None
            ),
            "observed_complete_source_effect_range":([lo,hi] if all_values else None),
            "hypothetical_fill_missing_among_selected_35_with_observed_source_range":conditional,
            "no_design_unbiased_mean_because_native_source_nonresponse":True,
            "no_qualified_global_animal_population_estimate":True
        }
    measure={}
    for status_name in ("q","daylength"):
        key="complete_source_native_q" if status_name=="q" else "complete_source_daylength"
        by={k:sum(r[key] for r in combined if r["source_stratum"]==k) for k in N}
        total=sum(by.values())
        # Multi-native source metadata finite frame now fully enumerated (11/11);
        # single-native metadata sample 24/109 is still only a source-SHA prefix.
        measure[status_name]={
            "source_complete_within_sample_count":total,
            "by_source_stratum":by,
            "sampled_fraction_single":24/109,
            "sampled_fraction_multi":11/11,
            "observed_complete_rate_single":by["single_native"]/24,
            "census_complete_rate_multi":by["multi_native"]/11,
            "pseudo_random_SHA_frame_weighted_measurability":
                (109*by["single_native"]/24+by["multi_native"])/120,
            "strict_finite_frame_identified_measurability_bounds":
                [total/120,(total+85)/120],
            "not_a_world_prevalence_or_estimate_of_causal_measurement":True,
            "SHA_hash_sampling_and_post_Wave1_quota_caveat":True,
        }

    return {
      "version":"RSDS_WAVE1_PLUS_WAVE2_35_ORIGINAL_SOURCE_STUDY_META_ADMISSION_V1",
      "research_stage":"POST_OLD_WAVE1_OUTCOME_WAVE2_SOURCE_SHA_SELECTED_BEFORE_NEW_EVENT_ACCESS",
      "original_registry_source_packages":437,
      "original_GPS_GNSS_sensor_packages":259,
      "original_metadata_candidate_frame":120,
      "original_source_study_units":35,
      "wave1_source_study_units":16,
      "wave2_source_study_units":19,
      "wave2_status_counts":dict(sorted(counts.items())),
      "combined_wave_source_status_counts":{
         wave:{key.split("::",1)[1]:num
               for key,num in sorted({
                   f"{label}::{status}":count for (label,status),count in source_status.items()
                   if label==wave}.items())}
         for wave in ("wave1","wave2")
      },
      "combined_original_registered_native_rows":
           sum(z["registered_native_count"] for z in combined),
      "combined_source_study_summary_rows":combined,
      "new_wave2_source_study_summary_rows":wave2,
      "source_complete_effects_and_nonresponse_sensitivity":estimands,
      "finite_120_metadata_frame_source_measurability":measure,
      "previous_48_source_sunlight_signed_external_not_replicated":True,
      "previous_12_ERA5_raw_SHA_HOLD":True,
      "human_location_source_no_access":True,
      "no_original_GPS_tracks_IDs_exact_dates_coordinates_or_bootstrap_draws_exported":True,
      "scientific_interpretation":"Exploratory, partially nonresponding source-meta analysis, distinct from population-level animal activity or causality; wave2 quota chosen after wave1 outcomes although original within-stratum SHA rank precedes them."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    data=[json.loads((a.directory/f"wave2_source_{i:02d}.json").read_text(encoding="utf-8"))
          for i in range(19)]
    result=summarize(data)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_WAVE2_35_SOURCE_META",json.dumps({
        "source_status":result["wave2_status_counts"],
        "registered_native_count":result["combined_original_registered_native_rows"],
        "source_measurability":result["finite_120_metadata_frame_source_measurability"],
        "respondent_only_effects":result["source_complete_effects_and_nonresponse_sensitivity"],
        "no_global_population_inference":True
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
