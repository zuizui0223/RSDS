"""Study-first aggregate of five paired matched-original-tag/day slope decompositions."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import json
import math
from pathlib import Path
from scripts.rsds_five_source_same_tag_date_month_decomposition_v1 import fixed,FIRST_RESULT_VERSION

VERSION="RSDS_FIVE_SOURCE_SAME_TAG_MONTH_MATCHED_META_FIRST_V1"


def aggregate(records):
    contract,orig,prior=fixed()
    if len(records)!=5:
        raise ValueError("All five previous original study source IDs mandatory")
    output=[]
    for idx,(r,source,priorrow) in enumerate(zip(records,orig["source_manifest"],prior)):
        if (r.get("version")!=FIRST_RESULT_VERSION or
            r.get("source_index")!=idx or
            r.get("original_source_key")!=source["source_key"] or
            r.get("original_DOI")!=source["doi"] or
            r.get("source_original_study_id")!=str(source["study_id"]) or
            r.get("source_native_taxon")!=source["taxon"] or
            r.get("all_source_registered_taxa")!=source["registered_taxa"] or
            r.get("historical_48_solar_heading_negative_external_verdict_unchanged") is not True or
            r.get("historical_12_ERA5_raw_SHA_HOLD_unchanged") is not True or
            r.get("source_original_GPS_IDs_coordinates_exact_dates_and_individual_effects_never_exported") is not True):
            raise ValueError("Source/DOI/native/history altered in matched-month comparison")
        status=r["source_status"]
        if status not in ("MATCHED_SUPPORT_PAIRED_SOURCE_ESTIMABLE",
                          "MATCHED_SUPPORT_NO_NATIVE_TAGS","SOURCE_HOLD"):
            raise ValueError("Unrecognized matched source outcome")
        p=r.get("matched_month_effect")
        if status=="SOURCE_HOLD" and p is not None:
            raise ValueError("Partial source HOLD original-data estimate leaked")
        if status!="SOURCE_HOLD" and not (
            r["old_q_and_daylength_original_baseline_replayed"] and
            r["first_month_FE_identical_sample_effect_replayed"] and
            r["original_file_UUID_MD5_bytes_verified"]):
            raise ValueError("Source month effects accepted without replayed original MD5/baseline")
        oldmonth=priorrow["month_native"]
        if p and p["source_original_eligible_tagged_animals"]!=oldmonth["original_tagged_animals_with_year_month_conditions"]:
            raise ValueError("Changed eligible independent original tag sample in matched source")
        effect=p["source_q_daylength_effects"] if p else {}
        for mode in ("local","reference"):
            keys=[mode+"_raw_same_dates",mode+"_month_FE_same_dates",mode+"_month_minus_raw"]
            if effect:
                for key in keys:
                    q=effect[key]
                    if q["status"]=="ESTIMABLE":
                        ci=q["tag_bootstrap95"]
                        if (not isinstance(q["mean"],(float,int)) or not math.isfinite(q["mean"]) or
                            len(ci)!=2 or not all(math.isfinite(x) for x in ci) or ci[0]>ci[1]):
                            raise ValueError("Invalid paired tag-bootstrap source effect")
                        if q["nominal_sign"]!=(
                            "POSITIVE" if ci[0]>0 else
                            "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"):
                            raise ValueError("Paired tag bootstrap source sign overwritten")
                    elif q["mean"] is not None:
                        raise ValueError("Insufficient paired native tag result must be null")
                if effect[keys[0]]["tag_n"]!=effect[keys[1]]["tag_n"]:
                    raise ValueError("Matched raw-month tagged original sample differs")
                if effect[keys[0]]["status"]=="ESTIMABLE":
                    diff=effect[keys[1]]["mean"]-effect[keys[0]]["mean"]
                    if not math.isclose(diff,effect[keys[2]]["mean"],rel_tol=0,abs_tol=1e-12):
                        raise ValueError("Algebraic month-minus-matched-raw tag mean differs")
        raw=effect.get("local_raw_same_dates",{})
        adj=effect.get("local_month_FE_same_dates",{})
        delta=effect.get("local_month_minus_raw",{})
        output.append({
            "source_index":idx,
            "source_key":source["source_key"],
            "original_DOI":source["doi"],
            "source_native_species":source["taxon"],
            "first_full_source_native_slope":source["baseline"]["daylength"],
            "first_month_FE_source_slope":oldmonth["local_original_solar_daylength"],
            "same_tag_date_unadjusted":raw,
            "same_tag_date_month_adjusted":adj,
            "same_tag_date_month_minus_raw":delta,
            "same_tag_date_fixed_reference_latitude_unadjusted":effect.get("reference_raw_same_dates"),
            "same_tag_date_fixed_reference_latitude_month_adjusted":effect.get("reference_month_FE_same_dates"),
            "same_tag_date_fixed_reference_latitude_month_minus_raw":effect.get("reference_month_minus_raw"),
            "source_status":status,
            "all_contrasts_same_original_tag_and_date_support":p["all_contrasts_use_identical_original_tag_year_month_date_support"] if p else False,
            "source_qualified_distinct_original_tags":p["source_original_eligible_tagged_animals"] if p else None,
            "source_hold_stage":r["source_hold_stage"],
            "source_hold_type":r["source_hold_error_type"],
        })
    by=defaultdict(list)
    for x in output:by[x["source_native_species"]].append(x)
    species={}
    for sp,source_rows in sorted(by.items()):
        species[sp]={
            "original_source_studies":len(source_rows),
            "paired_month_support_estimable_studies":sum(x["source_status"]=="MATCHED_SUPPORT_PAIRED_SOURCE_ESTIMABLE" for x in source_rows),
            "paired_adjustment_nominal_signs_by_source":{
                x["source_key"]:x["same_tag_date_month_minus_raw"].get("nominal_sign")
                for x in source_rows if x["same_tag_date_month_minus_raw"]},
            "within_month_source_nominal_signs_by_source":{
                x["source_key"]:x["same_tag_date_month_adjusted"].get("nominal_sign")
                for x in source_rows if x["same_tag_date_month_adjusted"]},
            "tag_population_overlap_and_sampling_comparability_unresolved":True,
        }
    return {
        "version":VERSION,
        "stage":"POSTHOC_SOURCE_NATIVE_MATCHED_ORIGINAL_TAG_DATE_SUPPORT_MONTH_DECOMPOSITION",
        "original_source_study_DOI_total":5,
        "original_species_target_count":2,
        "original_first_month_result_git_blob_unchanged":contract["prior_five_source_month_first_result_git_blob"],
        "source_status_counts":dict(sorted(Counter(z["source_status"] for z in output).items())),
        "all_five_original_study_source_native_matched_month_results_including_HOLD":output,
        "original_same_species_matched_support_summary":species,
        "all_original_source_tags_and_dates_held_identical_within_paired_comparison":True,
        "original_source_specific_tag_bootstrap_not_world_meta_CI":True,
        "post_outcome_method_diagnostic_not_causal_daylength_effect":True,
        "independent_populations_and_recruitment_not_verified":True,
        "previous_48_source_signed_solar_heading_negative_replication_UNCHANGED":True,
        "previous_12_source_ERA5_raw_SHA_HOLD_UNCHANGED":True,
        "no_original_tag_IDs_coordinates_dates_individual_slopes_or_bootstrap_draws_exported":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    arg=p.parse_args()
    files=list(arg.directory.glob("source_*.json"))
    if len(files)!=5:raise ValueError("Five original source artifacts including HOLD required")
    by={}
    for f in files:
        obj=json.loads(f.read_text(encoding="utf-8"))
        i=obj.get("source_index")
        if not isinstance(i,int) or i not in range(5) or i in by:
            raise ValueError("Duplicate or missing original study source index")
        by[i]=obj
    if len(by)!=5:raise ValueError("Incomplete 5 original source results")
    z=aggregate([by[i] for i in range(5)])
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(z,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_FIVE_MATCHED_ORIGINAL_TAG_DATE",json.dumps({
        "source_status":z["source_status_counts"],
        "species":z["original_same_species_matched_support_summary"],
        "source_month_estimator":"postoutcome_statistical_decomposition_not_causal"
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
