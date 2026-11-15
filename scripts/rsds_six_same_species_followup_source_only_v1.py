"""Source-only, post-discovery six taxon-target follow-up intake (NO new GPS).

Targets chosen after the 63-source first result; therefore NOT confirmation of
world-scale species regularities. All eligible source slots are kept regardless
of unknown future GPS q sign, provider file validity, or tagging support.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned
from scripts.rsds_original_source_native_same_species_repeatability_v1 import norm_doi

CONTRACT=HERE/"research/RSDS_SIX_SAME_SPECIES_FOLLOWUP_SOURCE_ONLY_V1_CONTRACT.json"
VERSION="RSDS_SAME_SPECIES_SIX_FOLLOWUP_ORIGINAL_SOURCE_METADATA_SELECTION_V1"

def freeze():
    c=json.loads(CONTRACT.read_text())
    if (c.get("version")!=VERSION or
        c.get("original_public_catalog_packages")!=437 or
        c.get("original_63_source_studies")!=63 or
        c.get("expected_new_six")!=6 or
        c.get("expected_target_names_with_candidate")!=5 or
        c.get("no_original_provider_file_list_queried") is not True or
        c.get("no_original_GPS_event_or_tag_data_read") is not True or
        c.get("no_outcome_based_replacement") is not True or
        c.get("old_solar_48_source_external_signed_direction")!="NOT_REPLICATED" or
        c.get("old_12_ERA5_raw_SHA_HOLD") is not True):
        raise ValueError("Original six-target source-only contract invalid")
    cat=pinned(c["catalog_path"],c["catalog_git_blob"])
    atlas=pinned(c["original_63_source_native_atlas_path"],
                 c["original_63_source_native_atlas_git_blob"])
    prior=pinned(c["original_prior_12_path"],c["original_prior_12_git_blob"])
    wave1=pinned(c["original_wave1_path"],c["original_wave1_git_blob"])
    wave2=pinned(c["original_wave2_path"],c["original_wave2_git_blob"])
    wave3=pinned(c["original_wave3_path"],c["original_wave3_git_blob"])
    if len(cat["packages"])!=437 or atlas["distinct_original_archived_source_DOIs_and_study_IDs"]!=63:
        raise ValueError("Source-only archive count changed")
    names={x["scientific_name"] for x in
           atlas["same_species_source_study_q_comparison_candidates"]}
    if names!=set(c["target_scientific_names_from_prior_post_outcome_repeatability"]):
        raise ValueError("Previously studied taxon target identity changed")
    old={norm_doi(x["original_study_doi"]) for x in
         prior["all_original_study_native_taxon_rows_including_HOLD"]}
    old|={norm_doi(x["source_DOI_label"]) for x in
          wave1["pre_event_selected_sources_in_rank_order"]}
    old|={norm_doi(x["original_DOI_labels"][0]) for x in
          wave2["all_original_wave2_sources_in_rank_order"]}
    old|={norm_doi(x["original_DOI_labels"][0]) for x in
          wave3["all_wave3_selected_source_metadata_in_SHA_rank_order"]}
    if len(old)!=63:
        raise ValueError("Previously evaluated old 12+51 original source DOIs changed")
    return c,cat["packages"],names,old

def evaluate():
    c,packages,names,old=freeze()
    rejection=[]
    selected=[]
    observed_ids=defaultdict(list)
    for p in packages:
        for sid in p["repository_study_ids"]:
            observed_ids[sid].append(p["package_uuid"])
    for p in packages:
        sp=sorted(set(p.get("scientific_names",[]))&names)
        if not sp or any(norm_doi(x) in old for x in p.get("dois") or []):
            continue
        issue=[]
        if not (set(p.get("sensor_names") or [])&{"GPS","GNSS"}):
            issue.append("NOT_GPS_GNSS")
        if p.get("licenses")!=["CC0 1.0 Universal"]:
            issue.append("NON_CC0_SOURCE_LICENSE")
        if (not isinstance(p.get("reported_animal_count"),int) or
            p["reported_animal_count"]<10 or
            not isinstance(p.get("reported_location_count"),int) or
            p["reported_location_count"]<10000):
            issue.append("SOURCE_CATALOG_TOO_FEW_REPORTED_ANIMALS_OR_LOCATIONS")
        if len(p.get("dois") or [])!=1 or len(p.get("repository_study_ids") or [])!=1:
            issue.append("AMBIGUOUS_SOURCE_ID_OR_DOI")
        elif len(observed_ids[p["repository_study_ids"][0]])!=1:
            issue.append("DUPLICATED_ORIGINAL_STUDY_ID")
        if not p.get("scientific_names") or len(set(p["scientific_names"]))!=len(p["scientific_names"]):
            issue.append("INVALID_NATIVE_TAXON_LABELS")
        row={
            "original_DOI":norm_doi(p["dois"][0]),
            "original_repository_study_id":p["repository_study_ids"][0],
            "original_public_package_uuid":p["package_uuid"],
            "original_scientific_names":p["scientific_names"],
            "target_matched_scientific_names":sp,
            "original_title":p["title"],
            "original_license_labels":p["licenses"],
            "original_reported_animals_metadata_only":p["reported_animal_count"],
            "original_reported_GPS_locations_metadata_only":p["reported_location_count"],
            "original_sensor_labels":p["sensor_names"],
            "source_prior_63_q_study":False,
            "provider_original_file_metadata_opened":False,
            "original_animal_GPS_events_opened":0,
        }
        if issue:
            row["source_only_status"]="SOURCE_CATALOG_NOT_ELIGIBLE_NO_SUBSTITUTION"
            row["source_only_reasons"]=issue
            rejection.append(row)
        else:
            row["source_only_status"]="ORIGINAL_SOURCE_METADATA_CANDIDATE_NOT_GPS_ADMITTED"
            row["source_only_reasons"]=[]
            selected.append(row)
    selected.sort(key=lambda x:(x["original_DOI"],x["original_public_package_uuid"]))
    rejection.sort(key=lambda x:(x["original_DOI"],x["original_public_package_uuid"]))
    actual={(x["original_DOI"],x["original_repository_study_id"],x["target_matched_scientific_names"][0])
            for x in selected}
    expected={(x["doi"],x["source_study_id"],x["native"])
              for x in c["fixed_new_six_original_candidates"]}
    if actual!=expected or len(selected)!=6:
        raise ValueError("Six original source-only target identities changed")
    if any(x["original_DOI"] in old for x in selected):
        raise ValueError("Previously seen original source was improperly selected")
    counter=Counter(sp for x in selected for sp in x["target_matched_scientific_names"])
    if len(counter)!=5 or counter["Anser albifrons"]!=2:
        raise ValueError("Original target source study composition changed")
    return {
        "version":VERSION,
        "stage":"AFTER_SOURCE_Q_DISCOVERY_BEFORE_SIX_NEW_PROVIDER_OR_GPS_READS",
        "all_original_public_source_packages":437,
        "original_previous_q_source_studies_EXCLUDED":63,
        "target_scientific_names_post_discovery":sorted(names),
        "original_source_metadata_eligible_six":len(selected),
        "eligible_by_target_taxon":dict(sorted(counter.items())),
        "original_target_metadata_candidates_eligible":selected,
        "all_other_original_target_metadata_candidates_noneligible":rejection,
        "other_noneligible_target_metadata_count":len(rejection),
        "target_names_with_no_new_source":[x for x in sorted(names) if x not in counter],
        "no_original_provider_file_metadata_or_GPS_opened":True,
        "individual_population_independence_verified":False,
        "other_RSDS_ODSP_prior_original_GPS_source_overlap_verified":False,
        "source_eligibility_targeted_after_original_q_outcomes":True,
        "no_global_species_confirmation_or_meta_probability_sample":True,
        "old_solar_48_source_external_signed_direction_replicated":False,
        "old_12_ERA5_raw_SHA_HOLD":True,
        "inference_boundary":c["inference_boundary"],
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    x=evaluate()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("SAME_SPECIES_SIX_FOLLOWUP_SOURCE_ONLY",json.dumps({
        "selected":x["original_source_metadata_eligible_six"],
        "source_dois":[y["original_DOI"] for y in x["original_target_metadata_candidates_eligible"]],
        "species":x["eligible_by_target_taxon"],
        "source_only_exclusions":x["other_noneligible_target_metadata_count"],
        "original_GPS_events_opened":0,
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
