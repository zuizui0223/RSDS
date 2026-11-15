"""Frozen six-source metadata-only prior solar/GPS-exposure audit.

This scans eight exact original JSON result surfaces only; not a full search of
all archives and not a proof of different animal tags or population cohorts.
No original animal telemetry, provider file inventory or GPS is requested.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned
from scripts.rsds_original_source_native_same_species_repeatability_v1 import norm_doi

SPEC=HERE/"research/RSDS_SIX_FOLLOWUP_ARCHIVED_SOURCE_EXPOSURE_V1_CONTRACT.json"
VERSION="RSDS_SIX_FOLLOWUP_ARCHIVED_SOURCE_EXPOSURE_AUDIT_V1"

def matching_original_source_nodes(node,doi,path="root"):
    hits=[]
    if isinstance(node,dict):
        ids=node.get("dois")
        if isinstance(ids,list) and any(norm_doi(x)==doi for x in ids):
            hits.append({
                "source_path":path,
                "original_source_state":node.get("source_state"),
                "provider_original_checksum_verified":(
                    bool((node.get("original_provider_checksum_verified") or {}).get(
                        "provider_checksum_verified",False))
                    if isinstance(node.get("original_provider_checksum_verified"),dict)
                    else False
                ),
                "original_source_coverage_rows_reported":(
                    (node.get("source_original_coverage") or {}).get("original_total_rows")
                    if isinstance(node.get("source_original_coverage"),dict)
                    else None
                ),
                "archive_hold_reason":node.get("reason"),
            })
        for key,value in node.items():
            hits+=matching_original_source_nodes(value,doi,path+"."+key)
    elif isinstance(node,list):
        for i,value in enumerate(node):
            hits+=matching_original_source_nodes(value,doi,path+"["+str(i)+"]")
    return hits

def evaluate():
    c=json.loads(SPEC.read_text())
    if (c.get("version")!=VERSION
        or len(c.get("eight_preexisting_solar_and_temporal_original_output_receipts") or [])!=8
        or c.get("frozen_six_source_manifest_git_blob")!="1315b435c40b59831a643e35e9b7f9063879673b"
        or c.get("no_provider_file_listing_or_GPS_reopen") is not True
        or c.get("prior_48_solar_direction_signed_external_nonreplicated") is not True
        or c.get("prior_12_ERA5_raw_SHA_HOLD") is not True):
        raise ValueError("Previous archive exposure source-contract changed")
    selected=pinned(c["frozen_six_source_manifest_path"],
                    c["frozen_six_source_manifest_git_blob"])
    sources=selected["original_target_metadata_candidates_eligible"]
    if len(sources)!=6 or not selected["no_original_provider_file_metadata_or_GPS_opened"]:
        raise ValueError("Original six source-only identity freeze changed")
    loaded={}
    for f in c["eight_preexisting_solar_and_temporal_original_output_receipts"]:
        loaded[f["path"]]=pinned(f["path"],f["git_blob"])
    result=[]
    for s in sources:
        doi=s["original_DOI"]
        matches=[]
        for path,raw in loaded.items():
            for h in matching_original_source_nodes(raw,doi):
                matches.append({"prior_original_result_file":path,**h})
        status="NO_MATCH_IN_EIGHT_PINNED_PRIOR_RESULT_SURFACES_NOT_EXHAUSTIVE"
        if any(x["original_source_state"]=="original_MD5_verified_source_analyzed"
               and x["provider_original_checksum_verified"] for x in matches):
            status="KNOWN_PREVIOUS_ORIGINAL_GPS_DATA_ANALYZED_FOR_SOLAR_DIRECTION"
        elif any(x["original_source_state"]=="SOURCE_HOLD" for x in matches):
            status="KNOWN_PREVIOUS_SOURCE_METADATA_HOLD_NO_GPS_FROM_THAT_ROUTE"
        result.append({
            "original_DOI":doi,
            "original_repository_study_id":s["original_repository_study_id"],
            "original_public_package_uuid":s["original_public_package_uuid"],
            "original_scientific_names":s["original_scientific_names"],
            "prior_eight_output_matches":matches,
            "prior_data_exposure_status":status,
            "prior_other_uninspected_outputs_and_tag_overlap_unresolved":True,
            "new_original_GPS_rows_reopened":0,
        })
    old_fox=next(x for x in result if x["original_DOI"]=="10.5441/001/1.72hh609t")
    old_duck=next(x for x in result if x["original_DOI"]=="10.5441/001/1.8dc0v84m")
    if (old_fox["prior_data_exposure_status"]!="KNOWN_PREVIOUS_ORIGINAL_GPS_DATA_ANALYZED_FOR_SOLAR_DIRECTION"
        or old_duck["prior_data_exposure_status"]!="KNOWN_PREVIOUS_SOURCE_METADATA_HOLD_NO_GPS_FROM_THAT_ROUTE"
        or sum(bool(x["prior_eight_output_matches"]) for x in result)!=2):
        raise ValueError("Prior original GPS observation exposure classification changed")
    return {
        "version":VERSION,
        "stage":"SOURCE_ONLY_EIGHT_ARCHIVE_PREVIOUS_EXPOSURE_DIAGNOSTIC",
        "previous_fixed_six_source_only_selection_git_blob":c["frozen_six_source_manifest_git_blob"],
        "number_of_original_prior_archive_surfaces_inspected":8,
        "original_six_study_IDs_retained":6,
        "known_previous_GPS_analysis_sources":1,
        "known_previous_original_provider_metadata_HOLD_sources":1,
        "no_match_in_eight_checked_previous_outputs_sources":4,
        "original_followup_source_exposure_rows":result,
        "global_prior_GPS_exposure_audit_completed":False,
        "same_animal_independence_verified":False,
        "no_original_source_GPS_animal_events_reopened":True,
        "prior_original_48_source_solar_direction_replicated":False,
        "prior_12_ERA5_raw_SHA_HOLD":True,
        "limited_scope_warning":c["audit_scope"],
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    x=p.parse_args()
    result=evaluate()
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("RSDS_SIX_SOURCE_PREVIOUS_ARCHIVE_EXPOSURE",json.dumps({
        "prior_GPS_analyzed":result["known_previous_GPS_analysis_sources"],
        "prior_metadata_HOLD":result["known_previous_original_provider_metadata_HOLD_sources"],
        "not_found_among_eight":result["no_match_in_eight_checked_previous_outputs_sources"],
        "source_DOI_status":[{"doi":x["original_DOI"],
                             "status":x["prior_data_exposure_status"]}
                            for x in result["original_followup_source_exposure_rows"]],
    },sort_keys=True))
if __name__=="__main__":
    main()
