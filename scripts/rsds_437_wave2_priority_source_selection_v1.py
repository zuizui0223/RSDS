"""Second, 19-source response-blind SHA priority wave from first Movebank frame.

The QUOTA was decided after Wave1 outcomes (explicit adaptive study-stage
selection); but the 16 single-source and 3 remaining multi-source identities
are chosen only by the ORIGINAL fixed, metadata-only Wave1 SHA ranking,
before any Wave2 source provider file metadata or GPS events are opened.
No new animal responses, species choice, place title or daylength effect enter.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
import hashlib
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import HERE, run as first_frame, pinned, git_blob

CONTRACT=HERE/"research/RSDS_437_SOURCE_FRAME_WAVE2_PRE_EVENT_CONTRACT_V1.json"
EXPECTED_VERSION="RSDS_437_SOURCE_FRAME_WAVE2_PRE_EVENT_SELECTION_V1"
BASE_SHA="1fbb6e4f7ef686aed63dfd4433ee425e751592ca"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"

def freeze():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (c.get("version")!=EXPECTED_VERSION or
        c.get("original_frame_git_blob")!=BASE_SHA or
        c.get("candidate_stage1_total")!=120 or
        c.get("single_frame_N")!=109 or c.get("multi_frame_N")!=11 or
        c.get("wave1_single")!=8 or c.get("wave1_multi")!=8 or
        c.get("wave2_single")!=16 or c.get("wave2_multi")!=3 or
        c.get("wave2_total")!=19 or c.get("combined_total")!=35 or
        c.get("stage1_animal_events_opened")!=0 or
        c.get("historical_48_external_sunlight")!=NEG or
        c.get("historical_ERA5_raw_SHA_HOLD") is not True):
        raise ValueError("Frozen second-wave source quotas/history changed")
    original_source=HERE/c["original_frame_script"]
    if git_blob(original_source.read_bytes())!=BASE_SHA:
        raise ValueError("Original Wave1 437-source frame sorting code changed")
    prior=pinned(c["prior_wave1_selection"],c["prior_wave1_selection_git_blob"])
    files=pinned(c["prior_wave1_file_metadata"],c["prior_wave1_file_metadata_git_blob"])
    first=pinned(c["prior_wave1_GPS_first_result"],c["prior_wave1_GPS_first_git_blob"])
    if (prior["original_registry_package_count"]!=437 or
        prior["meta_stage1_candidate_count"]!=120 or
        prior["meta_stage1_single_candidate_count"]!=109 or
        prior["meta_stage1_multi_candidate_count"]!=11 or
        len(prior["pre_event_selected_sources_in_rank_order"])!=16 or
        first["original_selection_wave_slots"]!=16 or
        first["stratified_source_complete_effect_sensitivity"]["q"]["first_wave_source_complete_response"]!=4 or
        files["original_package_slots"]!=16 or
        files["original_GPS_event_records_opened"]!=0):
        raise ValueError("Original Wave1 source selection, file HOLD or q first result changed")
    return c,prior

def select():
    c,prior=freeze()
    frame=first_frame()
    if (frame["summary"]["metadata_candidate_packages"]!=120 or
        frame["summary"]["metadata_candidate_single_native"]!=109 or
        frame["summary"]["metadata_candidate_multi_native"]!=11 or
        frame["summary"]["wave1_selected_packages"]!=16):
        raise ValueError("Frozen 120 candidates no longer present")
    rank={}
    for stratum,count in (("single_native",109),("multi_native",11)):
        ranks=sorted(
           (r for r in frame["all_437_metadata_status_rows"]
            if r["stage1_status"]=="FIRST_WAVE_METADATA_CANDIDATE"
            and r["selection_stratum"]==stratum),
           key=lambda r:(r["lottery_sha256"],r["package_uuid"]))
        if len(ranks)!=count:
            raise ValueError("Old wave selection stratum denominator changed")
        rank[stratum]=ranks
    prior_wave=prior["pre_event_selected_sources_in_rank_order"]
    expected=rank["single_native"][:8]+rank["multi_native"][:8]
    if [(x["package_uuid"],x["lottery_sha256"]) for x in expected]!=[
        (x["package_uuid"],frame_row["lottery_sha256"])
        for x,frame_row in zip(prior_wave,expected)]:
        raise ValueError("Wave1 original SHA rank prefix not replayed")
    if [x["package_uuid"] for x in expected]!=[z["package_uuid"] for z in prior_wave]:
        raise ValueError("Wave1 prior package identities no longer rank-prefix")
    base_selected={z["package_uuid"] for z in prior_wave}
    next_single=rank["single_native"][8:24]
    next_multi=rank["multi_native"][8:11]
    chosen=next_single+next_multi
    if (len(next_single)!=16 or len(next_multi)!=3 or len(chosen)!=19
        or any(r["package_uuid"] in base_selected for r in chosen)
        or len({r["repository_study_ids"][0] for r in chosen})!=19
        or len({r["doi_labels"][0] for r in chosen})!=19):
        raise ValueError("Wave2 ranks repeated previously evaluated source study")
    rows=[]
    for i,x in enumerate(chosen):
        sp=x["selection_stratum"]
        k=8 + i+1 if sp=="single_native" else 8+(i-16)+1
        rows.append({
            "wave2_index":i,"second_wave_ordinal":i+1,
            "stratum_rank_1_based":k,
            "package_uuid":x["package_uuid"],
            "original_DOI_labels":x["doi_labels"],
            "original_repository_study_id":x["repository_study_ids"][0],
            "scientific_name_labels":x["scientific_name_labels"],
            "original_sensor_labels":x["original_sensors"],
            "original_license_labels":x["original_license_labels"],
            "metadata_reported_animals":x["reported_animal_count_not_deduplicated"],
            "metadata_reported_locations":x["reported_location_count_not_qualified"],
            "stratum":sp,
            "source_lottery_SHA256":x["lottery_sha256"],
            "original_source_file_metadata_opened":False,
            "original_GPS_animal_events_opened":0,
        })
    return {
      "version":"RSDS_437_SHA_PRIORITY_WAVE2_19_SOURCE_SELECTION_V1",
      "research_stage":"POST_WAVE1_OUTCOME_QUOTA_DECISION_PRE_WAVE2_FILE_METADATA_AND_GPS_EVENTS",
      "original_source_catalog_SHA":c["original_catalog_git_blob"],
      "original_wave1_source_selection_git_blob":c["prior_wave1_selection_git_blob"],
      "original_wave1_biological_first_result_git_blob":c["prior_wave1_GPS_first_git_blob"],
      "sampling_method":"CONTINUE_IDENTICAL_METADATA_SHA_RANKS_SINGLE_9_TO_24_MULTI_9_TO_11",
      "original_candidate_source_frame":{"single_native":109,"multi_native":11},
      "wave1_source_counts":{"single_native":8,"multi_native":8},
      "wave2_new_source_counts":{"single_native":16,"multi_native":3},
      "two_wave_combined_counts":{"single_native":24,"multi_native":11},
      "wave2_new_original_source_count":19,
      "all_original_wave2_sources_in_rank_order":rows,
      "all_previously_unselected_metadata_candidates_remaining":85,
      "original_wave2_GPS_events_opened":0,
      "original_wave2_provider_files_opened":0,
      "old_48_signed_sunlight_external_replicated":False,
      "old_ERA5_original_12_raw_SHA_HOLD":True,
      "selection_boundaries":[
         "Source SHA rank was locked before Wave1 GPS outcomes; Wave2 quotas were selected afterward.",
         "SHA rank is metadata-pseudorandom conditional on the original frozen catalog, not a global wildlife probability sample.",
         "Both source-file and GPS nonresponse remain source-specific and cannot be dropped or recoded to zero.",
         "All three remaining multi-native source-study IDs are retained even if original records are unusable."
      ]
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",required=True,type=Path)
    a=parser.parse_args()
    d=select()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_WAVE2_SOURCE_SELECTION",json.dumps({
        "sources":19,"single":16,"multi":3,"combined_two_wave":35,
        "remaining":85,
        "new_sources":[
            {"index":x["wave2_index"],"DOI":x["original_DOI_labels"],
             "study_id":x["original_repository_study_id"],"stratum":x["stratum"]}
            for x in d["all_original_wave2_sources_in_rank_order"]
        ],"new_GPS_events":0
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
