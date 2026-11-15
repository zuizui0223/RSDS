"""Provider ORIGINAL file metadata ONLY for six frozen targeted follow-up source IDs.

The earlier mallard source was already a source-only metadata HOLD in a frozen
solar project and remains HOLD, with no provider requery. Another fox source
has previously been analyzed for solar direction; it cannot be called a
wholly untouched external GPS study even if q was not yet estimated.
Every original source remains in the denominator; no GPS bytes/events read.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned
from scripts.rsds_six_same_species_followup_source_only_v1 import evaluate as replay_selection
from scripts.rsds_six_prior_original_source_exposure_audit_v1 import evaluate as replay_exposure
from scripts.rsds_437_wave1_original_file_metadata_v1 import original_primary
from scripts.movebank_native_ecology_20261019 import original_files

SIX_SELECTION_PATH="research/RSDS_SIX_SAME_SPECIES_FOLLOWUP_ORIGINAL_SOURCE_ONLY_FIRST_RESULT_20261115.json"
SIX_SELECTION_BLOB="1315b435c40b59831a643e35e9b7f9063879673b"
EXPOSURE_PATH="research/RSDS_SIX_FOLLOWUP_PRIOR_ORIGINAL_SOURCE_EXPOSURE_FIRST_RESULT_20261115.json"
EXPOSURE_BLOB="51c078c0fc9f1353de4042ec3befdb3aa7b97f35"
VERSION="RSDS_SIX_FOLLOWUP_ORIGINAL_PROVIDER_FILE_METADATA_FIRST_V1"

def freeze():
    sel=pinned(SIX_SELECTION_PATH,SIX_SELECTION_BLOB)
    exp=pinned(EXPOSURE_PATH,EXPOSURE_BLOB)
    if (sel!=replay_selection() or exp!=replay_exposure()
        or sel["original_source_metadata_eligible_six"]!=6
        or not sel["no_original_provider_file_metadata_or_GPS_opened"]
        or exp["original_six_study_IDs_retained"]!=6
        or exp["known_previous_GPS_analysis_sources"]!=1
        or exp["known_previous_original_provider_metadata_HOLD_sources"]!=1
        or not exp["no_original_source_GPS_animal_events_reopened"]):
        raise ValueError("Six original source ID or exposure history changed before provider preflight")
    source=sel["original_target_metadata_candidates_eligible"]
    prior=exp["original_followup_source_exposure_rows"]
    if [x["original_DOI"] for x in source]!=[x["original_DOI"] for x in prior]:
        raise ValueError("Original source order or DOI changed after exposure audit")
    return source,prior

def inspect(src,prior):
    row={
        "original_source_DOI":src["original_DOI"],
        "original_repository_study_id":src["original_repository_study_id"],
        "original_public_package_uuid":src["original_public_package_uuid"],
        "registered_original_native_labels":src["original_scientific_names"],
        "target_original_native_labels":src["target_matched_scientific_names"],
        "prior_original_source_exposure_status":prior["prior_data_exposure_status"],
        "source_metadata_status":"SOURCE_METADATA_HOLD",
        "source_metadata_hold_reason":None,
        "original_file":None,
        "provider_original_file_list_newly_queried":False,
        "original_GPS_events_opened":0
    }
    if prior["prior_data_exposure_status"]=="KNOWN_PREVIOUS_SOURCE_METADATA_HOLD_NO_GPS_FROM_THAT_ROUTE":
        row["source_metadata_status"]="PREVIOUS_SOURCE_METADATA_HOLD_PRESERVED"
        row["source_metadata_hold_reason"]="PRIOR_ORIGINAL_SOLAR_SOURCE_FILE_AMBIGUITY_HOLD_NOT_REOPENED"
        return row
    if "Homo sapiens" in src["original_scientific_names"]:
        row["source_metadata_status"]="HUMAN_LOCATION_PRIVACY_HOLD"
        row["source_metadata_hold_reason"]="NO_HUMAN_LOCATION_SOURCE_ACCESS"
        return row
    try:
        result=original_files(src["original_public_package_uuid"])
        row["provider_original_file_list_newly_queried"]=True
        original,reason,n,rule=original_primary(result)
        row["source_original_file_candidates_after_filename_rules"]=n
        row["source_original_file_selection_rule"]=rule
        if original is None:
            row["source_metadata_hold_reason"]=reason
        else:
            row["original_file"]=original
            row["source_metadata_status"]="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
    except Exception as exc:
        row["source_metadata_hold_reason"]="PROVIDER_FILE_METADATA_ERROR_"+type(exc).__name__
    return row

def evaluate():
    sources,prior=freeze()
    with ThreadPoolExecutor(max_workers=4) as pool:
        tasks={pool.submit(inspect,s,p):i for i,(s,p) in enumerate(zip(sources,prior))}
        seen={tasks[f]:f.result() for f in as_completed(tasks)}
    rows=[seen[i] for i in range(6)]
    if (len(rows)!=6 or
        [x["original_source_DOI"] for x in rows]!=[x["original_DOI"] for x in sources] or
        sum(x["source_metadata_status"]=="PREVIOUS_SOURCE_METADATA_HOLD_PRESERVED" for x in rows)!=1 or
        sum(x["provider_original_file_list_newly_queried"] for x in rows)>5 or
        any(x["original_GPS_events_opened"]!=0 for x in rows)):
        raise ValueError("Six original source metadata/HOLD denominator changed")
    return {
      "version":VERSION,
      "stage":"AFTER_ALL_SIX_SOURCE_ID_EXPOSURE_FROZEN_BEFORE_TARGET_ORIGINAL_GPS_EVENTS",
      "frozen_six_source_only_identity_git_blob":SIX_SELECTION_BLOB,
      "frozen_prior_exposure_git_blob":EXPOSURE_BLOB,
      "original_followup_source_study_slots":6,
      "prior_original_source_study_ID_exclusions":63,
      "source_original_metadata_status_counts":dict(sorted(Counter(x["source_metadata_status"] for x in rows).items())),
      "original_source_metadata_rows":rows,
      "previous_original_source_GPS_analysis_from_solar_48_unverified_nonexhaustive":True,
      "prior_fox_original_GPS_read_known":True,
      "prior_mallard_source_metadata_HOLD_preserved":True,
      "all_new_original_animal_GPS_events_opened":0,
      "target_selected_after_63_archive_q_outcomes":True,
      "not_global_world_animal_population_sample":True,
      "old_48_source_solar_signed_direction_replicated":False,
      "old_12_ERA5_raw_SHA_HOLD":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z=evaluate()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("SIX_TARGET_ORIGINAL_FILE_METADATA_ONLY",json.dumps({
        "source_statuses":z["source_original_metadata_status_counts"],
        "sources":[{"doi":x["original_source_DOI"],"status":x["source_metadata_status"],
                    "hold":x["source_metadata_hold_reason"],
                    "provider_file_bytes":(x["original_file"] or {}).get("bytes")}
                   for x in z["original_source_metadata_rows"]],
        "original_GPS_opened":0
    },sort_keys=True))

if __name__=="__main__":
    main()
