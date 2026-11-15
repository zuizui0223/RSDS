"""Source-only original-file provider preflight for 19 frozen second-wave slots.

Original SHA rank source identities were committed before ANY Wave2 file
metadata access. Human mobility source is a privacy HOLD and is never queried.
No source GPS outcomes or exact coordinates are opened. Other 18 source
package file inventories are screened under previous byte/MD5 rules.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned
from scripts.rsds_437_wave2_priority_source_selection_v1 import select, NEG
from scripts.rsds_437_wave1_original_file_metadata_v1 import original_primary
from scripts.movebank_native_ecology_20261019 import original_files

SELECTION_PATH="research/RSDS_437_WAVE2_19_SOURCE_SELECTION_PRE_EVENT_FIRST_V1.json"
SELECTION_GIT_BLOB="516fafc85ade5f953b9f8c697bdcc472d453fd0c"

def freeze():
    doc=pinned(SELECTION_PATH,SELECTION_GIT_BLOB)
    replay=select()
    if (doc!=replay or len(doc["all_original_wave2_sources_in_rank_order"])!=19
        or doc["wave2_new_source_counts"]!={"single_native":16,"multi_native":3}
        or doc["original_wave2_provider_files_opened"]!=0
        or doc["original_wave2_GPS_events_opened"]!=0
        or doc["old_48_signed_sunlight_external_replicated"] is not False
        or doc["old_ERA5_original_12_raw_SHA_HOLD"] is not True):
        raise ValueError("Wave2 original source rank manifest changed after freezing")
    return doc["all_original_wave2_sources_in_rank_order"]

def inspect(s):
    result={
      "wave2_index":s["wave2_index"],
      "package_uuid":s["package_uuid"],
      "original_DOI_labels":s["original_DOI_labels"],
      "original_repository_study_id":s["original_repository_study_id"],
      "registered_native_taxa":s["scientific_name_labels"],
      "selection_stratum":s["stratum"],
      "original_file":None,
      "source_status":"SOURCE_METADATA_HOLD",
      "provider_original_file_list_obtained":False,
      "animal_GPS_events_opened":0
    }
    # An explicitly listed human trajectory is never read, even if public.
    if "Homo sapiens" in s["scientific_name_labels"]:
        result["source_status"]="HUMAN_LOCATION_PRIVACY_HOLD"
        result["source_hold_reason"]="NO_HUMAN_TRAJECTORY_SOURCE_ACCESS"
        return result
    try:
        candidates=original_files(s["package_uuid"])
        result["provider_original_file_list_obtained"]=True
        original,hold,n,rule=original_primary(candidates)
        result["candidate_file_count_after_filename_rules"]=n
        result["file_selection_rule"]=rule
        if original is None:
            result["source_hold_reason"]=hold
        else:
            result["source_status"]="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
            result["original_file"]=original
    except Exception as e:
        result["source_hold_reason"]="PROVIDER_FILE_METADATA_ERROR_"+type(e).__name__
    return result

def evaluate():
    wave=freeze()
    with ThreadPoolExecutor(max_workers=4) as pool:
        future={pool.submit(inspect,x):x["wave2_index"] for x in wave}
        by={future[f]:f.result() for f in as_completed(future)}
    rows=[by[i] for i in range(19)]
    if ([x["original_repository_study_id"] for x in rows]!=[
        x["original_repository_study_id"] for x in wave]
        or sum(x["source_status"]=="HUMAN_LOCATION_PRIVACY_HOLD" for x in rows)!=1
        or any(x["animal_GPS_events_opened"]!=0 for x in rows)):
        raise ValueError("Original fixed DOI/study roster or privacy filtering changed")
    return {
      "version":"RSDS_437_WAVE2_19_ORIGINAL_PROVIDER_FILE_METADATA_FIRST_V1",
      "stage":"AFTER_SHA_SOURCE_SELECTION_BEFORE_ORIGINAL_ANIMAL_GPS_OUTCOMES",
      "wave2_source_selection_git_blob":SELECTION_GIT_BLOB,
      "source_study_slots":19,
      "new_wave2_single_native":16,
      "new_wave2_multi_native":3,
      "source_status_counts":dict(sorted(Counter(x["source_status"] for x in rows).items())),
      "source_rows":rows,
      "provider_file_metadata_only":True,
      "human_GPS_or_position_data_retrieved":False,
      "all_new_wave2_animal_GPS_events_opened":0,
      "original_48_source_radiation_negative_external_replication_UNCHANGED":NEG,
      "old_12_ERA5_original_raw_SHA_HOLD_UNCHANGED":True,
      "limits":"File names, UUIDs and MD5 hints only; original GPS file bytes/schema/animal cadence/species taxonomy not yet inspected. Preserve all provider, budget, and human-location HOLDS."
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    x=evaluate()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_WAVE2_ORIGINAL_PROVIDER_METADATA",json.dumps({
      "statuses":x["source_status_counts"],
      "sources":[{"i":z["wave2_index"],"doi":z["original_DOI_labels"],
                  "status":z["source_status"],"hold":z.get("source_hold_reason")}
                 for z in x["source_rows"]],
      "original_GPS_events":0,
      "human_GPS_read":False
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
