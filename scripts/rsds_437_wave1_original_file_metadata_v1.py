"""Pre-event ORIGINAL FILE metadata audit for the SHA-frozen 16 source wave.

Reads public Movebank provider ORIGINAL bundle listings only; never downloads,
streams, parses or reads original GPS events/animal IDs/coordinates, and does
not inspect metadata source titles or ecological responses for file selection.
Stops with source-specific HOLD if ambiguous, oversized, missing checksum,
or provider metadata fail. All 16 slots, including HOLD, must be retained.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import original_files
from scripts.rsds_437_systematic_source_frame_v1 import (
    HERE, git_blob, run as build_frame
)

PINNED_SELECTION="research/RSDS_437_SOURCE_FRAME_WAVE1_PRE_EVENT_SELECTION_20261115.json"
SELECTION_BLOB_SHA="23d7489e4cf8347341cf67290499f278610651c6"
FRAME_SHA256="927e7bbd78d3cf7a5e0a6cfe144da6f47b173d4e543e08950b241a7fd1d34960"
SELECTED_SHA256="28c1eaa93c5ee903a65079bf20a6b362aa3784c04f0152575e41712bc77710e5"
MAX_BYTES=200_000_000

def freeze():
    raw=(HERE/PINNED_SELECTION).read_bytes()
    if git_blob(raw)!=SELECTION_BLOB_SHA:
        raise ValueError("Pre-event selected source manifest was modified")
    manifest=json.loads(raw)
    if (manifest["version"]!="RSDS_FIXED_WAVE1_437_SNAPSHOT_SOURCE_ONLY_RESULT_V1"
        or manifest["original_full_437_row_first_result_SHA256"]!=FRAME_SHA256
        or manifest["original_selected_16_normalized_row_SHA256"]!=SELECTED_SHA256
        or manifest["original_registry_package_count"]!=437
        or manifest["original_GPS_GNSS_candidate_count"]!=259
        or manifest["meta_stage1_candidate_count"]!=120
        or manifest["meta_stage1_single_candidate_count"]!=109
        or manifest["meta_stage1_multi_candidate_count"]!=11
        or manifest["animal_event_records_accessed"]!=0
        or manifest["old_48_signed_radiation_external_nonreplication"]!=
           "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
        or manifest["old_ERA5_raw_SHA_HOLD"] is not True):
        raise ValueError("Pre-event wave selection/censoring frame changed")
    frame=build_frame()
    serialized=json.dumps(frame,ensure_ascii=False,sort_keys=True,indent=2)+"\n"
    if hashlib.sha256(serialized.encode("utf-8")).hexdigest()!=FRAME_SHA256:
        raise ValueError("Replayed 437-package first source-frame SHA differs")
    selected=frame["wave1_selected_in_rank_order"]
    check=json.dumps(selected,sort_keys=True,separators=(",",":")).encode()
    if hashlib.sha256(check).hexdigest()!=SELECTED_SHA256:
        raise ValueError("Source rank or selected metadata changed")
    pinned=manifest["pre_event_selected_sources_in_rank_order"]
    if len(pinned)!=16 or len(selected)!=16:
        raise ValueError("Source wave is no longer sixteen slots")
    if len({r["original_repository_study_id"] for r in pinned})!=16:
        raise ValueError("Study ID duplicated in provider source wave")
    for i,(a,b) in enumerate(zip(pinned,selected),1):
        if (a["ordinal"]!=i or a["package_uuid"]!=b["package_uuid"]
            or a["source_DOI_label"]!=b["doi_labels"][0]
            or a["original_repository_study_id"]!=b["repository_study_ids"][0]
            or a["scientific_names"]!=b["scientific_name_labels"]
            or a["selection_stratum"]!=b["selection_stratum"]
            or b["stage1_status"]!="FIRST_WAVE_METADATA_CANDIDATE"):
            raise ValueError("Provider preflight attempted outcome-conditioned source replacement")
    return manifest,pinned

def original_primary(files):
    """Unique source-only candidate, subject to fixed filename/size rules."""
    viable=[]
    for x in files:
        name=str(x.get("name") or "")
        low=name.lower()
        if (not low.endswith((".csv",".csv.gz",".zip"))
            or any(w in low for w in ("reference-data","__macosx",
                                     "acceleromet","_acc","-acc","readme",
                                     "metadata","_meta","gps-dictionary"))
            or name.split("/")[-1].startswith("._")):
            continue
        viable.append(x)
    specific=[z for z in viable if "gps" in z["name"].lower()
              or "tracking" in z["name"].lower()]
    rule="ONE_ONLY_PRIMARY_SOURCE_FILE"
    if len(specific)==1:
        viable=specific
        rule="UNIQUE_GPS_OR_TRACKING_FILENAME_ONLY"
    if len(viable)==2 and all(z["name"].lower().endswith(".csv") for z in viable):
        large,small=sorted(viable,key=lambda z:z.get("bytes") or 0,reverse=True)
        if (isinstance(large.get("bytes"),int) and
            isinstance(small.get("bytes"),int)
            and large["bytes"]>=10_000_000
            and small["bytes"]>0
            and large["bytes"]>=10*small["bytes"]):
            viable=[large]
            rule="LARGEST_CSV_10X_SMALLER_COMPANION_ONLY"
    if len(viable)!=1:
        return None,"AMBIGUOUS_OR_MISSING_PRIMARY_SOURCE_FILE",len(viable),rule
    one=viable[0]
    n=one.get("bytes")
    if not isinstance(n,int) or not 0<n<=MAX_BYTES:
        return None,"FROZEN_200MB_SOURCE_FILE_BUDGET_OR_UNKNOWN_BYTES",1,rule
    c=one.get("checksum") or {}
    if (str(c.get("checkSumAlgorithm") or "").upper()!="MD5"
        or not isinstance(c.get("value"),str) or len(c["value"])!=32
        or any(z not in "0123456789abcdef" for z in c["value"].lower())
        or not one.get("bitstream_uuid")):
        return None,"SOURCE_PROVIDER_MD5_OR_BITSTREAM_UUID_MISSING",1,rule
    return {"filename":one["name"],"bitstream_uuid":one["bitstream_uuid"],
            "bytes":n,"checksum":{"algorithm":"MD5","value":c["value"].lower()}},None,1,rule

def inspect(row):
    report={
        "ordinal":row["ordinal"],
        "package_uuid":row["package_uuid"],
        "source_DOI_label":row["source_DOI_label"],
        "repository_study_id":row["original_repository_study_id"],
        "selection_stratum":row["selection_stratum"],
        "registered_native_labels":row["scientific_names"],
        "source_status":"SOURCE_METADATA_HOLD",
        "original_file":None,
        "provider_original_metadata_listed":False,
        "animal_event_records_opened":0,
    }
    try:
        original=original_files(row["package_uuid"])
        report["provider_original_metadata_listed"]=True
        file,hold,n,rule=original_primary(original)
        report["metadata_candidate_count_after_name_rules"]=n
        report["original_file_selection_rule"]=rule
        if file is not None:
            report["original_file"]=file
            report["source_status"]="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
        else:
            report["hold_reason"]=hold
    except Exception as e:
        report["hold_reason"]="SOURCE_PROVIDER_FILE_METADATA_ERROR_"+type(e).__name__
    return report

def evaluate():
    manifest,wave=freeze()
    from concurrent.futures import ThreadPoolExecutor,as_completed
    with ThreadPoolExecutor(max_workers=4) as pool:
        todo={pool.submit(inspect,r):r["ordinal"] for r in wave}
        rows={todo[f]:f.result() for f in as_completed(todo)}
    records=[rows[i] for i in range(1,17)]
    return {
        "version":"RSDS_WAVE1_16_ORIGINAL_FILE_METADATA_PREFLIGHT_V1",
        "stage":"POST_PRE_EVENT_SELECTION_BLOB_FROZEN_NO_GPS_EVENTS",
        "selection_manifest_git_blob":SELECTION_BLOB_SHA,
        "full_437_source_frame_sha256":FRAME_SHA256,
        "source_snapshot_date":"2026-10-19",
        "original_package_slots":16,
        "selected_single_native_slots":8,
        "selected_multi_native_slots":8,
        "source_status_counts":dict(sorted(Counter(r["source_status"] for r in records).items())),
        "provider_original_source_rows":[r for r in records],
        "original_GPS_event_records_opened":0,
        "prior_48_shortwave_direction_external_not_replicated":True,
        "prior_original_ERA5_12_source_HOLD":True,
        "scientific_claim_ceiling":"Source MD5/UUID/bytes hint only. Not a GPS schema, valid native taxonomy, temporal cadence, independent biological population, admissible ecology or global probability sample."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    result=evaluate()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                        encoding="utf-8")
    print("RSDS_WAVE1_SOURCE_METADATA",json.dumps({
        "total":result["original_package_slots"],
        "statuses":result["source_status_counts"],
        "studies":[
            {"i":r["ordinal"],"DOI":r["source_DOI_label"],
             "status":r["source_status"],"reason":r.get("hold_reason"),
             "bytes":(r["original_file"] or {}).get("bytes")}
            for r in result["provider_original_source_rows"]
        ],
        "GPS_event_rows":result["original_GPS_event_records_opened"]
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
