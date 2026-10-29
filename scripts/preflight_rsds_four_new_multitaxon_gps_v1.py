"""Pre-event source-only ORIGINAL bitstream inventory for four new field studies.

Reads exactly the fixed 437-package Movebank catalog and ORIGINAL provider
bundle metadata. Does NOT read a single original animal GPS event, date,
coordinate, raw track or behavior value. File choice uses filename/status/
size/MD5 only, never sign/magnitude of q or solar-daylight results.
"""
from __future__ import annotations

import argparse
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import HERE, original_files

CONTRACT=HERE/"research/RSDS_FOUR_NEW_MULTINATIVE_SOURCE_GPS_V1_CONTRACT.json"
CATALOG=HERE/"research/seasonal_database_catalog/2026-10-19/movebank/catalog.json"
OLD_EIGHT=HERE/"research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json"
MAX_BYTES=200_000_000


def source_registry():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    original=json.loads(CATALOG.read_text(encoding="utf-8"))
    old=json.loads(OLD_EIGHT.read_text(encoding="utf-8"))
    if (c["version"]!="RSDS_FOUR_NEW_MULTINATIVE_ORIGINAL_GPS_STUDY_TEMPORAL_V1"
            and c["version"]!="RSDS_FOUR_NEW_MULTINATIVE_ORIGINAL_GPS_STUDY_TEMPORAL_V1"):
        raise ValueError("Frozen four-source file/source taxonomy changed")
    if original["package_count"]!=437 or len(original["packages"])!=437 or len(c["sources"])!=4:
        raise ValueError("Original 437 source catalog or four preselected new studies changed")
    by={x["package_uuid"]:x for x in original["packages"]}
    previous={x["source_study_doi"] for x in old["all_sources_including_egret_original_metadata_HOLD"]}
    if len(previous)!=8:
        raise ValueError("Original eight temporal source studies collapsed")
    studies=set(); ids=set();species=set()
    for x in c["sources"]:
        p=by.get(x["package_uuid"])
        if not p:
            raise ValueError("Original prospectively frozen source not in public catalog")
        if (p["dois"]!=["doi:"+x["doi"]] or
            p["repository_study_ids"]!=[x["study_id"]] or
            p["reported_animal_count"]!=x["reported_animals"] or
            p["reported_location_count"]!=x["reported_fixes"] or
            p["scientific_names"]!=x["registered_taxa"] or
            len(set(x["registered_taxa"]))<2 or
            not any("GPS" in v.upper() for v in p["sensor_names"]) or
            p["licenses"]!=["CC0 1.0 Universal"]):
            raise ValueError("Original study DOI/native species/sensor/metadata eligibility inconsistent")
        if (x["doi"] in previous or x["doi"] in studies or
            x["study_id"] in ids):
            raise ValueError("Original study already used in temporal ecology or duplicated among 4 sources")
        studies.add(x["doi"]);ids.add(x["study_id"])
    # The chosen African vultures have another paper DOI sharing the same
    # original Movebank study ID. Never count that extra paper as a second study.
    dupe=next(x for x in original["packages"] if x["package_uuid"]==
              "5935dd02-ffa9-43a7-b091-580d260d89e5")
    if dupe["repository_study_ids"]!=["2919708"]:
        raise ValueError("Original African vulture duplicated study ID reference changed")
    return c,by


def source_only_one(x,catalog):
    report={
        "original_repository_package_uuid":x["package_uuid"],
        "published_DOI":x["doi"],
        "original_study_id":x["study_id"],
        "registered_native_taxa":x["registered_taxa"],
        "published_source_license":catalog["licenses"][0],
        "metadata_reported_animals":x["reported_animals"],
        "metadata_reported_original_fixes":x["reported_fixes"],
        "original_animal_events_opened_before_file_freeze":0,
        "source_original_file":None,
        "source_status":"ORIGINAL_SOURCE_METADATA_HOLD",
    }
    try:
        files=original_files(x["package_uuid"])
        viable=[]
        for z in files:
            name=str(z.get("name") or "")
            low=name.lower()
            if not low.endswith((".csv",".csv.gz",".zip")): continue
            if any(term in low for term in ("reference-data","__macosx","-acc","_acc",
                                             "acceleromet","_meta","readme","gps-dictionary")):
                continue
            if name.split("/")[-1].startswith("._"):continue
            viable.append(z)
        # Unique explicit GPS-tracking ORIGINAL file wins if otherwise
        # ORIGINAL contains companion files; deterministic on filename alone.
        specific=[z for z in viable if "gps" in z["name"].lower()
                  or "tracking" in z["name"].lower()]
        if len(specific)==1:
            viable=specific
        if len(viable)!=1:
            report["hold_reason"]="AMBIGUOUS_OR_MISSING_SINGLE_PRIMARY_TRACK_FILE"
            report["candidate_primary_original_files"]=len(viable)
            return report
        chosen=viable[0]
        if (not isinstance(chosen.get("bytes"),int) or
            chosen["bytes"]<=0 or chosen["bytes"]>MAX_BYTES):
            report["hold_reason"]="FROZEN_200MB_ORIGINAL_FILE_BUDGET_OR_MISSING_SIZE"
            report["candidate_original_size_bytes"]=chosen.get("bytes")
            return report
        checksum=chosen.get("checksum") or {}
        if checksum.get("checkSumAlgorithm","").upper() not in ("MD5","SHA256","SHA-256"):
            report["hold_reason"]="MISSING_OR_UNSUPPORTED_PROVIDER_SOURCE_CHECKSUM"
            return report
        if not checksum.get("value") or not chosen.get("bitstream_uuid"):
            report["hold_reason"]="MISSING_PROVIDER_BITSTREAM_UUID_OR_HASH"
            return report
        report["source_original_file"]={
            "filename":chosen["name"],
            "bitstream_uuid":chosen["bitstream_uuid"],
            "bytes":chosen["bytes"],
            "checksum":checksum
        }
        report["source_status"]="SOURCE_ONLY_ORIGINAL_FILE_ADMISSIBLE"
    except Exception as exc:
        report["hold_reason"]="SOURCE_PROVIDER_METADATA_ERROR_"+type(exc).__name__
    return report


def evaluate():
    c,registry=source_registry()
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(source_only_one,x,registry[x["package_uuid"]]):x["package_uuid"]
                 for x in c["sources"]}
        responses={futures[f]:f.result() for f in as_completed(futures)}
    records=[responses[x["package_uuid"]] for x in c["sources"]]
    return {
        "version":"RSDS_FOUR_NEW_MULTINATIVE_SOURCE_ONLY_ORIGINAL_GPS_V1",
        "date":"2026-10-29",
        "original_published_study_count":4,
        "registered_native_taxon_rows":sum(len(r["registered_native_taxa"]) for r in records),
        "source_only_admissible_original_files":sum(r["source_status"]=="SOURCE_ONLY_ORIGINAL_FILE_ADMISSIBLE" for r in records),
        "source_only_metadata_HOLDS":sum(r["source_status"]=="ORIGINAL_SOURCE_METADATA_HOLD" for r in records),
        "animal_GPS_events_read_by_source_only_preflight":0,
        "source_file_entries":records,
        "sample_was_frozen_before_animal_temporal_q_outcomes":True,
        "original_FOUR_studies_not_four_animal_populations_or_globally_random":True
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    result=evaluate()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_FOUR_NEW_MULTINATIVE_METADATA",json.dumps({
        "studies":result["original_published_study_count"],
        "taxa":result["registered_native_taxon_rows"],
        "admissible":result["source_only_admissible_original_files"],
        "HOLD":result["source_only_metadata_HOLDS"],
        "events":result["animal_GPS_events_read_by_source_only_preflight"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
