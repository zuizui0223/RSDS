"""Reproducible, outcome-blind source-selection frame from 437 Movebank packages.

Only original public repository METADATA from the immutable 2026-10-19
catalog and frozen historical DOI identity receipts. No new original GPS
rows, timestamps, coordinates, animal IDs or outcomes are accessed.

First-wave selection is a *fixed-snapshot pseudorandom priority sample*
WITHIN 8 single-native and 8 multi-native metadata strata, NOT a world
probability sample of populations. All unselected/duplicate/HOLD records
remain in the original 437-package accounting frame.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from typing import Any

HERE=Path(__file__).resolve().parents[1]
CONTRACT=HERE/"research/RSDS_437_SOURCE_FRAME_WAVE1_V1_CONTRACT.json"
VERSION="RSDS_437_SOURCE_FRAME_WAVE1_V1"

def git_blob(raw: bytes)->str:
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def pinned(path, sha):
    raw=(HERE/path).read_bytes()
    if git_blob(raw)!=sha:
        raise ValueError("Previously fixed source/blob changed: "+path)
    return json.loads(raw)

def normalize_doi(x):
    if not isinstance(x,str):
        raise ValueError("Nonstring source DOI")
    return x.strip().lower().removeprefix("doi:")

def inputs():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (c["version"]!=VERSION or c["original_registry_packages"]!=437
        or c["gps_gnss_sensor_candidates"]!=259
        or c["published_study_DOI_already_outcome_examined"]!=13
        or c["min_reported_animals"]!=10
        or c["min_reported_locations"]!=10000
        or c["wave1_single_native"]!=8 or c["wave1_multi_native"]!=8
        or c["existing_48_signed_shortwave"]!=
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
        or c["existing_ERA5_raw_JSON_12_HOLD"] is not True
        or c["event_level_data_opened"]!=0):
        raise ValueError("Metadata sampling contract changed")
    catalog=pinned(c["source_snapshot"],c["source_catalog_git_blob"])
    manifest=pinned(c["source_manifest"],c["source_manifest_git_blob"])
    previous=pinned(c["previous_12_source_first_result"],c["previous_12_source_git_blob"])
    ill=pinned(c["previous_illinois_result"],c["previous_illinois_git_blob"])
    if (catalog["package_count"]!=437
        or len(catalog["packages"])!=437
        or manifest["reported_total_elements"]!=437
        or manifest["unique_package_uuids"]!=437
        or manifest["page_count"]!=18
        or manifest["raw_pages_sha256_concatenated"]!=
           "86d0be80408b4ebb2531c9554634d5457928f8406e537d05315b40aaee90785d"):
        raise ValueError("Frozen public inventory metadata incomplete")
    old={r["original_study_doi"] for r in
        previous["all_original_study_native_taxon_rows_including_HOLD"]}
    old.add(ill["study_DOI"])
    if len(old)!=13:
        raise ValueError("Previous 12 + Illinois source DOI frame changed")
    if previous["external_48_source_negative_shortwave_direction_nonreplication_UNCHANGED"]!=c["existing_48_signed_shortwave"]:
        raise ValueError("Historical radiation nonreplication altered")
    if previous["original_12_ERA5_raw_JSON_source_HOLDS_UNCHANGED"] is not True:
        raise ValueError("Historical ERA5 HOLD altered")
    if len({x["package_uuid"] for x in catalog["packages"]})!=437:
        raise ValueError("Duplicate catalog package UUID")
    return c,catalog["packages"],old

def standard_date(p):
    # Repository publication date is for deterministic version priority,
    # NEVER interpreted as the biological monitoring date.
    date=(p.get("publication_dates") or [])
    return min(date) if date else "9999-99-99"

def study_id_versions(packages):
    versions=defaultdict(list)
    for p in packages:
        if (set(p.get("sensor_names") or [])&{"GPS","GNSS"}
            and len(p.get("repository_study_ids") or [])==1):
            versions[p["repository_study_ids"][0]].append(p)
    return {id:sorted(v,key=lambda p:(standard_date(p),p["package_uuid"]))
            for id,v in versions.items()}

def priority(salt,catalog_sha,p):
    # Every digit is determined from the source-only contract and catalog.
    # Do not rank by taxon, behavioral response, location or title.
    value="|".join([salt,catalog_sha,p["repository_study_ids"][0],p["package_uuid"]])
    return hashlib.sha256(value.encode("utf-8")).hexdigest()

def census(c,packages,old):
    versions=study_id_versions(packages)
    gps=[p for p in packages if
         set(p.get("sensor_names") or [])&{"GPS","GNSS"}]
    if len(gps)!=259:
        raise ValueError("Expected 259 GPS/GNSS source packages")
    prior_packages=[p for p in packages
                    if any(normalize_doi(d) in old for d in p.get("dois") or [])]
    if len(prior_packages)!=13:
        raise ValueError("Previous DOI package identity missing or ambiguous")
    prior_ids={id for p in prior_packages for id in p.get("repository_study_ids") or []}
    if len(prior_ids)!=13:
        raise ValueError("Previous source study IDs not distinct")

    out=[]
    for p in sorted(packages,key=lambda x:x["package_uuid"]):
        sensors=p.get("sensor_names") or []
        ids=p.get("repository_study_ids") or []
        dois=p.get("dois") or []
        names=p.get("scientific_names") or []
        animals=p.get("reported_animal_count")
        locations=p.get("reported_location_count")
        lic=p.get("licenses") or []
        s="NOT_GPS_GNSS_SENSOR"
        if set(sensors)&{"GPS","GNSS"}:
            if len(ids)!=1 or len(dois)!=1:
                s="SOURCE_ID_OR_DOI_METADATA_HOLD"
            elif ids[0] in prior_ids or normalize_doi(dois[0]) in old:
                s="PREVIOUSLY_OUTCOME_EXAMINED_SOURCE_OR_ALIAS"
            elif versions[ids[0]][0]["package_uuid"]!=p["package_uuid"]:
                s="DUPLICATE_STUDY_ID_NONCANONICAL_VERSION"
            elif lic!=["CC0 1.0 Universal"]:
                s="NON_CC0_LICENSE_MANUAL_REVIEW"
            elif (not isinstance(animals,int) or not isinstance(locations,int)
                  or animals<10 or locations<10000):
                s="INSUFFICIENT_REPORTED_METADATA_PRIORITY_ONLY"
            elif not names or len(set(names))!=len(names) or any(not str(n).strip() for n in names):
                s="MISSING_OR_AMBIGUOUS_TAXON_LABELS"
            else:
                s="FIRST_WAVE_METADATA_CANDIDATE"
        row={
          "package_uuid":p["package_uuid"],
          "doi_labels":dois,
          "repository_study_ids":ids,
          "original_sensors":sensors,
          "original_license_labels":lic,
          "scientific_name_labels":names,
          "reported_animal_count_not_deduplicated":animals,
          "reported_location_count_not_qualified":locations,
          "stage1_status":s,
          "candidate_file_header_or_GPS_event_read":False,
        }
        if s=="FIRST_WAVE_METADATA_CANDIDATE":
            row["selection_stratum"]="multi_native" if len(names)>1 else "single_native"
            row["lottery_sha256"]=priority(c["wave1_rank_salt"],c["source_catalog_git_blob"],p)
        out.append(row)
    candidates=[r for r in out if r["stage1_status"]=="FIRST_WAVE_METADATA_CANDIDATE"]
    ordered={k:sorted([r for r in candidates if r["selection_stratum"]==k],
                      key=lambda r:(r["lottery_sha256"],r["package_uuid"]))
             for k in ("single_native","multi_native")}
    for label,n in (("single_native",8),("multi_native",8)):
        if len(ordered[label])<n:
            raise ValueError("Insufficient response-blind metadata strata for fixed wave")
    chosen=ordered["single_native"][:8]+ordered["multi_native"][:8]
    selected_ids={r["package_uuid"] for r in chosen}
    if len(selected_ids)!=16 or len({r["repository_study_ids"][0] for r in chosen})!=16:
        raise ValueError("Source lottery repeated study IDs")
    for row in out:
        row["wave1_selected"]=row["package_uuid"] in selected_ids
    summary={
      "catalog_date":"2026-10-19",
      "all_original_packages":437,
      "gps_gnss_packages":len(gps),
      "previous_outcome_examined_DOI":len(old),
      "previous_outcome_examined_study_ids":len(prior_ids),
      "metadata_candidate_packages":len(candidates),
      "metadata_candidate_single_native":len(ordered["single_native"]),
      "metadata_candidate_multi_native":len(ordered["multi_native"]),
      "wave1_selected_packages":len(chosen),
      "wave1_selected_single_native":8,
      "wave1_selected_multi_native":8,
      "duplicate_study_ids_among_GPS_packages":
          sum(len(v)>1 for v in versions.values()),
      "original_study_geographic_independence_verified":0,
      "new_animal_GPS_events_or_effects_read":0,
      "status_counts":dict(sorted(Counter(r["stage1_status"] for r in out).items()))
    }
    return {
      "version":VERSION,
      "stage":"FROZEN_SNAPSHOT_METADATA_FRAME_PRE_NEW_ANIMAL_OUTCOME",
      "catalog_git_blob":c["source_catalog_git_blob"],
      "contract_method":"SALT_SHA256_STRATIFIED_METADATA_PRIORITY_NOT_GLOBALLY_RANDOM",
      "summary":summary,
      "all_437_metadata_status_rows":out,
      "wave1_selected_in_rank_order":chosen,
      "old_48_original_source_sunlight_direction_external_NOT_REPLICATED":True,
      "old_12_original_source_ERA5_raw_JSON_HOLD":True,
      "not_a_probability_sample_of_world_wildlife":True
    }

def run():
    c,packages,old=inputs()
    return census(c,packages,old)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",required=True,type=Path)
    a=parser.parse_args()
    result=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_437_META_WAVE1_SCREEN",json.dumps({
        **result["summary"],
        "selected_dois":[r["doi_labels"][0] for r in result["wave1_selected_in_rank_order"]]
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
