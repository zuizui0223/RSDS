"""Original Utah/Wyoming GPS-selected source-context ERA5 numeric identity audit.

Retrieves and verifies original published Movebank GPS solely in ephemeral CI
RAM/disk to RECREATE original fixed top-two q-support source 1-degree-gridyear
contexts. Validates original study MD5, native-taxon q and daylength source
results, then original qday/animal/solar-step support before ANY hourly ERA5
calls. Requests each unique source weather gridyear twice and fingerprints
UTC-hourly PHYSICAL values and serialized JSON separately.

The historical original v0 raw JSON bytes and their physical vectors were NOT
archived: comparison of old physical values is NOT IDENTIFIABLE. This never
substitutes new ERA5 for previously measured animal shortwave effects.
No original gridyear, animal tag, GPS position/date/time or raw weather is
written to the publicly committed study-group source output.
"""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.run_rsds_same_animal_solar_weather_v0 import (
    source_meta, replay_prior_taxon, original_q_gridyears,
    selected_daylight_bearings
)
from scripts.solar_azimuth_era5_wind_v0 import (
    year_request_parameters, get_original_hourly_response
)
from scripts.era5_canonical_physical_content_v1 import (
    compare_two_same_request_responses
)

CONTRACT = HERE / "research/RSDS_ERA5_PHYSICAL_CONTENT_VERSION_AUDIT_V1_CONTRACT.json"
BASELINE = HERE / "research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
STRICT_REPLAY = HERE / "research/RSDS_SAME_TAG_CONTEXT_PAIR_V1_FIRST_RESULT_20261026.json"
NEGATIVE = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def freeze():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    original=json.loads(BASELINE.read_text(encoding="utf-8"))
    earlier=json.loads(STRICT_REPLAY.read_text(encoding="utf-8"))
    if (spec["version"]!="RSDS_ERA5_PHYSICAL_CONTENT_VERSION_AUDIT_V1"
            or original["fixed_original_study_count"]!=2
            or original["total_source_native_taxon_rows"]!=6
            or original["frozen_external_radiation_nonreplication_UNCHANGED"]!=NEGATIVE
            or earlier["native_study_taxon_rows"]!=6
            or earlier["heldout_negative_radiation_nonreplication_UNCHANGED"]!=NEGATIVE):
        raise ValueError("Source-old/new physical weather validation contract inconsistent")
    weather_holds=[c for row in earlier["all_native_source_taxon_results"]
                   for c in row["original_contexts"]]
    if (len(weather_holds)!=12 or
            any(c["status"]!="ERA5_VERSION_HOLD" or
                c.get("ERA5_source_raw_SHA256_replayed") is not False
                for c in weather_holds)):
        raise ValueError("Historic twelve original weather raw JSON holds changed")
    return spec,{s["region"]:s for s in original["sources"]}


def context_guard(now, before):
    known={
        "selected_gridyear_original_coverage_rank":now["rank"],
        "original_qualified_q_days":now["source_q_days"],
        "original_q_qualified_animals":now["source_q_animals"],
        "original_moving_daylight_steps":now["preweather_solar_steps"],
    }
    for k,v in known.items():
        if before.get(k)!=v:
            raise ValueError("Preweather source-original q/gridyear/animal/solar-step eligibility replay changed")
    if not isinstance(before.get("ERA5_original_provider_response_sha256"),str):
        raise ValueError("Previously frozen raw provider SHA not recorded")
    return known


def audit_source(region,out):
    _,prior=freeze()
    source,taxa,refs=source_meta(region)
    if region not in prior or prior[region]["original_source_doi"]!=source["doi"]:
        raise ValueError("Original Movebank field study DOI changed")
    original=prior[region]
    results={
        "version":"RSDS_ERA5_ORIGINAL_GPS_SOURCE_CURRENT_PHYSICAL_CONTENT_AUDIT_V1",
        "date":"2026-10-28", "region":region,
        "source_original_study_DOI":source["doi"],
        "source_original_GPS_MD5_replayed_before_weather":False,
        "all_native_taxon_q_and_daylength_replayed_before_weather":False,
        "native_source_taxa":len(taxa), "taxa":[],
        "status":"ORIGINAL_SOURCE_HOLD",
        "historical_original_12_ERA5_raw_SHA_holds_UNCHANGED":True,
        "previous_source_radiation_external_nonreplication_UNCHANGED":NEGATIVE,
        "historical_v0_physical_content_equivalence":"NOT_IDENTIFIABLE_NO_OLD_HOURLY_NUMERIC_ARRAYS",
        "original_GPS_or_grid_location_republished":False,
        "old_biological_weather_responses_reclassified":False,
    }
    original_tmp=out/("private_original_weather_audit_"+region.lower()+
                      (".zip" if region=="Utah" else ".csv"))
    try:
        candidates=[x for x in original_files(source["package_uuid"])
                    if x["bitstream_uuid"]==source["original_bitstream_uuid"]]
        if (len(candidates)!=1 or candidates[0]["bytes"]!=source["original_bytes"]
                or candidates[0]["checksum"]["value"].lower()!=source["original_md5"].lower()):
            raise ValueError("Original GPS file source UUID/MD5/byte count changed")
        original_verify=checksum_and_download(candidates[0],original_tmp)
        if not original_verify["provider_checksum_verified"]:
            raise ValueError("Original source MD5 failed")
        results["source_original_GPS_MD5_replayed_before_weather"]=True
        if region=="Utah":
            animals, counts, unknown=extract_time_only_original_source(original_tmp,taxa)
        else:
            animals, counts, unknown=streaming_time_selected_csv(original_tmp,taxa)
        if unknown<0:
            raise ValueError("Invalid native source taxon provenance")
        # First pass only: q and q~astronomical L are exactly replayed for
        # ALL original source taxa before a single ERA5 physical request.
        replayed={}
        for taxon in taxa:
            q=replay_prior_taxon(taxon,animals[taxon],refs[taxon])
            oldtaxon=next(z for z in original["taxon_results"] if z["taxon"]==taxon)
            for a,b in (
                ("source_original_q_animals","source_original_q_animals"),
                ("source_original_daylength_animals","source_original_daylength_animals")
            ):
                if q[a]!=oldtaxon[b]:
                    raise ValueError("Old source original animal denominator changed")
            for k in ("source_original_q","source_original_daylength_slope"):
                if abs(q[k]-oldtaxon[k])>1e-12:
                    raise ValueError("Old first study source temporal q biological result changed")
            replayed[taxon]=q
        results["all_native_taxon_q_and_daylength_replayed_before_weather"]=True
        print("RSDS_ORIGINAL_GPS_PREWEATHER_REPLAY",region,len(taxa),flush=True)

        # Cache by PRIVATE source-context grid-year. No private keys, values,
        # original coordinate, GPS year or timestamps survive result writing.
        weather_cache={}
        for taxon in taxa:
            oldtaxon=next(z for z in original["taxon_results"] if z["taxon"]==taxon)
            selected,_=original_q_gridyears(animals[taxon])
            selected=selected_daylight_bearings(animals[taxon],selected)
            if len(selected)!=len(oldtaxon["selected_contexts"]) or len(selected)!=2:
                raise ValueError("Original study q selected source gridyear frame changed")
            row={"taxon":taxon,
                 "original_temporal_q_reproduced":True,
                 "original_selected_contexts":2,
                 "context_results":[]}
            for c,old in zip(selected,oldtaxon["selected_contexts"]):
                context_guard(c,old)
                grid=c["_grid"]
                if grid not in weather_cache:
                    try:
                        params=year_request_parameters(*grid)
                        raw1=get_original_hourly_response(params)
                        time.sleep(1.0)
                        raw2=get_original_hourly_response(params)
                        weather_cache[grid]=(
                            compare_two_same_request_responses(raw1,raw2,grid[2]),
                            None)
                    except Exception as exc:
                        weather_cache[grid]=(None,type(exc).__name__)
                pair,error=weather_cache[grid]
                public={
                    "original_source_q_coverage_rank":c["rank"],
                    "original_q_qualified_source_animals":c["source_q_animals"],
                    "original_qualified_q_days":c["source_q_days"],
                    "original_source_sunlight_moving_steps":c["preweather_solar_steps"],
                    "old_raw_ERA5_source_SHA_still_unverified":True,
                    "historical_old_weather_physical_array_not_archived":True,
                }
                if error is not None:
                    public.update(status="NEW_PROVIDER_PHYSICAL_SOURCE_HOLD",
                                  new_provider_error_type=error)
                else:
                    public.update({
                        "status":"CURRENT_VERSION_TWO_REQUESTS_PHYSICAL_AUDITED",
                        "first_new_raw_SHA_equals_old_raw_SHA":
                            pair["first_new_raw_SHA256"]==
                            old["ERA5_original_provider_response_sha256"],
                        "second_new_raw_SHA_equals_old_raw_SHA":
                            pair["second_new_raw_SHA256"]==
                            old["ERA5_original_provider_response_sha256"],
                        **pair,
                    })
                row["context_results"].append(public)
            results["taxa"].append(row)
            animals[taxon].clear()
        results["unique_source_gridyear_weather_requests"]=len(weather_cache)
        allpairs=[c for t in results["taxa"] for c in t["context_results"]]
        results["source_original_contexts"]=len(allpairs)
        results["source_current_physical_UTC_hourly_same"]=sum(
            c.get("physical_hourly_numeric_SHA_identical") is True for c in allpairs)
        results["source_current_raw_JSON_SHA_changed"]=sum(
            c.get("raw_response_SHA_identical") is False for c in allpairs)
        results["source_old_raw_SHA_first_exact_replays"]=sum(
            c.get("first_new_raw_SHA_equals_old_raw_SHA") is True for c in allpairs)
        results["source_weather_provider_HOLD_contexts"]=sum(
            c["status"]=="NEW_PROVIDER_PHYSICAL_SOURCE_HOLD" for c in allpairs)
        results["status"]="ORIGINAL_GPS_SOURCE_PHYSICAL_VERSION_AUDIT_COMPLETED"
    except Exception as ex:
        results["status"]="ORIGINAL_SOURCE_HOLD"
        results["source_error_type"]=type(ex).__name__
        if not results["taxa"]:
            results["taxa"]=[{"taxon":t,"status":"SOURCE_HOLD"} for t in taxa]
    finally:
        original_tmp.unlink(missing_ok=True)
    out.mkdir(parents=True,exist_ok=True)
    dest=out/("source_"+region.lower()+"_physical_content.json")
    dest.write_text(json.dumps(results,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                    encoding="utf-8")
    print("RSDS_GPS_SOURCE_ERA5_CONTENT_AUDIT",json.dumps({
        "region":region,"status":results["status"],
        "source_taxa":len(results["taxa"]),
        "contexts":results.get("source_original_contexts",0),
        "current_physical_same":results.get("source_current_physical_UTC_hourly_same"),
        "old_original_source_weather_holds_unchanged":True
    },sort_keys=True),flush=True)
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region",required=True,choices=("Utah","Wyoming"))
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    audit_source(a.region,a.out)


if __name__=="__main__":
    main()
