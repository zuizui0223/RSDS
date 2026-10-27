"""Two independent same-request ERA5 hourly responses at public NO-GPS sites.

Purpose: determine whether provider raw serialized SHA instability implies a
physical hourly weather-value change. These two fixed *synthetic benchmark*
geographies do NOT reveal original GPS locations and are NOT used to fit
animal-behavior effects. The old twelve animal-source raw JSON mismatches
remain HOLD regardless of this test.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import HERE
from scripts.solar_azimuth_era5_wind_v0 import (
    get_original_hourly_response, year_request_parameters
)
from scripts.era5_canonical_physical_content_v1 import (
    compare_two_same_request_responses
)

CONTRACT=HERE/"research/RSDS_ERA5_PHYSICAL_CONTENT_VERSION_AUDIT_V1_CONTRACT.json"


def run():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c["version"] != "RSDS_ERA5_PHYSICAL_CONTENT_VERSION_AUDIT_V1" or \
            c["prior_old_raw_12"]["old_raw_JSON_SHA_mismatch_contexts"] != 12:
        raise ValueError("Old ERA5 first-result source context frame changed")
    rows=[]
    for meta in c["two_strata"]["public_nonGPS_benchmark"]:
        result={
            "benchmark_name":meta["name"],
            "source_type":"PRECHOSEN_PUBLIC_NO_ORIGINAL_ANIMAL_GPS_SITE",
            "year":meta["year"],
            "status":"BENCHMARK_PROVIDER_HOLD",
        }
        try:
            params=year_request_parameters(
                meta["latitude_floor"],meta["longitude_floor"],meta["year"])
            raw1=get_original_hourly_response(params)
            time.sleep(1)
            raw2=get_original_hourly_response(params)
            result.update(compare_two_same_request_responses(raw1,raw2,meta["year"]))
            result["status"]="BENCHMARK_TWO_REQUESTS_VALIDATED"
        except Exception as e:
            result["provider_error_type"]=type(e).__name__
        rows.append(result)
    return {
        "version":"RSDS_ERA5_UTC_HOURLY_CANONICAL_BENCHMARK_V1_FIRST_RESULT",
        "date":"2026-10-28",
        "benchmarks":rows,
        "public_no_GPS_source_benchmarks":len(rows),
        "two_valid_requests_per_site_completed":sum(
            r["status"]=="BENCHMARK_TWO_REQUESTS_VALIDATED" for r in rows),
        "raw_JSON_SHA_disagreements":sum(
            r.get("raw_response_SHA_identical") is False for r in rows),
        "stable_hourly_physical_content_despite_raw_JSON_mismatch":sum(
            r.get("raw_SHA_disagreement_with_identical_physical_values_and_grid") is True
            for r in rows),
        "unarchived_old_12_original_GPS_context_physical_data_equivalence":
            "NOT_IDENTIFIABLE_WITH_ONLY_OLD_RAW_RESPONSE_SHA",
        "old_12_original_source_ERA5_version_HOLDS_UNCHANGED":True,
        "old_48_source_external_signed_radiation_NONREPLICATED_UNCHANGED":
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE",
        "claims":[
            "This benchmark tests the reproducibility of NOW-to-NOW physical hourly numeric values, not OLD-original-versus-NOW source climate values.",
            "An identical canonical hash may coexist with a different raw JSON SHA because JSON metadata/formatting differ. Grid metadata is flagged separately.",
            "Benchmark sites are source-independent fixed coordinates, not animal source study sites, and no animal behavior or migration effect is estimated.",
            "These results cannot turn the old 12 GPS-linked weather SHA HOLDS into supported biological weather-replication effects."
        ],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,default=HERE/"research/RSDS_ERA5_CANONICAL_BENCHMARK_FIRST_RESULT_20261028.json")
    a=p.parse_args()
    d=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2,ensure_ascii=False,sort_keys=True)+"\n",
                     encoding="utf-8")
    print("RSDS_ERA5_PHYSICAL_BENCHMARK",json.dumps({
        "valid":d["two_valid_requests_per_site_completed"],
        "raw_sha_diff":d["raw_JSON_SHA_disagreements"],
        "semantic_equal_raw_changed":d["stable_hourly_physical_content_despite_raw_JSON_mismatch"],
        "old_v0_physical_comparison":d["unarchived_old_12_original_GPS_context_physical_data_equivalence"]
    },sort_keys=True),flush=True)


if __name__ == "__main__":
    main()
