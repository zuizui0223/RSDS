"""Summarize real original GPS-selected ERA5 canonical physical audits by study.

Study + native taxon are the only output units. Prior historical raw JSON
SHA disagreements and unarchived physical arrays remain HOLD, never biological
nulls or confirmed old-now weather matching.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import HERE

ROOT=HERE
CONTRACT=HERE/"research/RSDS_ERA5_PHYSICAL_CONTENT_VERSION_AUDIT_V1_CONTRACT.json"
BASELINE=HERE/"research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
HOLD=HERE/"research/RSDS_SAME_TAG_CONTEXT_PAIR_V1_FIRST_RESULT_20261026.json"
OUT=HERE/"research/RSDS_ERA5_ORIGINAL_SOURCES_CANONICAL_FIRST_RESULT_20261028.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def aggregate(utah,wyoming):
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    base=json.loads(BASELINE.read_text(encoding="utf-8"))
    historic=json.loads(HOLD.read_text(encoding="utf-8"))
    studies={"Utah":utah,"Wyoming":wyoming}
    reference={x["region"]:x for x in base["sources"]}
    if (spec["version"]!="RSDS_ERA5_PHYSICAL_CONTENT_VERSION_AUDIT_V1"
            or historic["native_study_taxon_rows"]!=6
            or base["fixed_original_study_count"]!=2
            or historic["heldout_negative_radiation_nonreplication_UNCHANGED"]!=NEGATIVE):
        raise ValueError("Original source version audit denominator changed")
    full=[]
    cats=[]
    records=[]
    for region in ("Utah","Wyoming"):
        s=studies[region]
        if (s.get("region")!=region or s.get("source_original_study_DOI")!=reference[region]["original_source_doi"]
                or s.get("version")!="RSDS_ERA5_ORIGINAL_GPS_SOURCE_CURRENT_PHYSICAL_CONTENT_AUDIT_V1"
                or s.get("previous_source_radiation_external_nonreplication_UNCHANGED")!=NEGATIVE
                or s.get("old_biological_weather_responses_reclassified") is not False):
            raise ValueError("Invalid original source study DOI, version or old nonreplication")
        expected={t["taxon"]:t for t in reference[region]["taxon_results"]}
        found=s["taxa"]
        if len(found)!=len(expected) or {x["taxon"] for x in found}!=set(expected):
            raise ValueError("Original native study taxonomy not retained")
        if (s["status"]!="ORIGINAL_SOURCE_HOLD" and
                (not s.get("source_original_GPS_MD5_replayed_before_weather")
                or not s.get("all_native_taxon_q_and_daylength_replayed_before_weather"))):
            raise ValueError("Weather was opened before original GPS biology source checks")
        records.append({
            "region":region,"source_DOI":s["source_original_study_DOI"],
            "state":s["status"],
            "GPS_MD5_replayed":s["source_original_GPS_MD5_replayed_before_weather"],
            "original_native_temporal_results_replayed_before_weather":
                s["all_native_taxon_q_and_daylength_replayed_before_weather"],
            "native_taxon_count":len(found),
            "unique_physical_weather_source_queries":s.get("unique_source_gridyear_weather_requests"),
        })
        for x in found:
            contexts=x.get("context_results",[])
            if x.get("original_selected_contexts")==2 and len(contexts)!=2:
                raise ValueError("Source contexts vanish despite original q/coverage replay")
            if s["status"]!="ORIGINAL_SOURCE_HOLD" and (not x.get("original_temporal_q_reproduced") or len(contexts)!=2):
                raise ValueError("The original source q and two selected contexts did not replay")
            entry={"region":region,"taxon":x["taxon"],
                   "source_original_weather_contexts":len(contexts),
                   "previous_v0_physical_content_provable":False,
                   "context_outcomes":contexts,
                   "status":"SOURCE_HOLD" if not contexts else
                            ("CURRENT_TWO_REQUESTS_PHYSICAL_AUDITED"
                             if all(c["status"]=="CURRENT_VERSION_TWO_REQUESTS_PHYSICAL_AUDITED" for c in contexts)
                             else "SOME_CONTEXT_PROVIDER_HOLD")}
            full.append(entry)
            cats.extend(contexts)
    if len(records)!=2 or len(full)!=6:
        raise ValueError("Not the two original GPS field studies/six native taxa")
    grouped=Counter(c.get("status","SOURCE_HOLD") for c in cats)
    comparisons=[c for c in cats if c.get("status")=="CURRENT_VERSION_TWO_REQUESTS_PHYSICAL_AUDITED"]
    return {
        "version":"RSDS_ORIGINAL_GPS_CONTEXT_CURRENT_ERA5_PHYSICAL_V1_RESULT",
        "date":"2026-10-28",
        "status":"NEW_PROVIDER_HOURLY_CONTENT_AUDIT_NOT_OLD_PHYSICAL_REPLICATION",
        "original_studies":2,
        "original_study_taxon_groups":6,
        "original_selected_weather_contexts_expected":12,
        "original_selected_weather_contexts_reconstructed":len(cats),
        "current_provider_contexts_two_response_comparisons":len(comparisons),
        "context_provider_result_status_counts":dict(grouped),
        "within_current_session_raw_JSON_sha_mismatch":sum(
            c.get("raw_response_SHA_identical") is False for c in comparisons),
        "within_current_session_canonical_UTC_hourly_values_equal":sum(
            c.get("physical_hourly_numeric_SHA_identical") is True for c in comparisons),
        "within_current_session_raw_json_diff_but_physical_and_grid_equal":sum(
            c.get("raw_SHA_disagreement_with_identical_physical_values_and_grid") is True
            for c in comparisons),
        "CURRENT_raw_SHA_matches_historic_v0_raw_SHA_first_response":sum(
            c.get("first_new_raw_SHA_equals_old_raw_SHA") is True for c in comparisons),
        "baseline_raw_JSON_SHA_historical_replay_HOLD_contexts":12,
        "baseline_physical_content_OLD_vs_NEW":"NOT_IDENTIFIABLE_OLD_HOURLY_NUMERIC_VALUES_NOT_ARCHIVED",
        "new_original_source_weather_numeric_snapshot_does_not_change_old_biological_effects":True,
        "source_study_results":records,
        "full_original_taxon_context_results":full,
        "historical_original_12_source_version_HOLDs_UNCHANGED":True,
        "original_48_study_external_negative_shortwave_replication_UNCHANGED":NEGATIVE,
        "interpretation":[
            "Two current provider requests at the same GPS-selected gridyear test current physical-data reproducibility, not the unarchived old numeric arrays.",
            "A raw JSON SHA mismatch can reflect transient response-generation metadata, formatting or numbers, so only independently frozen canonical hourly numeric hashes can establish physical equality for NEW versions.",
            "No GPS exact coordinates, grid identifiers, calendar years, timestamps, tag IDs, individual tracking records or original weather arrays are published.",
            "The old twelve v0-v1 raw SHA differences are real as file identities but OLD-vs-NEW physical equivalence is unknowable from SHA fingerprints alone.",
            "The previous source-specific animal effects and failed external solar radiation replication remain untouched; no causal ERA5 wildlife response is inferred."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input",type=Path,required=True)
    cli.add_argument("--out",type=Path,default=OUT)
    args=cli.parse_args()
    files={}
    for region in ("utah","wyoming"):
        f=args.input/("source_"+region+"_physical_content.json")
        if not f.is_file():
            raise ValueError("Previously frozen original-source region missing: "+region)
        files[region]=json.loads(f.read_text(encoding="utf-8"))
    d=aggregate(files["utah"],files["wyoming"])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print("RSDS_ORIGINAL_ERA5_CANONICAL_SOURCE_VERSION",json.dumps({
        "studies":d["original_studies"],"native_taxa":d["original_study_taxon_groups"],
        "selected_contexts":d["original_selected_weather_contexts_reconstructed"],
        "provider_physical_equal":d["within_current_session_canonical_UTC_hourly_values_equal"],
        "raw_SHA_now_mismatched":d["within_current_session_raw_JSON_sha_mismatch"],
        "original_historical_holds":d["baseline_raw_JSON_SHA_historical_replay_HOLD_contexts"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
