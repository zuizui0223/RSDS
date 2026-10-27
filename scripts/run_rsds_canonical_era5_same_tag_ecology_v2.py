"""Re-analyze SAME original tagged animal GPS q and Sun-heading shortwave effects.

New semantic ERA5 version only: verifies original Movebank provider MD5 and all
native study q/daylight results first, exact old q-selected original contexts,
then checks each current ERA5 full UTC-hourly PHYSICAL fingerprint matches the
NEW source-locked 12-context first result. Never silently substitutes an old
serialized JSON response or revises old response claims.

Uses the EXISTING fixed old v0 within animal,month,1-hour apparent clock +
solar altitude, azimuth, temperature and wind partial-regression estimator.
Keeps all original native taxa and response/source holds. Private original
GPS coordinates, tags, times and individual estimates remain ephemeral.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.run_rsds_same_animal_solar_weather_v0 import (
    source_meta, replay_prior_taxon, original_q_gridyears,
    selected_daylight_bearings, evaluate_context, boot_summary
)
from scripts.run_rsds_same_tag_context_pair_v1 import (
    paired_context_effects
)
from scripts.era5_canonical_physical_content_v1 import physical_identity
from scripts.solar_azimuth_era5_wind_v0 import (
    year_request_parameters, get_original_hourly_response,
    validate_and_parse_weather
)

CONTRACT=HERE/"research/RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_CONTRACT.json"
FIRST=HERE/"research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
AUDIT=HERE/"research/RSDS_ERA5_ORIGINAL_SOURCES_CANONICAL_FIRST_RESULT_20261028.json"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
AXES=("sun_cos","sun_sin","q")


def frozen():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    def lock(path,expected):
        raw=path.read_bytes()
        sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
        if sha != expected:
            raise ValueError("Original source biological or new physical version SHA changed")
        return json.loads(raw.decode("utf-8"))
    legacy=lock(FIRST,spec["historic_first_biological_result_git_blob_sha"])
    audit=lock(AUDIT,spec["parent_version_audit_git_blob_sha"])
    if spec["version"]!="RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2":
        raise ValueError("The frozen canonical version scientific route has changed")
    if (legacy["fixed_original_study_count"]!=2
            or legacy["total_source_native_taxon_rows"]!=6
            or audit["original_studies"]!=2 or audit["original_study_taxon_groups"]!=6
            or audit["original_selected_weather_contexts_reconstructed"]!=12
            or audit["current_provider_contexts_two_response_comparisons"]!=12
            or audit["within_current_session_canonical_UTC_hourly_values_equal"]!=12
            or audit["original_48_study_external_negative_shortwave_replication_UNCHANGED"]!=NEG
            or legacy["frozen_external_radiation_nonreplication_UNCHANGED"]!=NEG
            or audit["baseline_raw_JSON_SHA_historical_replay_HOLD_contexts"]!=12):
        raise ValueError("First source numeric version frame or negative old outcome changed")
    return spec,{x["region"]:x for x in legacy["sources"]},audit


def compare_context_eligibility(selected, first, audit):
    if (selected["rank"]!=first["selected_gridyear_original_coverage_rank"]
            or selected["source_q_days"]!=first["original_qualified_q_days"]
            or selected["source_q_animals"]!=first["original_q_qualified_animals"]
            or selected["preweather_solar_steps"]!=first["original_moving_daylight_steps"]):
        raise ValueError("Original q-day, animal and daylight direction eligible context changed")
    if selected["rank"]!=audit["original_source_q_coverage_rank"] or \
       audit["status"]!="CURRENT_VERSION_TWO_REQUESTS_PHYSICAL_AUDITED":
        raise ValueError("Original provider current UTC-hourly physical version not previously frozen")
    if (audit["original_qualified_q_days"]!=selected["source_q_days"]
            or audit["original_q_qualified_source_animals"]!=selected["source_q_animals"]
            or audit["original_source_sunlight_moving_steps"]!=selected["preweather_solar_steps"]
            or not audit["physical_hourly_numeric_SHA_identical"]
            or not audit["returned_weather_grid_metadata_identical"]
            or audit["historical_old_weather_physical_array_not_archived"] is not True):
        raise ValueError("Original source canonical climate/animal preflight changed")


def analyse_taxon_weather(taxon, tracks, old, audited, cache):
    contexts,_=original_q_gridyears(tracks)
    contexts=selected_daylight_bearings(tracks,contexts)
    if len(contexts)!=2 or len(old["selected_contexts"])!=2 or len(audited["context_outcomes"])!=2:
        raise ValueError("Original native taxon has wrong two source context support")
    out={"taxon":taxon,"original_selected_source_contexts":2,
         "original_source_native_q_replay":True,
         "new_current_provider_version": "SEPARATE_CANONICAL_PHYSICAL_V2",
         "previous_raw_ERA5_version_HOLD_UNCHANGED":True,
         "context_results":[],"source_animal_aggregate":{},
         "within_same_tag_rank2_minus_rank1":None}
    all_tag=defaultdict(lambda:defaultdict(list))
    perctx=[]
    for c,orig,current in zip(contexts,old["selected_contexts"],
                              audited["context_outcomes"]):
        compare_context_eligibility(c,orig,current)
        key=c["_grid"]
        if key not in cache:
            try:
                raw=get_original_hourly_response(year_request_parameters(*key))
                newphysical=physical_identity(raw,key[2])
                reference=current["first_new_hourly_physical_SHA256"]
                if newphysical["physical_UTC_hourly_SHA256"]!=reference:
                    raise ValueError("CURRENT_SOURCE_NUMERIC_VERSION_CHANGED")
                weather,receipt=validate_and_parse_weather(raw,key[2])
                cache[key]=(weather,"CANONICAL_PHYSICAL_EXACT_SHA_MATCH",newphysical)
            except Exception as exc:
                cache[key]=(None,type(exc).__name__,None)
        weather,state,numeric=cache[key]
        row={
            "original_q_coverage_context_rank":c["rank"],
            "old_source_q_qualified_days":c["source_q_days"],
            "old_source_q_animals":c["source_q_animals"],
            "old_sun_moving_steps":c["preweather_solar_steps"],
            "new_provider_physical_content_status":state,
            "current_hourly_numeric_SHA_expected_from_pre_estimator_audit":
                current["first_new_hourly_physical_SHA256"],
            "original_raw_ERA5_v0_physical_equivalence":
                "NOT_IDENTIFIABLE_NO_OLD_NUMERIC_ARRAYS",
            "new_source_effects":None,
            "original_first_v0_effects_REFERENCE_ONLY":{
                k:orig.get(k,{}).get("source_mean")
                for k in ("direction_sun_cos","direction_sun_sin","daypart_q_vs_irradiance")
            }
        }
        if weather is None:
            row["status"]="NEW_CANONICAL_PROVIDER_NUMERIC_HOLD"
            perctx.append(None)
        else:
            estimate,effects=evaluate_context(c,weather)
            row["new_source_effects"]={
                "source_state":estimate["ecology_state"],
                "sun_cos":estimate.get("direction_sun_cos"),
                "sun_sin":estimate.get("direction_sun_sin"),
                "daypart_q":estimate.get("daypart_q_vs_irradiance"),
                "eligible_direction_source_animals":
                    estimate.get("direction_animals_with_identifiable_slope",0),
                "eligible_joint_q_and_sun_source_animals":
                    estimate.get("joint_daypart_and_direction_animals",0),
                "ERA5_GPS_hourly_link_fraction":
                    estimate.get("hourly_ERA5_movement_join_fraction")
            }
            row["status"]="NEW_CANONICAL_WEATHER_BIOLOGY_EVALUATED"
            perctx.append(effects)
            for tag,axes in effects.items():
                for field in AXES:
                    if axes.get(field) is not None:
                        all_tag[tag][field].append(axes[field])
        out["context_results"].append(row)
    agg={
        tag:{axis:mean(values) for axis,values in byax.items()}
        for tag,byax in all_tag.items()
    }
    out["source_animal_aggregate"]={
        "daypart_q_vs_irradiance":boot_summary([r["q"] for r in agg.values() if "q" in r]),
        "sun_cos_vs_irradiance":boot_summary([r["sun_cos"] for r in agg.values() if "sun_cos" in r]),
        "sun_sin_vs_irradiance":boot_summary([r["sun_sin"] for r in agg.values() if "sun_sin" in r]),
        "original_tag_ids_not_exported":True
    }
    if len(perctx)==2 and all(x is not None for x in perctx):
        out["within_same_tag_rank2_minus_rank1"]=paired_context_effects(
            perctx[0],perctx[1])
    out["status"]="CURRENT_VERSION_SOURCE_CONTEXTS_COMPLETE" \
        if all(x["status"]=="NEW_CANONICAL_WEATHER_BIOLOGY_EVALUATED" for x in out["context_results"]) \
        else "CANONICAL_WEATHER_HOLD"
    return out


def run_source(region,out_dir):
    spec,first,audit=frozen()
    regions={"Utah":spec["study_coverage"]["utah"],
             "Wyoming":spec["study_coverage"]["wyoming"]}
    source,taxa,reference=source_meta(region)
    if (region not in regions or source["doi"]!=regions[region]["doi"]
            or taxa!=regions[region]["native_taxa"]
            or first[region]["original_source_doi"]!=source["doi"]):
        raise ValueError("Original Movebank source region native taxonomy changed")
    audit_rows={(r["region"],r["taxon"]):r for r in audit["full_original_taxon_context_results"]}
    previous={r["taxon"]:r for r in first[region]["taxon_results"]}
    output={
        "version":"RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_ORIGINAL_SOURCE_RESULT",
        "region":region,"source_study_DOI":source["doi"],
        "original_q_and_daylength_replayed_before_weather":False,
        "original_GPS_provider_MD5_replayed":False,
        "source_native_taxa":list(taxa),"taxa":[],
        "status":"ORIGINAL_SOURCE_HOLD",
        "separate_new_physical_weather_version":True,
        "previous_ERA5_raw_SHA_12_context_HOLDS_UNCHANGED":True,
        "prior_48_source_radiation_external_NONREPLICATED_UNCHANGED":NEG,
        "new_ERA5_shortwave_effect_is_NOT_replication_of_old_weather_dataset":True
    }
    private=out_dir/("private_original_new_canonical_"+region.lower()+
                     (".zip" if region=="Utah" else ".csv"))
    try:
        native=[p for p in original_files(source["package_uuid"])
                if p["bitstream_uuid"]==source["original_bitstream_uuid"]]
        if (len(native)!=1 or native[0]["bytes"]!=source["original_bytes"]
                or native[0]["checksum"]["value"].lower()!=source["original_md5"].lower()):
            raise ValueError("Original source native GPS bitstream MD5/size changed")
        verification=checksum_and_download(native[0],private)
        if not verification["provider_checksum_verified"]:
            raise ValueError("Provider original source MD5 was not verified")
        output["original_GPS_provider_MD5_replayed"]=True
        if region=="Utah":
            tracks,cnts,other=extract_time_only_original_source(private,taxa)
        else:
            tracks,cnts,other=streaming_time_selected_csv(private,taxa)
        for taxon in taxa:
            saved=replay_prior_taxon(taxon,tracks[taxon],reference[taxon])
            old=previous[taxon]
            if (saved["source_original_q_animals"]!=old["source_original_q_animals"]
                    or saved["source_original_daylength_animals"]!=old["source_original_daylength_animals"]
                    or abs(saved["source_original_q"]-old["source_original_q"])>1e-12
                    or abs(saved["source_original_daylength_slope"]-old["source_original_daylength_slope"])>1e-12):
                raise ValueError("Original native GPS q and daylength stage changed")
        output["original_q_and_daylength_replayed_before_weather"]=True
        print("RSDS_NEW_CANONICAL_SAME_TAG_ORIGINAL_GPS_REPLAY",region,len(taxa),flush=True)
        cache={}
        for taxon in taxa:
            try:
                group=analyse_taxon_weather(taxon,tracks[taxon],previous[taxon],
                                             audit_rows[(region,taxon)],cache)
                output["taxa"].append(group)
            except Exception as exc:
                output["taxa"].append({"taxon":taxon,"status":"ORIGINAL_CONTEXT_OR_NUMERIC_HOLD",
                                      "error_type":type(exc).__name__,
                                      "source_animal_aggregate":None,
                                      "within_same_tag_rank2_minus_rank1":None})
            tracks[taxon].clear()
        output["unique_private_ERA5_gridyears_requested"]=len(cache)
        output["unique_canonical_weather_numeric_SHA_matching_frozen_audit"] = sum(
            v[1]=="CANONICAL_PHYSICAL_EXACT_SHA_MATCH" for v in cache.values())
        output["status"]="ORIGINAL_SOURCE_REANALYSIS_COMPLETED"
    except Exception as exc:
        output["status"]="ORIGINAL_SOURCE_HOLD"
        output["source_error_type"]=type(exc).__name__
        if not output["taxa"]:
            output["taxa"]=[{"taxon":t,"status":"ORIGINAL_SOURCE_HOLD"} for t in taxa]
    finally:
        private.unlink(missing_ok=True)
    out_dir.mkdir(parents=True,exist_ok=True)
    dest=out_dir/("source_"+region.lower()+"_new_physical_biology_v2.json")
    dest.write_text(json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                    encoding="utf-8")
    print("RSDS_CANONICAL_ERA5_NEW_ECOLOGY",json.dumps({
        "region":region,"status":output["status"],
        "original_source_q_replay":output["original_q_and_daylength_replayed_before_weather"],
        "supported_taxa":sum(r.get("status")=="CURRENT_VERSION_SOURCE_CONTEXTS_COMPLETE" for r in output["taxa"]),
        "canonical_sha_matches":output.get("unique_canonical_weather_numeric_SHA_matching_frozen_audit",0),
        "previous_12_ERA5_json_version_holds_unchanged":True
    },sort_keys=True),flush=True)
    return output


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--region",choices=("Utah","Wyoming"),required=True)
    cli.add_argument("--out",type=Path,required=True)
    a=cli.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    run_source(a.region,a.out)


if __name__=="__main__":
    main()
