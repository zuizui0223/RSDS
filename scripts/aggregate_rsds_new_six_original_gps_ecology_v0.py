"""Source-complete original-GPS time and Sun-heading external ecological atlas.

Exactly six new source-only-frozen Movebank field studies. Every source and
undefined / broken original file is retained. Aggregated source-effect signs
are descriptive, not a world-wide probability sample or a climate mechanism.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

# Aggregation intentionally requires NO numpy or GPS/astronomy code. It reads
# only the six group-only JSON artifacts and source-only frozen metadata.
ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "research/RSDS_NEW_SIX_ORIGINAL_GPS_V0_CONTRACT.json"
MANIFEST = ROOT / "research/RSDS_NEW_SIX_ORIGINAL_GPS_SOURCE_METADATA_V0_20261026.json"
SHA_RECEIPT = ROOT / "research/RSDS_NEW_SIX_ORIGINAL_GPS_SOURCE_METADATA_BLOB_V0.json"
NEGATIVE = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
OUTPUT = ROOT / "research/RSDS_NEW_SIX_ORIGINAL_GPS_ECOLOGY_FIRST_RESULT_20261026.json"


def load_frozen():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    receipt = json.loads(SHA_RECEIPT.read_text(encoding="utf-8"))
    raw = MANIFEST.read_bytes()
    sha = hashlib.sha1(("blob " + str(len(raw)) + "\0").encode("ascii") + raw).hexdigest()
    if (sha != receipt.get("original_metadata_git_blob_sha")
            or receipt.get("original_six_source_metadata_commit_phase") != "BEFORE_NEW_GPS_EVENTS_OPENED"):
        raise ValueError("Source-only original GPS file snapshot SHA is not frozen")
    manifest = json.loads(raw.decode("utf-8"))
    sources = contract["original_fixed_sources"]
    entries = manifest["source_entries"]
    if (contract["version"] != "RSDS_NEW_SIX_ORIGINAL_PUBLIC_GPS_JOINT_TEMPORAL_DIRECTION_V0"
            or contract["frozen_prior_source_shortwave_nonreplication"] != NEGATIVE
            or manifest["version"] != "RSDS_NEW_SIX_ORIGINAL_GPS_FILE_METADATA_V0"
            or len(sources) != 6 or len(entries) != 6
            or manifest["original_animal_GPS_rows_read_by_preflight"] != 0):
        raise ValueError("Six original source studies altered or the pre-event frame was lost")
    for s, e in zip(sources, entries):
        if (s["id"] != e["package_uuid"] or
                [s["doi"]] != e["original_dois"] or
                [s["taxon"]] != e["original_species_labels"] or
                [s["study_id"]] != e["source_study_ids"]):
            raise ValueError("Original study source metadata identity is not frozen")
    return sources, entries


def supported(value, point_key, interval_key):
    if value is None:
        return {"status":"UNAVAILABLE","n_animals":0,"mean":None,"bootstrap95":None}
    n = value["n_individuals"]
    point = value.get(point_key)
    ci = value.get(interval_key)
    if not isinstance(n,int) or n<0:
        raise ValueError("Invalid original tagged animal source count")
    if n < 5:
        if point is not None or ci is not None:
            raise ValueError("Source endpoint reported without five original tagged animals")
        return {"status":"UNAVAILABLE","n_animals":n,"mean":None,"bootstrap95":None}
    if point is None or ci is None or len(ci)!=2:
        raise ValueError("Admissible source effect absent or bootstrap unavailable")
    if not all(isinstance(v,(int,float)) and math.isfinite(v) for v in [point]+list(ci)):
        raise ValueError("Nonfinite original source effect")
    if not ci[0] <= point <= ci[1]:
        raise ValueError("Bootstrap interval inconsistent with source group outcome")
    sign = "POSITIVE_NOMINAL_BOOTSTRAP" if ci[0]>0 else \
           "NEGATIVE_NOMINAL_BOOTSTRAP" if ci[1]<0 else "INCLUDES_ZERO"
    return {"status":sign,"n_animals":n,"mean":point,"bootstrap95":ci}


def summarize(results):
    studies, preflight = load_frozen()
    if len(results)!=6 or {str(i) for i in range(6)}!=set(results):
        raise ValueError("Six fixed independent original GPS study results missing")
    rows = []
    count = {"q":Counter(),"daylength":Counter(),"sun_raw":Counter(),"sun_joint":Counter(),
             "paired":Counter(),"source":Counter()}
    for index,(s,pr) in enumerate(zip(studies,preflight)):
        item = results[str(index)]
        if (item.get("native_taxon") != s["taxon"]
            or item.get("original_study_id") != s["study_id"]
            or item.get("source_doI") != s["doi"]
            or item.get("source_metadata_status") != pr["source_metadata_status"]
            or item.get("old_negative_external_shortwave_replication_UNCHANGED") != NEGATIVE):
            raise ValueError("Original fixed source/taxon/outcome status was changed")
        status=item["source_status"]
        count["source"][status]+=1
        response=item["all_outcome_axes"]
        if set(response)!={"day_night_q","daylength_slope","sunward_raw",
                            "sunward_joint","same_tag_two_contexts"}:
            raise ValueError("Incomplete frozen original source endpoint response axes")
        q=supported(response["day_night_q"],"mean","animal_bootstrap_95_percentile")
        light=supported(response["daylength_slope"],"mean","animal_bootstrap_95_percentile")
        raw=supported(response["sunward_raw"],"mean","bootstrap95")
        joint=supported(response["sunward_joint"],"mean","bootstrap95")
        paired=response["same_tag_two_contexts"]
        paired_state=paired["source_status"] if paired else "UNAVAILABLE"
        if paired is not None:
            if paired_state not in ("PAIRED_TAG_ESTIMABLE","PAIRED_TAG_UNAVAILABLE"):
                raise ValueError("Unrecognized same-tag context comparison")
            p=paired["paired_context_effect"]
            if paired_state=="PAIRED_TAG_ESTIMABLE" and \
               (p is None or p["same_original_tagged_animals_all_axes"]<5):
                raise ValueError("Source claims unsupported paired original tagged animals")
        else:
            p=None
        for k,v in (("q",q),("daylength",light),("sun_raw",raw),("sun_joint",joint)):
            count[k][v["status"]]+=1
        count["paired"][paired_state]+=1
        rows.append({
            "fixed_original_source_rank":index+1,
            "study_doi":s["doi"],"native_taxon":s["taxon"],
            "guild_metadata_hint":s["guild"],"region_metadata_hint":s["region"],
            "original_metadata_state":pr["source_metadata_status"],
            "original_GPS_outcome_status":status,
            "source_q_night_minus_day":q,
            "source_q_daylength_per_hour":light,
            "source_solar_relative_heading_cos":raw,
            "source_joint_heading_sun_covariance":joint,
            "within_same_original_tag_two_q_coverage_contexts":{
                "status":paired_state,
                "actual_paired_tagged_animal_support":(
                    p["same_original_tagged_animals_all_axes"] if p else None),
                "source_coverage_candidate_pairs":(
                    paired.get("eligible_context_pairs") if paired else None),
                "source_original_eligible_q_gridyears":(
                    paired.get("eligible_q_coverage_gridyears") if paired else None),
                "source_paired_context_time_and_heading_effects":p
            },
            "source_missingness_or_hold_reason":item.get("source_hold_reason")
                or item.get("source_hold_error_type")
        })
    return {
        "version":"RSDS_NEW_SIX_ORIGINAL_PUBLIC_GPS_TIME_DIRECTION_V0_FIRST_RESULT",
        "date":"2026-10-26",
        "research_stage":"PRESELECTED_SIX_ORIGINAL_STUDIES_POST_EVENT_DESCRIPTIVE",
        "frozen_source_original_study_count":6,
        "unit_of_source_replication":"source original DOI/study ID; six purposively selected studies are NOT representative of world ecology",
        "six_full_source_study_taxon_results":rows,
        "number_of_studies_by_source_status":dict(count["source"]),
        "source_sample_sign_support_counts":{
            axis:dict(count[axis]) for axis in ("q","daylength","sun_raw","sun_joint","paired")
        },
        "original_metadata_preflight":str(MANIFEST.relative_to(ROOT)),
        "previous_external_48_study_solar_radiation_nonreplication_UNCHANGED":NEGATIVE,
        "previous_12_hourly_ERA5_JSON_SHA_HOLDS_UNCHANGED":True,
        "claim_ceiling":[
            "Source sampling was fixed before new GPS events, but was purposively selected on public metadata and study title, not random ecosystems or a prospectively preregistered response hypothesis.",
            "Six original source datasets can contain non-independent regions and life-history states, and may have different GPS duty cycles and taxon-level field capture biases.",
            "Original q is NIGHT−DAY GPS minimum straight-line displacement speed, not actual active time allocation. It has different sampling from solar-bearing days and/or individuals in the same source.",
            "Raw cos(GPS travel bearing−solar azimuth) can arise from geographical route marginals, not biological Sun compass navigation.",
            "Source q and raw heading averages may be estimated from different sets of tagged animals. Their group-level comparison is not within-animal covariance.",
            "Two-context paired results are selected by ORIGINAL eligible q-tag support, not by sunlight/radiation/temperature/season levels. Some regions may have no identifiable paired contexts.",
            "Distinct unsigned or zero-crossing source bootstraps represent conditional sampling uncertainty, not equivalent species effects or global ecological null results.",
            "All unsupported original source MD5, file format, row coverage, and weather version holds remain distinct from zero biological response; no new weather was requested.",
            "No original tagged animal IDs, exact event times/coordinates, raw GPS files or individual source slopes were published."
        ]
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--out",type=Path,default=OUTPUT)
    args=parser.parse_args()
    full={}
    for i in range(6):
        file=args.input/("original_independent_study_"+str(i)+".json")
        if not file.is_file():
            raise ValueError("Frozen study result missing; no source replacement: "+str(i))
        full[str(i)]=json.loads(file.read_text(encoding="utf-8"))
    result=summarize(full)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_EXTERNAL_SIX_ORIGINAL_GPS_RESULT",json.dumps({
        "study_total":result["frozen_source_original_study_count"],
        "source_statuses":result["number_of_studies_by_source_status"],
        "axis_support":result["source_sample_sign_support_counts"],
        "older_external_radiation":result["previous_external_48_study_solar_radiation_nonreplication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
