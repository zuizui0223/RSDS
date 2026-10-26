"""Prospective six-original-study GPS time and heading ecology, first animal outcomes.

The six source study UUIDs/DOIs/taxa and thresholds were source-frozen before
opening new original GPS positions. Only source-only manifest-verified provider
bitstreams are accepted. The same native observation/geometry code is reused
for every source; no result-dependent replacements or new weather requests.

Only source/taxon and paired-tag aggregate statistics leave runner memory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import HERE, original_files, checksum_and_download
from scripts.global_solar_azimuth_12_new_original_events_v1 import (
    load_source_taxa, StrictSourceHold
)
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.solar_azimuth_macro_v0 import analyze_taxon as solar_taxon
from scripts.run_rsds_same_tag_gps_geometry_v2 import analyze_taxon as paired_geometry_taxon

CONTRACT = HERE / "research/RSDS_NEW_SIX_ORIGINAL_GPS_V0_CONTRACT.json"
MANIFEST = HERE / "research/RSDS_NEW_SIX_ORIGINAL_GPS_SOURCE_METADATA_V0_20261026.json"
SHA_RECEIPT = HERE / "research/RSDS_NEW_SIX_ORIGINAL_GPS_SOURCE_METADATA_BLOB_V0.json"
NEGATIVE = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def blob_sha(raw: bytes) -> str:
    return hashlib.sha1(("blob " + str(len(raw)) + "\0").encode("ascii") + raw).hexdigest()


def load_frozen():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    receipt = json.loads(SHA_RECEIPT.read_text(encoding="utf-8"))
    raw = MANIFEST.read_bytes()
    if (blob_sha(raw) != receipt["original_metadata_git_blob_sha"]
            or receipt["original_six_source_metadata_commit_phase"] != "BEFORE_NEW_GPS_EVENTS_OPENED"):
        raise ValueError("Frozen original study-only source manifest hash not verified")
    manifest = json.loads(raw.decode("utf-8"))
    studies = contract["original_fixed_sources"]
    entries = manifest["source_entries"]
    if (contract["version"] != "RSDS_NEW_SIX_ORIGINAL_PUBLIC_GPS_JOINT_TEMPORAL_DIRECTION_V0"
            or manifest["version"] != "RSDS_NEW_SIX_ORIGINAL_GPS_FILE_METADATA_V0"
            or contract["frozen_prior_source_shortwave_nonreplication"] != NEGATIVE
            or len(studies) != 6 or len(entries) != 6
            or manifest["original_animal_GPS_rows_read_by_preflight"] != 0
            or not manifest["no_source_replacement_or_effect_selection"]):
        raise ValueError("Frozen six study source frame not intact")
    for p, e in zip(studies, entries):
        if (e["package_uuid"] != p["id"] or e["original_dois"] != [p["doi"]]
                or e["original_species_labels"] != [p["taxon"]]
                or e["source_study_ids"] != [p["study_id"]]
                or e["registry_reported_animals"] != p["reported_animals"]
                or e["registry_reported_locations"] != p["reported_fixes"]):
            raise ValueError("Source preflight metadata changed after selection")
    return studies, entries


def source_analysis(source, preflight, working):
    outcome = {
        "version": "RSDS_NEW_SIX_ORIGINAL_GPS_FIRST_STUDY_RESULT_V0",
        "source_doI": source["doi"],
        "original_study_id": source["study_id"],
        "study_region_hint": source["region"],
        "guild_hint": source["guild"],
        "native_taxon": source["taxon"],
        "source_status": "ORIGINAL_SOURCE_HOLD",
        "source_metadata_status": preflight["source_metadata_status"],
        "all_outcome_axes": {
            "day_night_q": None,
            "daylength_slope": None,
            "sunward_raw": None,
            "sunward_joint": None,
            "same_tag_two_contexts": None,
        },
        "claim_ceiling": "Original GPS source study observational temporal displacement and solar-relative heading, not activity duration, solar compass, causal weather or independent biome inference.",
        "old_negative_external_shortwave_replication_UNCHANGED": NEGATIVE,
    }
    if preflight["source_metadata_status"] != "ORIGINAL_FILE_METADATA_ADMISSIBLE":
        outcome["source_status"] = "ORIGINAL_METADATA_HOLD"
        outcome["source_hold_reason"] = preflight.get("metadata_hold_reason", "ORIGINAL_METADATA_UNAVAILABLE")
        return outcome
    original = preflight["original_file"]
    ext = (".zip" if original["filename"].lower().endswith(".zip") else
           ".csv.gz" if original["filename"].lower().endswith(".csv.gz") else ".csv")
    local = working / ("private_original_" + source["id"] + ext)
    try:
        sources = [x for x in original_files(source["id"]) if x["bitstream_uuid"] == original["bitstream_uuid"]]
        if (len(sources) != 1 or sources[0]["bytes"] != original["bytes"]
                or sources[0]["checksum"] != original["checksum"]):
            raise StrictSourceHold("Original provider source MD5 / bytes / UUID no longer matches frozen metadata")
        outcome["original_provider_checksum_verification"] = checksum_and_download(sources[0], local)
        if (outcome["original_provider_checksum_verification"]["bytes"] >
                200_000_000):
            raise StrictSourceHold("Source exceeds predeclared original file budget")
        native, coverage, total = load_source_taxa(local, original["filename"], [source["taxon"]])
        outcome["source_native_observation_coverage"] = coverage[source["taxon"]]
        outcome["source_native_original_rows"] = total
        original_animals = native[source["taxon"]]
        if not original_animals:
            outcome["source_status"] = "NO_SOURCE_NATIVE_GPS_FIXES"
            return outcome
        # Reuse the source-blind observation rules from RSDS before any result
        # signs or intervals from this independent source are opened.
        q, _ = one_taxon(source["taxon"], original_animals)
        outcome["all_outcome_axes"]["day_night_q"] = q["day_night_speed_contrast"]
        outcome["all_outcome_axes"]["daylength_slope"] = q["daylength_association"]
        outcome["original_GPS_q_source_support"] = {
            "original_animals": q["total_original_animals"],
            "q_eligible_tagged_animals": q["admitted_animals_diel"],
            "daylength_eligible_tagged_animals": q["admitted_animals_photoperiod_slope"],
            "original_adjacent_hour_steps": q["original_adjacent_hour_steps"],
            "qualified_daypart_step_observations": q["qualified_daypart_step_observations"],
        }
        solar = solar_taxon({animal: list(by_hour.values()) for animal,by_hour in original_animals.items()})
        outcome["all_outcome_axes"]["sunward_raw"] = solar["results"]["raw_toward_sun_mean_cos"]
        outcome["all_outcome_axes"]["sunward_joint"] = solar["results"]["compass_independence_excess_cos"]
        outcome["original_solar_heading_support"] = {
            "source_animals": solar["source_animals"],
            "supported_animals": solar["supported_individuals"],
            "qualified_solar_daylight_moving_steps": solar["daylight_directional_qualified_steps"],
            "original_heading_status": solar["status"],
        }
        paired = paired_geometry_taxon(source["taxon"], original_animals)
        # Only public source-level aggregates; individual tags/grid keys private.
        outcome["all_outcome_axes"]["same_tag_two_contexts"] = paired
        outcome["source_status"] = "ORIGINAL_SOURCE_GPS_ECOLOGY_EVALUATED"
    except Exception as exc:
        # Earlier supported axes are kept with source failure, without inventing
        # missing downstream biological responses or adjusting thresholds.
        outcome["source_status"] = "ORIGINAL_GPS_SOURCE_HOLD"
        outcome["source_hold_error_type"] = type(exc).__name__
    finally:
        local.unlink(missing_ok=True)
    return outcome


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--index", type=int, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    fixed, source_preflight = load_frozen()
    if not 0 <= args.index < 6:
        p.error("source index must belong to six frozen studies [0,5]")
    args.out.mkdir(parents=True, exist_ok=True)
    result = source_analysis(fixed[args.index], source_preflight[args.index], args.out)
    target = args.out / ("original_independent_study_" + str(args.index) + ".json")
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print("RSDS_NEW_INDEPENDENT_SOURCE_ORIGINAL_FIRST", json.dumps({
        "source_index": args.index, "doi": result["source_doI"],
        "taxon": result["native_taxon"], "source_status": result["source_status"],
        "q_tag_support": result.get("original_GPS_q_source_support",{}).get("q_eligible_tagged_animals"),
        "sun_tag_support": result.get("original_solar_heading_support",{}).get("supported_animals"),
        "paired_tag_support": ((result["all_outcome_axes"]["same_tag_two_contexts"] or {})
                              .get("paired_context_effect") or {}).get("same_original_tagged_animals_all_axes")
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
