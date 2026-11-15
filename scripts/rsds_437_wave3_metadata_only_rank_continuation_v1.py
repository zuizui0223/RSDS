"""Frozen metadata-only Wave3 continuation of the 2026-10-19 Movebank source ranking.

NO new provider file listings, source telemetry, animal/tag IDs, coordinates,
timestamps, source q, or daylength outcomes are used for source selection.
Wave3 quota was chosen post-Wave2 (not a prospective biological replication).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import (
    HERE,
    git_blob,
    pinned,
    run as fixed_original_frame,
)

CONTRACT = HERE / "research/RSDS_437_WAVE3_METADATA_ONLY_SOURCE_RANKS25_40_CONTRACT_V1.json"
VERSION = "RSDS_437_WAVE3_METADATA_ONLY_RANK_CONTINUATION_V1"
OLD_NEG = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def _verify_scientific_blob_without_reading_outcome(path: str, sha: str) -> None:
    """Pin historical results without loading their scientific coefficients."""
    raw = (HERE / path).read_bytes()
    if git_blob(raw) != sha:
        raise ValueError("Historical science file changed: " + path)


def freeze() -> dict:
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    required = {
        "version": VERSION,
        "date": "2026-11-15",
        "fixed_catalog_git_blob": "9822ec80ed98cab1e3d17a853d0ce5e63e47298c",
        "frozen_frame_script_git_blob": "1fbb6e4f7ef686aed63dfd4433ee425e751592ca",
        "original_public_packages": 437,
        "original_GPS_GNSS_packages": 259,
        "eligible_metadata_sources": 120,
        "candidate_single_native": 109,
        "candidate_multi_native": 11,
        "previously_selected_single_native": 24,
        "previously_selected_multi_native": 11,
        "wave3_single_native_start_rank_inclusive": 25,
        "wave3_single_native_end_rank_inclusive": 40,
        "wave3_multi_native_added": 0,
        "wave3_new_original_source_studies": 16,
        "after_wave3_selected_original_source_studies": 51,
        "after_wave3_unselected_single_native": 69,
        "original_meta_rank_salt": "RSDS_META_CENSUS_PRE_EVENT_WAVE1_20261115",
        "no_GPS_or_provider_file_read_for_this_route": True,
        "new_original_GPS_events_opened": 0,
        "new_provider_file_metadata_opened": False,
        "old_48_source_signed_solar_radiation": OLD_NEG,
        "old_12_ERA5_raw_SHA_HOLD": True,
    }
    for key, expected in required.items():
        if c.get(key) != expected:
            raise ValueError(f"Wave3 scientific freeze/denominator changed: {key}")
    if git_blob((HERE / c["frozen_frame_script"]).read_bytes()) != c["frozen_frame_script_git_blob"]:
        raise ValueError("Original metadata-only ranking algorithm changed")
    _verify_scientific_blob_without_reading_outcome(
        c["frozen_35_source_q_path"], c["frozen_35_source_q_git_blob"]
    )
    _verify_scientific_blob_without_reading_outcome(
        c["frozen_35_source_sign_sensitivity_path"],
        c["frozen_35_source_sign_sensitivity_git_blob"],
    )
    return c


def select() -> dict:
    c = freeze()
    wave1 = pinned(c["wave1_selection_path"], c["wave1_selection_git_blob"])
    wave2 = pinned(c["wave2_selection_path"], c["wave2_selection_git_blob"])
    frame = fixed_original_frame()

    summary = frame["summary"]
    if (summary["all_original_packages"] != 437
        or summary["gps_gnss_packages"] != 259
        or summary["metadata_candidate_packages"] != 120
        or summary["metadata_candidate_single_native"] != 109
        or summary["metadata_candidate_multi_native"] != 11):
        raise ValueError("Original source-only Movebank frame no longer matches")
    if (wave1["original_registry_package_count"] != 437
        or wave1["meta_stage1_candidate_count"] != 120
        or len(wave1["pre_event_selected_sources_in_rank_order"]) != 16
        or len(wave2["all_original_wave2_sources_in_rank_order"]) != 19
        or wave2["original_wave2_GPS_events_opened"] != 0
        or wave2["original_wave2_provider_files_opened"] != 0
        or wave2["all_previously_unselected_metadata_candidates_remaining"] != 85):
        raise ValueError("Historical source selection or privacy frame changed")

    ranked = {}
    for stratum, count in (("single_native", 109), ("multi_native", 11)):
        rows = [
            row for row in frame["all_437_metadata_status_rows"]
            if row["stage1_status"] == "FIRST_WAVE_METADATA_CANDIDATE"
            and row["selection_stratum"] == stratum
        ]
        ranked[stratum] = sorted(rows, key=lambda r: (r["lottery_sha256"], r["package_uuid"]))
        if len(ranked[stratum]) != count:
            raise ValueError("Original metadata candidate stratum changed")

    w1 = wave1["pre_event_selected_sources_in_rank_order"]
    w2 = wave2["all_original_wave2_sources_in_rank_order"]
    if [r["package_uuid"] for r in ranked["single_native"][:8]] != [
        r["package_uuid"] for r in w1[:8]
    ]:
        raise ValueError("Wave1 single-native source prefix changed")
    if [r["package_uuid"] for r in ranked["multi_native"][:8]] != [
        r["package_uuid"] for r in w1[8:]
    ]:
        raise ValueError("Wave1 multi-native source prefix changed")
    if [r["package_uuid"] for r in ranked["single_native"][8:24]] != [
        r["package_uuid"] for r in w2[:16]
    ]:
        raise ValueError("Wave2 single-native source continuation changed")
    if [r["package_uuid"] for r in ranked["multi_native"][8:11]] != [
        r["package_uuid"] for r in w2[16:]
    ]:
        raise ValueError("Wave2 multi-native source continuation changed")

    new_rows = ranked["single_native"][24:40]
    old_rows = ranked["single_native"][:24] + ranked["multi_native"]
    if len(new_rows) != 16 or len(old_rows) != 35:
        raise ValueError("Wave3 source quota changed")
    all_rows = old_rows + new_rows
    if len({r["package_uuid"] for r in all_rows}) != 51:
        raise ValueError("Wave3 duplicate package UUID")
    if len({r["repository_study_ids"][0] for r in all_rows}) != 51:
        raise ValueError("Wave3 duplicate original study ID")
    if len({r["doi_labels"][0] for r in all_rows}) != 51:
        raise ValueError("Wave3 duplicate DOI")
    if any(len(set(r["scientific_name_labels"])) != 1 for r in new_rows):
        raise ValueError("Wave3 is single-native only")

    selected = []
    for i, row in enumerate(new_rows):
        names = row["scientific_name_labels"]
        selected.append({
            "wave3_index": i,
            "source_rank_1_based": 25 + i,
            "original_repository_study_id": row["repository_study_ids"][0],
            "original_DOI_labels": row["doi_labels"],
            "package_uuid": row["package_uuid"],
            "scientific_name_labels": names,
            "original_sensor_labels": row["original_sensors"],
            "original_license_labels": row["original_license_labels"],
            "reported_animal_count_metadata_only": row["reported_animal_count_not_deduplicated"],
            "reported_location_count_metadata_only": row["reported_location_count_not_qualified"],
            "source_lottery_SHA256": row["lottery_sha256"],
            "human_privacy_file_and_GPS_HOLD": "Homo sapiens" in names,
            "new_provider_file_metadata_opened": False,
            "new_GPS_animal_events_opened": 0,
        })

    return {
        "version": VERSION,
        "stage": "SOURCE_ONLY_WAVE3_IDENTITIES_FROZEN_NO_NEW_PROVIDER_OR_GPS_ACCESS",
        "original_fixed_2026_10_19_catalog_git_blob": c["fixed_catalog_git_blob"],
        "old_wave1_selection_git_blob": c["wave1_selection_git_blob"],
        "old_wave2_selection_git_blob": c["wave2_selection_git_blob"],
        "old_wave1_plus_wave2_35_source_science_git_blob": c["frozen_35_source_q_git_blob"],
        "source_frame_packages": 437,
        "source_frame_GPS_packages": 259,
        "source_frame_metadata_priority_candidates": 120,
        "prior_selected_source_studies": 35,
        "wave3_new_original_source_studies": 16,
        "wave3_selected_single_native_rank_range_inclusive": [25, 40],
        "selected_original_source_studies_after_wave3": 51,
        "remaining_metadata_priority_single_native_sources": 69,
        "all_wave3_selected_source_metadata_in_SHA_rank_order": selected,
        "wave3_new_original_animal_GPS_events_opened": 0,
        "wave3_new_original_provider_file_metadata_opened": False,
        "previous_48_source_solar_direction_external_nonreplication_UNCHANGED": True,
        "previous_12_source_ERA5_raw_SHA_HOLD_UNCHANGED": True,
        "scientific_boundary": c["inference_boundary"],
        "quota_selected_after_two_prior_waves": True,
        "global_biological_probability_sample": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    result = select()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    print("RSDS_WAVE3_METADATA_ONLY_LOCKED", json.dumps({
        "selected": result["wave3_new_original_source_studies"],
        "ranks": result["wave3_selected_single_native_rank_range_inclusive"],
        "total_source_studies": result["selected_original_source_studies_after_wave3"],
        "remaining": result["remaining_metadata_priority_single_native_sources"],
        "provider_file_metadata_opened": False,
        "GPS_events_opened": 0,
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
