"""Wave3 ORIGINAL bundle source-file metadata ONLY, after frozen new source identities.

Checks public provider original-file list, UUID, byte size, and MD5 hint.
Does NOT stream, parse, download or open original animal GPS/time/position.
The selected 16 are fixed source identities; any failed source stays HOLD.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned
from scripts.rsds_437_wave3_metadata_only_rank_continuation_v1 import (
    select as frozen_metadata_only_select,
)
from scripts.rsds_437_wave1_original_file_metadata_v1 import original_primary
from scripts.movebank_native_ecology_20261019 import original_files

SELECTION_PATH = "research/RSDS_437_WAVE3_16_SOURCE_METADATA_ONLY_FIRST_RESULT_20261115.json"
SELECTION_GIT_BLOB = "a5a11c75a08df99755d3f6d560d111ab4904375b"
NEG = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
VERSION = "RSDS_437_WAVE3_16_ORIGINAL_PROVIDER_FILE_METADATA_FIRST_V1"


def freeze() -> list[dict]:
    doc = pinned(SELECTION_PATH, SELECTION_GIT_BLOB)
    replay = frozen_metadata_only_select()
    if (doc != replay
        or doc["version"] != "RSDS_437_WAVE3_METADATA_ONLY_RANK_CONTINUATION_V1"
        or doc["stage"] != "SOURCE_ONLY_WAVE3_IDENTITIES_FROZEN_NO_NEW_PROVIDER_OR_GPS_ACCESS"
        or len(doc["all_wave3_selected_source_metadata_in_SHA_rank_order"]) != 16
        or doc["prior_selected_source_studies"] != 35
        or doc["selected_original_source_studies_after_wave3"] != 51
        or doc["wave3_selected_single_native_rank_range_inclusive"] != [25, 40]
        or doc["wave3_new_original_animal_GPS_events_opened"] != 0
        or doc["wave3_new_original_provider_file_metadata_opened"] is not False
        or doc["previous_48_source_solar_direction_external_nonreplication_UNCHANGED"] is not True
        or doc["previous_12_source_ERA5_raw_SHA_HOLD_UNCHANGED"] is not True):
        raise ValueError("Original Wave3 metadata-source selection manifest changed")
    rows = doc["all_wave3_selected_source_metadata_in_SHA_rank_order"]
    if ([r["source_rank_1_based"] for r in rows] != list(range(25, 41))
        or [r["wave3_index"] for r in rows] != list(range(16))
        or len({r["original_repository_study_id"] for r in rows}) != 16):
        raise ValueError("Original frozen Wave3 rank/study identity changed")
    return rows


def inspect(source: dict) -> dict:
    report = {
        "wave3_index": source["wave3_index"],
        "original_single_native_source_rank_1_based": source["source_rank_1_based"],
        "package_uuid": source["package_uuid"],
        "original_repository_study_id": source["original_repository_study_id"],
        "original_DOI_labels": source["original_DOI_labels"],
        "registered_native_taxa": source["scientific_name_labels"],
        "selection_stratum": "single_native",
        "original_file": None,
        "source_status": "SOURCE_METADATA_HOLD",
        "provider_original_file_list_obtained": False,
        "animal_GPS_events_opened": 0,
    }
    if "Homo sapiens" in source["scientific_name_labels"]:
        report["source_status"] = "HUMAN_LOCATION_PRIVACY_HOLD"
        report["source_hold_reason"] = "NO_HUMAN_TRAJECTORY_SOURCE_ACCESS"
        return report
    try:
        original_file_inventory = original_files(source["package_uuid"])
        report["provider_original_file_list_obtained"] = True
        original, hold, n, rule = original_primary(original_file_inventory)
        report["candidate_file_count_after_filename_rules"] = n
        report["file_selection_rule"] = rule
        if original is None:
            report["source_hold_reason"] = hold
        else:
            report["source_status"] = "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
            report["original_file"] = original
    except Exception as e:
        report["source_hold_reason"] = "PROVIDER_FILE_METADATA_ERROR_" + type(e).__name__
    return report


def evaluate() -> dict:
    fixed = freeze()
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(inspect, r): r["wave3_index"] for r in fixed}
        by_index = {futures[f]: f.result() for f in as_completed(futures)}
    rows = [by_index[i] for i in range(16)]
    if ([r["original_repository_study_id"] for r in rows] !=
        [r["original_repository_study_id"] for r in fixed]):
        raise ValueError("Provider file preflight substituted original source")
    if ([r["package_uuid"] for r in rows] != [r["package_uuid"] for r in fixed]
        or any(r["animal_GPS_events_opened"] != 0 for r in rows)
        or any(r["source_status"] == "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
               and r["original_file"] is None for r in rows)):
        raise ValueError("Provider metadata source identity/file/GPS integrity error")
    if any(("Homo sapiens" in old["scientific_name_labels"]) !=
           (new["source_status"] == "HUMAN_LOCATION_PRIVACY_HOLD")
           for old, new in zip(fixed, rows)):
        raise ValueError("Human-location privacy violated")
    return {
        "version": VERSION,
        "stage": "AFTER_SHA_FROZEN_WAVE3_SOURCE_SELECTION_BEFORE_ORIGINAL_GPS_EVENTS",
        "original_wave3_selection_git_blob": SELECTION_GIT_BLOB,
        "original_wave3_selection_source_study_slots": 16,
        "source_study_slots": 16,
        "new_wave3_single_native": 16,
        "new_wave3_multi_native": 0,
        "old_selected_source_studies": 35,
        "total_selected_original_source_studies_after_wave3": 51,
        "source_status_counts": dict(sorted(Counter(r["source_status"] for r in rows).items())),
        "source_rows": rows,
        "provider_file_metadata_only": True,
        "human_GPS_or_position_data_retrieved": False,
        "all_new_wave3_animal_GPS_events_opened": 0,
        "original_48_source_radiation_negative_external_replication_UNCHANGED": NEG,
        "old_12_ERA5_original_raw_SHA_HOLD_UNCHANGED": True,
        "scientific_claim_ceiling": (
            "Original provider file UUID/bytes/MD5/filename only. Source events, GPS schema, "
            "animal biology, response completeness, and population ecology have not been "
            "observed. Every source and technical/privacy HOLD retains its frozen rank."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    a = parser.parse_args()
    result = evaluate()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("RSDS_WAVE3_PROVIDER_ORIGINAL_METADATA", json.dumps({
        "source_study_slots": result["source_study_slots"],
        "source_status_counts": result["source_status_counts"],
        "sources": [
            {"i": r["wave3_index"], "doi": r["original_DOI_labels"][0],
             "status": r["source_status"], "reason": r.get("source_hold_reason"),
             "original_bytes": (r["original_file"] or {}).get("bytes")}
            for r in result["source_rows"]
        ],
        "original_GPS_event_records": 0,
        "human_GPS_or_position_retrieved": False,
    }, ensure_ascii=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
