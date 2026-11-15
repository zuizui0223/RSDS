"""First original GPS q/daylength on SHA-frozen Wave3 source IDs.

No source identities are chosen here: original 16 source IDs and 11 file MD5
candidates were committed before any Wave3 animal GPS outcome read.
Source-wide original GPS errors fail closed; HOLD and null are not q=0.
Output source-native aggregate results only, never original tracks or tags.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.rsds_437_systematic_source_frame_v1 import pinned
from scripts.rsds_437_wave3_original_file_metadata_v1 import freeze as frozen_wave3_source_selection

SPEC = HERE / "research/RSDS_437_WAVE3_ORIGINAL_GPS_FIRST_EFFECTS_V1_CONTRACT.json"
NEG = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
VERSION = "RSDS_437_WAVE3_16_ORIGINAL_GPS_FIRST_EFFECTS_CONTRACT_V1"
FIRST_RESULT_VERSION = "RSDS_WAVE3_16_ORIGINAL_GPS_FIRST_SOURCE_Q_DAYLENGTH_V1"


def freeze_original():
    c = json.loads(SPEC.read_text(encoding="utf-8"))
    required = {
        "version": VERSION,
        "date": "2026-11-15",
        "wave3_source_slots": 16,
        "wave3_file_metadata_admissible": 11,
        "wave3_source_file_metadata_HOLD": 5,
        "original_file_max_bytes": 200_000_000,
        "original_read_max_CSV_rows": 2_000_000,
        "new_wave3_GPS_event_rows_opened_prior_to_contract": 0,
        "previous_48_signed_irradiance_external": NEG,
        "previous_12_ERA5_original_raw_SHA_HOLD": True,
        "source_only_selection_git_blob": "a5a11c75a08df99755d3f6d560d111ab4904375b",
        "source_provider_file_metadata_git_blob": "02253b1743500bc51fc166b249b7d3660ed9610a",
        "original_previous_35_source_result_git_blob": "9b02774f7c31d29b2eec5039a80431a895dccb36",
    }
    for k, expected in required.items():
        if c.get(k) != expected:
            raise ValueError("Source Wave3 pre-outcome contract altered: "+k)
    old = pinned(
        c["original_previous_35_source_result_path"],
        c["original_previous_35_source_result_git_blob"],
    )
    if (old["original_source_study_units"] != 35
        or len(old["combined_source_study_summary_rows"]) != 35
        or old["previous_48_source_sunlight_signed_external_not_replicated"] is not True
        or old["previous_12_ERA5_raw_SHA_HOLD"] is not True):
        raise ValueError("Original Wave1+Wave2 science changed before Wave3")
    selection = frozen_wave3_source_selection()
    src = pinned(
        c["source_provider_file_metadata_path"],
        c["source_provider_file_metadata_git_blob"],
    )
    sources = src["source_rows"]
    if (len(selection) != 16 or len(sources) != 16
        or src["source_status_counts"] != {
            "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE": 11,
            "SOURCE_METADATA_HOLD": 5,
        }
        or src["all_new_wave3_animal_GPS_events_opened"] != 0
        or src["human_GPS_or_position_data_retrieved"] is not False
        or src["original_wave3_selection_git_blob"] != c["source_only_selection_git_blob"]):
        raise ValueError("Original Wave3 file/gate receipt was modified")
    for i, (a, b) in enumerate(zip(selection, sources)):
        if (a["wave3_index"] != b["wave3_index"] or
            a["wave3_index"] != i or
            a["source_rank_1_based"] != b["original_single_native_source_rank_1_based"] or
            a["package_uuid"] != b["package_uuid"] or
            a["original_repository_study_id"] != b["original_repository_study_id"] or
            a["original_DOI_labels"] != b["original_DOI_labels"] or
            a["scientific_name_labels"] != b["registered_native_taxa"] or
            len(b["registered_native_taxa"]) != 1):
            raise ValueError("Wave3 source/file/taxon changed after SHA pre-event freeze")
    prior_ids = {s["source_original_repository_study_id"] for s in
                 old["combined_source_study_summary_rows"]}
    if len(prior_ids) != 35 or any(
        r["original_repository_study_id"] in prior_ids for r in sources
    ):
        raise ValueError("New source study ID aliases previously seen biological study")
    return c, sources


def process(index: int) -> dict:
    _, sources = freeze_original()
    if type(index) is not int or not 0 <= index < 16:
        raise ValueError("Only sixteen fixed Wave3 original source slots")
    source = sources[index]
    taxa = source["registered_native_taxa"]
    result = {
        "version": FIRST_RESULT_VERSION,
        "wave3_index": index,
        "wave3_original_source_rank_1_based": source["original_single_native_source_rank_1_based"],
        "source_package_uuid": source["package_uuid"],
        "source_original_repository_study_id": source["original_repository_study_id"],
        "source_original_DOI_labels": source["original_DOI_labels"],
        "source_metadata_stratum": source["selection_stratum"],
        "original_registered_native_taxa": taxa,
        "original_file_source_status": source["source_status"],
        "original_file_MD5_bytes_verified": False,
        "source_science_status": "SOURCE_GPS_HOLD",
        "source_hold_stage": None,
        "source_hold_reason": None,
        "native_group_results": [
            {"taxon": sp, "q": None, "daylength": None, "status": "SOURCE_UNEVALUATED"}
            for sp in taxa
        ],
        "native_joint_calendar_covariance": {
            "status": "NOT_APPLICABLE_SINGLE_NATIVE_SOURCE"
        },
        "original_source_GPS_event_rows_opened": 0,
        "human_mobility_GPS_data_opened": False,
        "old_48_signed_solar_external_nonreplication_UNCHANGED": NEG,
        "old_12_ERA5_source_raw_SHA_HOLD_UNCHANGED": True,
        "original_raw_GPS_dates_tag_IDs_coordinates_exported": False,
    }
    if "Homo sapiens" in taxa:
        if source["source_status"] != "HUMAN_LOCATION_PRIVACY_HOLD":
            raise ValueError("Original human source no-access gate reclassified")
        result["source_science_status"] = "HUMAN_LOCATION_PRIVACY_HOLD"
        result["source_hold_stage"] = "HUMAN_SOURCE_PRIVACY_NO_ACCESS"
        result["source_hold_reason"] = "NO_HUMAN_TRAJECTORY_SOURCE_ACCESS"
        return result
    if source["source_status"] != "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":
        result["source_science_status"] = "SOURCE_FILE_METADATA_HOLD"
        result["source_hold_stage"] = "SOURCE_ONLY_ORIGINAL_FILE_METADATA"
        result["source_hold_reason"] = source.get("source_hold_reason")
        return result

    meta = source["original_file"]
    extension = (
        ".zip" if meta["filename"].lower().endswith(".zip") else
        ".csv.gz" if meta["filename"].lower().endswith(".csv.gz") else ".csv"
    )
    stage = "ORIGINAL_UUID_BYTES_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_wave3_private_gps_") as directory:
            filepath = Path(directory) / ("original" + extension)
            current = [
                f for f in original_files(source["package_uuid"])
                if f["bitstream_uuid"] == meta["bitstream_uuid"]
            ]
            if (len(current) != 1
                or current[0]["bytes"] != meta["bytes"]
                or current[0].get("checksum", {}).get("checkSumAlgorithm", "").upper() != "MD5"
                or current[0]["checksum"]["value"].lower() != meta["checksum"]["value"]):
                raise ValueError("Original source bitstream file changed after metadata freeze")
            confirmed = checksum_and_download(current[0], filepath)
            if confirmed.get("provider_checksum_verified") is not True:
                raise ValueError("Original provider file MD5 not verified")
            result["original_file_MD5_bytes_verified"] = True

            stage = "ORIGINAL_NATIVE_GPS_Q_AND_DAYLENGTH"
            tracks, coverage, report = load_source_taxa(filepath, meta["filename"], taxa)
            result["original_source_GPS_event_rows_opened"] = report["original_total_rows"]
            native = []
            for sp in taxa:
                data, _ = one_taxon(sp, tracks[sp])
                native.append({
                    "taxon": sp,
                    "original_animal_count": data["total_original_animals"],
                    "original_hour_selected_fixes": coverage[sp]["hour_selected_fixes"],
                    "original_adjacent_intervals": data["original_adjacent_hour_steps"],
                    "qualified_same_daypart_adjacent_intervals": data["qualified_daypart_step_observations"],
                    "n_q_tagged_animals": data["admitted_animals_diel"],
                    "n_daylength_tagged_animals": data["admitted_animals_photoperiod_slope"],
                    "q": data["day_night_speed_contrast"],
                    "daylength": data["daylength_association"],
                    "status": (
                        "NATIVE_Q_ESTIMABLE"
                        if data["day_night_speed_contrast"] is not None
                        else "NATIVE_Q_UNAVAILABLE"
                    ),
                })
            result["native_group_results"] = native
            result["source_science_status"] = (
                "SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE"
                if any(x["q"] is not None for x in native)
                else "SOURCE_NO_NATIVE_Q_ESTIMABLE"
            )
            tracks.clear()
    except Exception as exc:
        result["source_science_status"] = "SOURCE_GPS_HOLD"
        result["source_hold_stage"] = stage
        result["source_hold_reason"] = type(exc).__name__
        result["native_group_results"] = [
            {"taxon": sp, "q": None, "daylength": None, "status": "SOURCE_WIDE_FAIL_CLOSED"}
            for sp in taxa
        ]
        result["native_joint_calendar_covariance"] = {"status": "SOURCE_WIDE_HOLD"}
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    result = process(args.index)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print("RSDS_WAVE3_GPS_ORIGINAL_SOURCE", json.dumps({
        "index": args.index,
        "original_source_rank": result["wave3_original_source_rank_1_based"],
        "source_science_status": result["source_science_status"],
        "native_q_supported": sum(x["q"] is not None for x in result["native_group_results"]),
        "source_MD5_verified": result["original_file_MD5_bytes_verified"],
        "original_event_rows_opened": result["original_source_GPS_event_rows_opened"],
        "hold_stage": result["source_hold_stage"],
        "hold_type": result["source_hold_reason"],
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
