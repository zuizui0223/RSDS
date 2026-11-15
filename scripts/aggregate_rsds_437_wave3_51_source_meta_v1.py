"""Frozen study-ID first aggregation of Wave3 original GPS sources (35 + 16).

One original Movebank repository study ID is a replication unit; not a
tag, event, native species or GPS fix. All 16 Wave3 source slots are required.
Preserves old 35 source observations and the historical solar nonreplication.
This route gives SOURCE-RESPONDENT descriptions, not global ecology estimates.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from statistics import mean

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned
from scripts.run_rsds_437_wave3_original_gps_first_effects_v1 import (
    freeze_original,
    FIRST_RESULT_VERSION,
    NEG,
)

def native_mean_if_supported(effect):
    if effect is None:
        return None
    if not isinstance(effect, dict):
        raise ValueError("Invalid original native effect object")
    n = effect.get("n_individuals")
    mu = effect.get("mean")
    ci = effect.get("animal_bootstrap_95_percentile")
    if (not isinstance(n, int) or n < 5
        or not isinstance(mu, (int, float)) or not math.isfinite(mu)
        or not isinstance(ci, list) or len(ci) != 2
        or not all(isinstance(x, (float, int)) and math.isfinite(x) for x in ci)
        or ci[0] > ci[1]):
        raise ValueError("Invalid original source-native sample/effect/bootstrap")
    return float(mu)


def build(first_results: list[dict]) -> dict:
    if len(first_results) != 16:
        raise ValueError("Exactly sixteen original Wave3 sources including HOLD required")
    c, source_roster = freeze_original()
    old = pinned(
        c["original_previous_35_source_result_path"],
        c["original_previous_35_source_result_git_blob"],
    )
    original = old["combined_source_study_summary_rows"]
    if len(original) != 35 or old["original_source_study_units"] != 35:
        raise ValueError("Historical original 35 studies missing")
    seen = {r["source_original_repository_study_id"] for r in original}
    if len(seen) != 35:
        raise ValueError("Original 35 source-ID duplicates")
    new = []
    for i, (s, row) in enumerate(zip(source_roster, first_results)):
        if (
            row.get("version") != FIRST_RESULT_VERSION or
            row.get("wave3_index") != i or
            row.get("wave3_original_source_rank_1_based") != i+25 or
            row.get("source_package_uuid") != s["package_uuid"] or
            row.get("source_original_repository_study_id") != s["original_repository_study_id"] or
            row.get("source_original_DOI_labels") != s["original_DOI_labels"] or
            row.get("source_metadata_stratum") != "single_native" or
            row.get("original_registered_native_taxa") != s["registered_native_taxa"] or
            row.get("original_file_source_status") != s["source_status"] or
            row.get("old_48_signed_solar_external_nonreplication_UNCHANGED") != NEG or
            row.get("old_12_ERA5_source_raw_SHA_HOLD_UNCHANGED") is not True or
            row.get("human_mobility_GPS_data_opened") is not False or
            row.get("original_raw_GPS_dates_tag_IDs_coordinates_exported") is not False
        ):
            raise ValueError("Wave3 original metadata/native/history/privacy changed")
        sid = s["original_repository_study_id"]
        if sid in seen:
            raise ValueError("Original source-ID pseudoreplication")
        seen.add(sid)
        status = row["source_science_status"]
        if status not in (
            "SOURCE_GPS_HOLD", "SOURCE_FILE_METADATA_HOLD", 
            "HUMAN_LOCATION_PRIVACY_HOLD", "SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE",
            "SOURCE_NO_NATIVE_Q_ESTIMABLE",
        ):
            raise ValueError("Unknown Wave3 source status")
        if ((s["source_status"] == "SOURCE_METADATA_HOLD") !=
            (status == "SOURCE_FILE_METADATA_HOLD")):
            raise ValueError("Original source-file HOLD was reclassified")
        native = row.get("native_group_results")
        if not isinstance(native, list) or len(native) != 1:
            raise ValueError("Original Wave3 native registered source count changed")
        if native[0].get("taxon") != s["registered_native_taxa"][0]:
            raise ValueError("Native taxon changed after event access")
        q = native_mean_if_supported(native[0].get("q"))
        day = native_mean_if_supported(native[0].get("daylength"))
        if ((q is not None) != (status == "SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE")):
            raise ValueError("Source q support/classification is inconsistent")
        if (status in ("SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE","SOURCE_NO_NATIVE_Q_ESTIMABLE")
            and row["original_file_MD5_bytes_verified"] is not True):
            raise ValueError("Original GPS parsed without verified source MD5")
        if status in ("SOURCE_FILE_METADATA_HOLD","HUMAN_LOCATION_PRIVACY_HOLD"):
            if row["original_file_MD5_bytes_verified"] or row["original_source_GPS_event_rows_opened"] != 0:
                raise ValueError("No-contact source HOLD opened provider GPS")
        if row["source_science_status"] == "SOURCE_GPS_HOLD" and (
            q is not None or day is not None
        ):
            raise ValueError("Source-wide fail-closed retained partial biology")
        new.append({
            "wave": "wave3",
            "source_original_repository_study_id": sid,
            "source_DOI_labels": s["original_DOI_labels"],
            "source_stratum": "single_native",
            "original_source_status": status,
            "registered_native_count": 1,
            "native_q_supported": int(q is not None),
            "native_daylength_supported": int(day is not None),
            "complete_source_native_q": q is not None,
            "complete_source_daylength": day is not None,
            "source_q_mean_equal_registered_native": q,
            "source_daylength_mean_equal_registered_native": day,
            "joint_native_month_covariance_status": row["native_joint_calendar_covariance"]["status"],
            "source_MD5_confirmed": row["original_file_MD5_bytes_verified"],
            "original_source_event_row_count": row["original_source_GPS_event_rows_opened"],
            "source_hold_stage": row["source_hold_stage"],
            "source_hold_error_type": row["source_hold_reason"],
        })
    combined = original + new
    assert len(combined) == 51 and len(seen) == 51
    q_old = [r["source_q_mean_equal_registered_native"] for r in original
             if r["complete_source_native_q"]]
    q_new = [r["source_q_mean_equal_registered_native"] for r in new
             if r["complete_source_native_q"]]
    q_both = q_old + q_new
    day_old = [r["source_daylength_mean_equal_registered_native"] for r in original
               if r["complete_source_daylength"]]
    day_new = [r["source_daylength_mean_equal_registered_native"] for r in new
               if r["complete_source_daylength"]]
    if len(q_old) != 10 or len(day_old) != 8:
        raise ValueError("Original 35-source q/daylength measurability changed")
    src_complete_q = sum(r["complete_source_native_q"] for r in combined)
    src_complete_day = sum(r["complete_source_daylength"] for r in combined)
    count_sign = lambda arr: {
        "positive": sum(x > 0 for x in arr),
        "negative": sum(x < 0 for x in arr),
        "exact_zero": sum(x == 0 for x in arr)
    }
    result = {
        "version": "RSDS_437_51_SOURCE_WAVE3_ORIGINAL_GPS_FIRST_META_V1",
        "stage": "FIRST_SOURCE_NATIVE_ORIGINAL_GPS_OUTCOMES_NOT_GENERALIZABLE_WORLD_META",
        "original_35_source_science_git_blob_unchanged": c["original_previous_35_source_result_git_blob"],
        "wave3_original_file_metadata_git_blob_unchanged": c["source_provider_file_metadata_git_blob"],
        "original_catalog_packages": 437,
        "original_GPS_GNSS_packages": 259,
        "original_metadata_priority_frame_sources": 120,
        "original_source_study_units": 51,
        "original_previous_source_study_units": 35,
        "new_wave3_source_study_units": 16,
        "original_source_study_unit_all_registered_native_rows": 66,
        "new_wave3_study_status_counts": dict(sorted(Counter(x["original_source_status"] for x in new).items())),
        "source_complete_all_native_q": src_complete_q,
        "source_complete_all_native_daylength": src_complete_day,
        "source_incomplete_q_count": 51 - src_complete_q,
        "original_35_source_study_rows_unmodified": original,
        "new_wave3_original_source_study_rows": new,
        "all_51_source_study_rows": combined,
        "source_respondent_only_q_descriptive": {
            "previous_35_source_respondents": len(q_old),
            "new_wave3_source_respondents": len(q_new),
            "all_51_source_respondents": len(q_both),
            "previous_35_source_sign_counts": count_sign(q_old),
            "new_wave3_source_sign_counts": count_sign(q_new),
            "all_51_source_sign_counts": count_sign(q_both),
            "new_wave3_unweighted_source_mean_q_if_any": mean(q_new) if q_new else None,
            "all_51_unweighted_respondent_source_mean_q_if_any": mean(q_both) if q_both else None,
            "not_source_frame_population_effect_or_night_activity_global_mean": True,
            "wave3_source_quota_chosen_after_previous_outcomes": True,
        },
        "no_original_GPS_tracks_IDs_exact_dates_coordinates_or_individual_effects_exported": True,
        "previous_48_source_sunlight_signed_external_not_replicated": True,
        "previous_12_ERA5_raw_SHA_HOLD_unchanged": True,
        "scientific_ceilings": [
            "Time observations are original GPS straight-line displacement, not duration of activity.",
            "Wave3 source quota was adaptively chosen after previous outcome and fixed rank identities before current outcome; global wildlife sample is not randomized.",
            "Respondent source q distribution is selected by original measurement process and cannot be unbiased extrapolation to all 51 or 120 or world wildlife.",
            "Original 35 source output is pinned and unchanged. No source, native taxon, failed source, or unavailable effect was silently omitted.",
            "Source overlaps in actual geography and managed/domestic status may remain; native study ID independence is not global ecological independence.",
            "No new solar-relative heading effects or independent sunlight direction replication were evaluated.",
        ],
    }
    return result


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--directory", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    files = sorted(a.directory.glob("wave3_source_*.json"))
    if len(files) != 16:
        raise ValueError("All sixteen Wave3 result files including holds must be present")
    rows = [json.loads(x.read_text(encoding="utf-8")) for x in files]
    output = build(rows)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(
        json.dumps(output, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print("RSDS_WAVE3_51_SOURCE_META", json.dumps({
        "all_source_studies": 51,
        "new_wave3_source_studies": 16,
        "new_status_counts": output["new_wave3_study_status_counts"],
        "complete_q": output["source_complete_all_native_q"],
        "complete_daylength": output["source_complete_all_native_daylength"],
        "new_sign_counts": output["source_respondent_only_q_descriptive"]["new_wave3_source_sign_counts"],
        "all_sign_counts": output["source_respondent_only_q_descriptive"]["all_51_source_sign_counts"],
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
