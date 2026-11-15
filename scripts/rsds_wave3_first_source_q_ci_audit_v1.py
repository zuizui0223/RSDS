"""Post-first-outcome source-level q bootstrap-interval audit from pinned artifacts.

NEVER opens raw GPS or Movebank provider files: reads existing 16 source-
aggregate results from first run only, and the immutable first 51-source
study-level result committed to GitHub. This is descriptive, post-outcome.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path
from statistics import mean

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned

CONTRACT = HERE / "research/RSDS_WAVE3_FIRST_ORIGINAL_SOURCE_Q_CI_AUDIT_V1_CONTRACT.json"
VERSION = "RSDS_WAVE3_ORIGINAL_SOURCE_Q_CI_POSTOUTCOME_AUDIT_V1"
ORIGINAL_NEW = "RSDS_WAVE3_16_ORIGINAL_GPS_FIRST_SOURCE_Q_DAYLENGTH_V1"


def sign(value):
    return "POSITIVE" if value > 0 else "NEGATIVE" if value < 0 else "EXACT_ZERO"


def bootstrap_interval_status(effect):
    if effect is None:
        return {"available": False, "mean": None, "n_individuals": 0,
                "bootstrap_95_percentile": None, "nominal_point_sign": "UNAVAILABLE",
                "nominal_bootstrap_sign": "UNAVAILABLE"}
    n = effect.get("n_individuals")
    x = effect.get("mean")
    ci = effect.get("animal_bootstrap_95_percentile")
    if (not isinstance(n, int) or n < 5
        or not isinstance(x, (float, int)) or not math.isfinite(x)
        or not isinstance(ci, list) or len(ci) != 2
        or not all(isinstance(y, (float, int)) and math.isfinite(y) for y in ci)
        or ci[0] > ci[1]):
        raise ValueError("Incorrect original individual bootstrap effect or source support")
    interval_sign = (
        "NOMINAL_POSITIVE" if ci[0] > 0 else
        "NOMINAL_NEGATIVE" if ci[1] < 0 else
        "INCLUDES_ZERO"
    )
    return {"available": True, "mean": float(x), "n_individuals": n,
            "bootstrap_95_percentile": [float(ci[0]), float(ci[1])],
            "nominal_point_sign": sign(float(x)),
            "nominal_bootstrap_sign": interval_sign}


def analyze(directory):
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (c.get("version") != VERSION
        or c.get("first_successful_wave3_GPS_action_run_id") != 43702456789
        or c.get("original_wave3_first_51_meta_git_blob") !=
           "4644b8ed38ae74dae345df8d90ea182f5dce4ee3"
        or c.get("new_16_source_studies") != 16
        or c.get("new_source_native_q_supported") != 6
        or c.get("no_GPS_source_event_rerun_or_download") is not True):
        raise ValueError("Original source first-artifact contract changed")
    d = pinned(c["original_wave3_first_51_meta_path"],
               c["original_wave3_first_51_meta_git_blob"])
    meta = d["new_wave3_original_source_study_rows"]
    if (len(meta) != 16 or d["original_source_study_units"] != 51
        or d["source_complete_all_native_q"] != 16
        or d["source_complete_all_native_daylength"] != 12
        or d["original_35_source_science_git_blob_unchanged"] != c["original_35_result_git_blob"]
        or d["previous_48_source_sunlight_signed_external_not_replicated"] is not True):
        raise ValueError("Previously original first meta and endpoint changed")
    expected_artifacts = c["original_16_source_artifacts"]
    if (len(expected_artifacts) != 16
        or len({x["artifact_id"] for x in expected_artifacts}) != 16
        or len({x["digest"] for x in expected_artifacts}) != 16):
        raise ValueError("Original 16 first artifact receipts incomplete")
    files = sorted(directory.glob("wave3_source_*.json"))
    if len(files) != 16:
        raise ValueError("Require all sixteen first original source artifacts")
    actual = {}
    for p in files:
        row = json.loads(p.read_text(encoding="utf-8"))
        i = row.get("wave3_index")
        if not isinstance(i, int) or not 0 <= i < 16 or i in actual:
            raise ValueError("Duplicate/out-of-source study index")
        actual[i] = row
    if set(actual) != set(range(16)):
        raise ValueError("Original first source slots incomplete")
    rows = []
    for index, expected in enumerate(meta):
        raw = actual[index]
        if (raw.get("version") != ORIGINAL_NEW
            or raw.get("wave3_original_source_rank_1_based") != index + 25
            or raw.get("source_original_repository_study_id") != expected["source_original_repository_study_id"]
            or raw.get("source_original_DOI_labels") != expected["source_DOI_labels"]
            or raw.get("source_science_status") != expected["original_source_status"]
            or raw.get("original_file_MD5_bytes_verified") != expected["source_MD5_confirmed"]
            or raw.get("original_source_GPS_event_rows_opened") != expected["original_source_event_row_count"]
            or raw.get("original_raw_GPS_dates_tag_IDs_coordinates_exported") is not False
            or raw.get("human_mobility_GPS_data_opened") is not False
            or len(raw.get("native_group_results") or []) != 1):
            raise ValueError("Original first source effect mismatches frozen source study")
        native = raw["native_group_results"][0]
        if native["taxon"] != raw["original_registered_native_taxa"][0]:
            raise ValueError("Original source native species changed")
        q = bootstrap_interval_status(native.get("q"))
        daylight = bootstrap_interval_status(native.get("daylength"))
        for key, x, positive in (
            ("q",q,"native_q_supported"),
            ("daylength",daylight,"native_daylength_supported"),
        ):
            frozen_value = expected["source_q_mean_equal_registered_native"] if key=="q" else expected["source_daylength_mean_equal_registered_native"]
            if ((x["available"] != bool(expected[positive]))
                or (x["mean"] is not None and (
                    frozen_value is None or abs(x["mean"]-frozen_value)>1e-12))
                or (x["mean"] is None and frozen_value is not None)):
                raise ValueError("Frozen original source effect changed")
        rows.append({
            "wave3_index": index,
            "original_repository_study_id": expected["source_original_repository_study_id"],
            "source_DOI_labels": expected["source_DOI_labels"],
            "original_registered_native_taxa": raw["original_registered_native_taxa"],
            "original_source_status": expected["original_source_status"],
            "original_file_MD5_verified": expected["source_MD5_confirmed"],
            "native_q": q,
            "native_daylength": daylight,
            "both_q_and_daylength_available": q["available"] and daylight["available"],
            "q_daylength_point_signs_differ": (
                sign(q["mean"]) != sign(daylight["mean"])
                if q["available"] and daylight["available"] else None
            ),
        })
    nq = [r["native_q"] for r in rows if r["native_q"]["available"]]
    nd = [r["native_daylength"] for r in rows if r["native_daylength"]["available"]]
    if len(nq) != 6 or len(nd) != 4:
        raise ValueError("Wave3 source-native first q/daylight support changed")
    qvals = [x["mean"] for x in nq]
    loo = [mean(qvals[:i]+qvals[i+1:]) for i in range(len(qvals))]
    both = [r for r in rows if r["both_q_and_daylength_available"]]
    if len(both) != 4:
        raise ValueError("Same-source original q/daylight pair count changed")
    return {
        "version": VERSION,
        "stage": "POST_FIRST_ORIGINAL_SOURCE_OUTCOMES_CI_SENSITIVITY_NO_NEW_GPS",
        "first_source_run_id": c["first_successful_wave3_GPS_action_run_id"],
        "previous_51_meta_git_blob_unchanged": c["original_wave3_first_51_meta_git_blob"],
        "source_original_study_units": 51,
        "new_original_study_units": 16,
        "native_original_q_respondents": len(nq),
        "native_original_daylength_respondents": len(nd),
        "source_q_point_signs": dict(sorted(Counter(x["nominal_point_sign"] for x in nq).items())),
        "source_q_bootstrap_nominal_signs": dict(sorted(Counter(x["nominal_bootstrap_sign"] for x in nq).items())),
        "source_daylength_point_signs": dict(sorted(Counter(x["nominal_point_sign"] for x in nd).items())),
        "source_daylength_bootstrap_nominal_signs": dict(sorted(Counter(x["nominal_bootstrap_sign"] for x in nd).items())),
        "source_joint_q_daylength_pairs": len(both),
        "source_joint_q_daylength_point_sign_disagreement_count": sum(x["q_daylength_point_signs_differ"] for x in both),
        "new_respondent_unweighted_mean_q": mean(qvals),
        "new_respondent_leave_one_source_out_mean_range": [min(loo),max(loo)],
        "new_respondent_leave_one_source_out_all_negative": all(x<0 for x in loo),
        "new_16_source_original_native_ci_and_support_rows": rows,
        "descriptive_only_no_biological_world_prevalence_or_causal_claim": True,
        "source_bootstrap_not_multisource_adjusted": True,
        "first_artifact_digest_roster": expected_artifacts,
        "first_combined_meta_artifact_receipt": c["original_51_meta_first_artifact"],
        "original_48_source_external_sunlight_heading_direction_replicated": False,
        "older_12_ERA5_raw_SHA_hold_unchanged": True,
        "no_original_GPS_source_rows_animal_IDs_coordinates_or_exact_times_exported": True,
        "limitations": [
            "Original animal bootstrap 95% intervals describe each source only and have no source-selection/missingness or multiple-testing correction.",
            "Wave3 quota selected after observing Wave1 and Wave2 data; this is explicitly post-first-outcome exploratory analysis.",
            "Observations are straight-line step displacements, not time spent active or causal sunlight avoidance.",
            "Source 51 studies / 16 new selection is not a probability-sample world ecosystem, and managed/domestic source domains may remain."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    j=analyze(a.directory)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(j,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_WAVE3_SOURCE_BOOTSTRAP_CI_AUDIT",json.dumps({
        "new_q_n":j["native_original_q_respondents"],
        "nominal_q_signs":j["source_q_bootstrap_nominal_signs"],
        "new_daylight_n":j["native_original_daylength_respondents"],
        "nominal_daylight_signs":j["source_daylength_bootstrap_nominal_signs"],
        "q_daylength_pairs":j["source_joint_q_daylength_pairs"],
        "q_daylight_point_sign_disagreement":j["source_joint_q_daylength_point_sign_disagreement_count"],
        "new_leave_one_source_out_range":j["new_respondent_leave_one_source_out_mean_range"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
