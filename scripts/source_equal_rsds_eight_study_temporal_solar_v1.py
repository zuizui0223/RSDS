"""Original source-study equal-weight macro sensitivity, not a global parameter.

Reads only prior eight source study / 12 native taxon FIRST ecological aggregate
JSON. Every source, including original egret 200MB metadata HOLD, is retained.
Comparison separates (1) study×taxon row equal weights (pseudoreplication
risk) and (2) original field study DOI equal weights, conditional on source
availability. Study-block bootstrap is descriptive conditional on seven/eight
PURPOSIVELY selected original field studies; no global p-value or causal
ecological adaptation inferred.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json"
CONTRACT = ROOT / "research/RSDS_EIGHT_STUDY_SOURCE_EQUAL_MACRO_V1_CONTRACT.json"
RESULT = ROOT / "research/RSDS_EIGHT_STUDY_SOURCE_EQUAL_MACRO_SENSITIVITY_FIRST_RESULT_20261026.json"
EXPECTED_SHA = "2c4670b54fa5b52737993ebf6babc3a7a58b166a"
SEED, DRAW = 20261026, 5000
VARIABLES = {
    "night_day_q": "q",
    "q_daylength_per_hour": "daylength_slope",
    "raw_solar_relative_cos": "sun_raw",
    "joint_circular_covariance": "sun_joint",
}
HISTORIC = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def load_original(root=ROOT):
    p = root / INPUT.relative_to(ROOT)
    raw = p.read_bytes()
    sha = hashlib.sha1(("blob " + str(len(raw)) + "\0").encode() + raw).hexdigest()
    frozen = json.loads((root / CONTRACT.relative_to(ROOT)).read_text(encoding="utf-8"))
    if (sha != EXPECTED_SHA
            or frozen["original_git_blob_sha"] != sha
            or frozen["version"] != "RSDS_ORIGINAL_STUDY_EQUAL_WEIGHT_MACRO_SENSITIVITY_V1"):
        raise ValueError("Original ecological 8-study first result changed or contract not pinned")
    return json.loads(raw.decode("utf-8"))


def interval(arr):
    seq = sorted(arr)
    return [seq[int(.025 * len(seq))],
            seq[min(len(seq)-1, int(.975 * len(seq))-1)]]


def equal_source_bootstrap(values, seed=SEED, draws=DRAW):
    if len(values) < 3:
        return None
    rng = random.Random(seed)
    n = len(values)
    return interval([mean(values[rng.randrange(n)] for _ in range(n))
                     for _ in range(draws)])


def summarize_axis(studies, var, i):
    grouped = []
    frame = []
    for source in studies:
        vals = [row[VARIABLES[var]]["mean"] for row in source["original_native_taxa"]
                if row[var]["mean"] is not None]
        if len(vals) != source["original_taxa_supported_per_axis"][var]:
            raise ValueError("Source supported/unknown effect count changed")
        record = {"study_doi":source["study_doi"],
                  "taxon_rows_in_original_study":source["original_taxon_count"],
                  "taxa_with_defined_effect":len(vals),
                  "study_equal_taxon_mean":mean(vals) if vals else None,
                  "study_mean_sign":"positive" if vals and mean(vals)>0 else
                                    "negative" if vals and mean(vals)<0 else
                                    "zero" if vals else "UNAVAILABLE"}
        grouped.append(record)
        if vals:
            frame.append((source["study_doi"], mean(vals), vals))
    src = [v for _,v,_ in frame]
    naïve = [x for _,_,xs in frame for x in xs]
    source_n = len(src)
    if source_n < 2:
        raise ValueError("Insufficient fixed original study units for source-equal sensitivity")
    actual = mean(src)
    naive = mean(naïve)
    loo = [mean(src[:j] + src[j+1:]) for j in range(source_n)]
    return {
        "fixed_original_studies_in_denominator":len(studies),
        "original_studies_with_supported_axis":source_n,
        "original_source_taxon_rows_with_supported_axis":len(naïve),
        "equal_source_mean":actual,
        "equal_native_taxon_row_mean_NOT_INDEPENDENT_STUDIES":naive,
        "equal_source_minus_equal_taxon_row":actual-naive,
        "equal_source_95percentile_study_resampling_CONDITIONAL_NONPROBABILITY":
            None if var == "joint_circular_covariance" else
            equal_source_bootstrap(src, seed=SEED + 19*i),
        "source_resampling_interval_available":
            var != "joint_circular_covariance",
        "uncertainty_ceiling":
            ("OLD_SIX_JOINT_ESTIMATES_ARE_POINT_ONLY_NO_TAG_INTERVAL"
             if var == "joint_circular_covariance" else
             "SOURCE_RESAMPLING_CONDITIONAL_ON_PURPOSIVELY_SELECTED_FIELD_STUDIES_ONLY"),
        "leave_one_original_study_out_equal_source_mean_range":[min(loo),max(loo)],
        "supported_original_source_signs":{
            "positive":sum(x>0 for x in src),
            "negative":sum(x<0 for x in src),
            "exact_zero":sum(x==0 for x in src),
        },
        "first_original_study_QA":grouped,
    }


def run(root=ROOT):
    data = load_original(root)
    if (data["original_study_denominator"] != 8
            or data["source_taxon_rows"] != 12
            or data["source_taxa_with_q"] != 11
            or data["source_taxa_with_solar_raw"] != 10
            or data["old_48_radiation_study_external_negative_direction_test_UNCHANGED"] != HISTORIC
            or not data["old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED"]):
        raise ValueError("Original registered ecological source frame/heldout history changed")
    sources = defaultdict(list)
    for r in data["all_sources_including_egret_original_metadata_HOLD"]:
        sources[r["source_study_doi"]].append(r)
    if len(sources) != 8:
        raise ValueError("DOI identity no longer unique across eight original field studies")
    study = []
    for doi in sorted(sources):
        taxa = sources[doi]
        if len({t["taxon"] for t in taxa}) != len(taxa):
            raise ValueError("Duplicate native taxon in original field study")
        number = {key:sum(t[field]["mean"] is not None for t in taxa)
                  for key,field in VARIABLES.items()}
        study.append({
            "study_doi":doi, "original_taxon_count":len(taxa),
            "original_native_taxa":taxa,
            "original_taxa_supported_per_axis":number,
        })
    if sorted(s["original_taxon_count"] for s in study) != [1,1,1,1,1,1,2,4]:
        raise ValueError("Original study taxonomic replication frame changed")
    axes = {key:summarize_axis(study,key,i)
            for i,key in enumerate(VARIABLES)}
    if (axes["night_day_q"]["original_studies_with_supported_axis"]!=7
            or axes["raw_solar_relative_cos"]["original_studies_with_supported_axis"]!=6):
        raise ValueError("Earlier incomplete source endpoints not preserved")
    # Source equal input claims reflect differing q/direction eligibility
    # even within the same tagged study. No cross-axis microbehavior inference.
    return {
        "version":"RSDS_EIGHT_STUDY_EQUAL_SOURCE_WEIGHT_ORIGINAL_GPS_SENSITIVITY_V1",
        "date":"2026-10-26",
        "analysis_stage":"POST_FIRST_RESULTS_SOURCE_GRAIN_DESCRIPTIVE_SENSITIVITY",
        "original_source_studies":8,
        "original_source_taxon_rows":12,
        "source_equal_method":"one original study DOI = one unit, equal available native taxon means within each DOI; no GPS fix, animal tag or number-of-taxa weighting",
        "fixed_original_source_sample_frame":{
            "original_studies":len(study),
            "study_native_taxon_counts":{s["study_doi"]:s["original_taxon_count"] for s in study},
            "source_with_only_metadata_HOLD":[s["study_doi"] for s in study
                 if all(v==0 for v in s["original_taxa_supported_per_axis"].values())],
        },
        "source_equal_axes":axes,
        "old_external_shortwave_direction_nonreplication_UNCHANGED":HISTORIC,
        "old_12_ERA5_provider_source_SHA_HOLDS_UNCHANGED":True,
        "scientific_ceiling":[
            "This route is a post-outcome change in descriptive sampling grain, not an independent new animal GPS field study.",
            "Seven/eight original public studies were chosen purposively, not from a random population of world ecosystems. Study-bootstrap percentiles are study-resampling sensitivity only, not global 95% confidence intervals.",
            "A shared original field study with multiple taxa contributes one equal original-study mean; even distinct study IDs can share actual ecological populations.",
            "The nighttime/daytime q proxy estimates observed GPS step displacement speed, not true animal active time.",
            "Changing study vs taxon weights does not causally establish species, guild, temperature, daylight adaptation or geographical heterogeneity mechanisms.",
            "Joint circular covariance point-only old source estimates have no source-specific source-animal CI. Study-resampling of their point estimates does not create missing tag-level uncertainty.",
            "All excluded source endpoints and egret size cap/source HOLD remain explicitly UNAVAILABLE, not ecological zeros.",
            "Historical 48 original study signed radiation nonreplication and 12 ERA5 raw provider SHA holds are unchanged.",
            "No original animal telemetry tracks, tag IDs, exact GPS positions/times or per-tag ecological slopes are redistributed."
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=RESULT)
    args=parser.parse_args()
    result = run()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_SOURCE_EQUAL_MACRO_V1",json.dumps({
        "studies":result["original_source_studies"],
        "taxon_rows":result["original_source_taxon_rows"],
        "q_source_equal":result["source_equal_axes"]["night_day_q"]["equal_source_mean"],
        "q_taxon_equal":result["source_equal_axes"]["night_day_q"]["equal_native_taxon_row_mean_NOT_INDEPENDENT_STUDIES"],
        "source_signed":result["source_equal_axes"]["night_day_q"]["supported_original_source_signs"],
        "prior_nonrep":result["old_external_shortwave_direction_nonreplication_UNCHANGED"],
    },sort_keys=True),flush=True)


if __name__ == "__main__":
    main()
