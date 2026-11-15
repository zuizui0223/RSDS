"""DOI-grain diagnostic discrimination of original q measurability.

Uses the FIRST frozen 35-source original GPS observation-gate aggregate.
NO raw GPS event, no q-effect re-estimation. AUC is the empirical probability
that a source selected as q-complete has a larger *measurement* score than an
incomplete source. This is post-outcome descriptive discrimination among
selected parsed sources, NOT external predictive accuracy or a biological
effect. Source-group labels are outcomes of the preexisting observation rule.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from statistics import mean,median

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned

FIRST=("research/RSDS_35_SOURCE_GPS_OBSERVATION_FUNNEL_FIRST_RESULT_20261115.json",
       "49dcf5c798f8eee788c8001906a1fe0e2c195086")

METRICS={
 "original_adjacent_intervals":
   "source_total_original_adjacent_source_fix_intervals",
 "accepted_speed_intervals":
   "source_total_accepted_finite_speed_max60kmh_intervals",
 "time_0p5_to_2h_fraction":
   "source_native_equal_original_time_0p5_to_2h_over_adjacent",
 "same_solar_date_daypart_given_time_fraction":
   "source_native_equal_same_solar_date_daypart_over_time_eligible",
 "full_original_segment_acceptance_fraction":
   "source_native_equal_total_accepted_over_adjacent",
 "UTC_even_hour_fraction":
   "source_native_equal_UTC_even_hour_fraction",
 "original_tagged_animals":
   "source_total_original_tags",
 "qualified_q_tags_definition_aware":
   "source_total_tags_q_eligible",
}


def rank_probability(group):
    """Nonparametric AUC at the ORIGINAL source study-ID grain."""
    plus=[float(r["score"]) for r in group if r["q_complete"]]
    minus=[float(r["score"]) for r in group if not r["q_complete"]]
    if len(plus)<2 or len(minus)<2:
        return {"status":"UNAVAILABLE_SMALL_SOURCE_STRATUM",
                "q_complete_sources":len(plus),
                "q_incomplete_sources":len(minus),
                "pairwise_source_comparisons":len(plus)*len(minus),
                "AUC_HIGHSCORE_PREDICTS_MEASURABLE":None}
    n=len(plus)*len(minus)
    successes=sum(float(p>m)+.5*float(p==m) for p in plus for m in minus)
    return {
       "status":"SOURCE_LEVEL_POST_OUTCOME_DESCRIPTIVE",
       "q_complete_sources":len(plus),
       "q_incomplete_sources":len(minus),
       "pairwise_source_comparisons":n,
       "AUC_HIGHSCORE_PREDICTS_MEASURABLE":successes/n,
       "median_score_q_complete":median(plus),
       "median_score_q_incomplete":median(minus),
       "score_range_q_complete":[min(plus),max(plus)],
       "score_range_q_incomplete":[min(minus),max(minus)]
    }


def run():
    data=pinned(*FIRST)
    rows=data["individual_original_source_diagnostic_rows"]
    if (data["source_original_registered_native_rows"]!=50
        or len(rows)!=35 or data["original_source_funnel_replayed"]!=25
        or data["first_source_study_q_complete"]!=10
        or data["original_human_GPS_source_opened"] is not False
        or data["original_48_solar_negative_external_NONREPLICATION_UNCHANGED"]!=
          "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
        or data["old_12_ERA5_raw_source_versions_HOLD_UNCHANGED"] is not True):
        raise ValueError("Original 35-source observed admissibility or null changed")
    parsed=[r for r in rows if r["native_observation_funnel_summary"] is not None]
    if len(parsed)!=25 or len({r["source_study_ID"] for r in parsed})!=25:
        raise ValueError("Parsed original 25 source IDs/independence changed")
    if sum(bool(r["complete_native_q"]) for r in parsed)!=10:
        raise ValueError("Native complete-source q eligibility changed")
    by_metric={}
    for label,field in METRICS.items():
        raw=[]
        for r in parsed:
            v=r["native_observation_funnel_summary"][field]
            if v is None:
                continue
            if not isinstance(v,(int,float)) or not math.isfinite(v) or v<0:
                raise ValueError("Invalid original source-level q measurement score")
            raw.append({
               "source_DOI":r["source_DOI"],
               "source_study_ID":r["source_study_ID"],
               "wave":r["wave"],
               "stratum":r["stratum"],
               "q_complete":bool(r["complete_native_q"]),
               "score":v,
            })
        overall=rank_probability(raw)
        waves={w:rank_probability([r for r in raw if r["wave"]==w])
               for w in ("wave1","wave2")}
        by_metric[label]={
           "original_source_technical_count_score_key":field,
           "source_study_N_with_defined_score":len(raw),
           "source_unique_study_IDs":len({r["source_study_ID"] for r in raw}),
           "all_parsed_sources":overall,
           "within_original_wave_descriptive":waves,
           "not_validated_out_of_source_prediction_or_selection_correction":True
        }
    return {
        "version":"RSDS_GPS_SOURCE_Q_MEASURABILITY_RANK_DISCRIMINATION_V1",
        "stage":"POST_GPS_FIRST_OUTCOME_SOURCE_LEVEL_SENSITIVITY_ONLY",
        "original_35_source_frame":35,
        "GPS_parsed_original_sources":25,
        "GPS_native_source_q_complete":10,
        "original_replayed_source_native_rows":34,
        "original_measurement_scores":by_metric,
        "important_negative_control":"qualified_q_tags_definition_aware is an algebraically close component of the endpoint and is NOT an independent predictor; do not promote its AUC into mechanistic predictive transfer",
        "scientific_claim_ceiling":[
           "AUC is descriptive within 25 previously source-selected GPS-parsed DOIs; its pairwise comparisons do not represent 10 x 15 independent population replications.",
           "The q completeness outcome is defined by original tag/year sampling; tag eligibility itself is tautologically predictive.",
           "Original UTC time gate ratio and total accepted interval ratio are source measurements, with geography, habitat, schedule and biological geometry still intertwined.",
           "No observational category is randomly assigned, no out-of-study validation performed, and one human original source remains unaccessed.",
           "Earlier 48-source shortwave direction external negative result NOT REPLICATED and 12 original ERA5 raw-source HOLDS unchanged."
        ],
        "original_48_sunlight_external_negative_nonreplication_UNCHANGED":True,
        "original_12_ERA5_raw_JSON_HOLD_UNCHANGED":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    r=run()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print("RSDS_Q_MEASURABILITY_SOURCE_AUC",json.dumps({
       k:{
        "n":v["source_study_N_with_defined_score"],
        "overall":v["all_parsed_sources"].get("AUC_HIGHSCORE_PREDICTS_MEASURABLE"),
        "wave1":v["within_original_wave_descriptive"]["wave1"].get("AUC_HIGHSCORE_PREDICTS_MEASURABLE"),
        "wave2":v["within_original_wave_descriptive"]["wave2"].get("AUC_HIGHSCORE_PREDICTS_MEASURABLE")
       } for k,v in r["original_measurement_scores"].items()
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
