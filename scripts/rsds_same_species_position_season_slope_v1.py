"""Post-outcome source-native temporal-position vs seasonal slope comparison.

Only reads two immutable aggregated original-GPS output JSONs. Distinct DOI
does not prove independent animal populations. Nominal original tagged-animal
bootstrap source intervals are not a test of source heterogeneity.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from collections import Counter
from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned

SPEC=HERE/"research/RSDS_SAME_SPECIES_POSITION_SEASON_SLOPE_V1_CONTRACT.json"
VERSION="RSDS_SAME_SPECIES_TIME_NICHE_POSITION_VS_SEASONAL_SLOPE_POSTOUTCOME_V1"

def supported_sign(endpoint):
    if endpoint is None or endpoint.get("available") is not True:
        return None
    s=endpoint.get("nominal_source_CI_sign")
    if s=="POSITIVE" or s=="NEGATIVE":
        return s
    if s=="INCLUDES_ZERO":
        return None
    raise ValueError("Unknown original source-native interval sign")

def frozen_sources():
    c=json.loads(SPEC.read_text(encoding="utf-8"))
    if (c.get("version")!=VERSION or
        c.get("prior_63_species_source_atlas_git_blob")!="0cb517ae48d80e503cd2ed5acf150abc047a4352" or
        c.get("frozen_six_target_first_q_git_blob")!="8b254ea8a155098adeac30e70f990575d3b624c2" or
        c.get("first_source_original_q_Action_id")!=43708512175 or
        c.get("first_source_original_q_freeze_Action_id")!=43709239579 or
        c.get("original_followup_study_IDs")!=6 or
        c.get("first_q_supported_sources")!=3 or
        c.get("first_daylength_supported_sources")!=3 or
        c.get("historical_48_source_solar_direction_external_nonreplication") is not True or
        c.get("historical_12_ERA5_raw_SHA_source_HOLD") is not True or
        c.get("target_species_selected_after_prior_results") is not True or
        c.get("no_raw_GPS_access_or_reinterpretation") is not True or
        c.get("confirmatory_inference") is not False):
        raise ValueError("Original historic source pin or claim ceiling changed")
    prior=pinned(c["prior_63_species_source_atlas"],
                 c["prior_63_species_source_atlas_git_blob"])
    new=pinned(c["frozen_six_target_first_q"],c["frozen_six_target_first_q_git_blob"])
    if (prior["distinct_original_archived_source_DOIs_and_study_IDs"]!=63 or
        prior["all_registered_original_source_native_rows"]!=88 or
        new["original_followup_source_DOIs_all"]!=6 or
        new["new_original_source_q_supported"]!=3 or
        new["new_original_source_daylength_supported"]!=3 or
        new["new_source_effect_repeatability_nominal_verdict_counts"]!={"NOMINAL_SUPPORTED_SAME_SIGN":3,"UNAVAILABLE_OR_ORIGINAL_SOURCE_HOLD":3} or
        new["previous_original_48_solar_direction_signed_heldout_nonreplication"] is not True):
        raise ValueError("Original 63-source/first six source biological output changed")
    return c,prior,new

def analyze():
    c,prior,new=frozen_sources()
    by_species={r["scientific_name"]:r for r in
                prior["same_species_source_study_q_comparison_candidates"]}
    if len(by_species)!=6:
        raise ValueError("Old repeated source-native species results lost")
    out=[]
    for row in new["all_six_original_source_taxon_results_including_HOLD"]:
        sp=row["scientific_name"]
        old=by_species[sp]["source_native_rows"]
        prev=[r for r in old if r["original_native_q"]["available"]]
        if len(prev)<2 or row["original_source_DOI"] in {r["original_DOI"] for r in old}:
            raise ValueError("New original source not previously unseen 63-source DOI")
        newq=row["original_new_source_q"]
        newday=row["original_new_source_q_daylength_slope"]
        prior_q_sign=[supported_sign(r["original_native_q"]) for r in prev]
        old_days=[r for r in old if r["original_native_daylength"]["available"]]
        prior_d_sign=[supported_sign(r["original_native_daylength"]) for r in old_days]
        n_q=supported_sign(newq)
        n_d=supported_sign(newday)
        position_agree=(n_q is not None and
                        len(prev)>=2 and all(x==n_q for x in prior_q_sign))
        opposite=(n_d is not None and
                  any(x in ("POSITIVE","NEGATIVE") and x!=n_d for x in prior_d_sign))
        same_slope=(n_d is not None and
                    any(x==n_d for x in prior_d_sign))
        out.append({
            "original_source_DOI":row["original_source_DOI"],
            "scientific_name":sp,
            "new_original_source_status":row["source_status"],
            "old_source_native_q_available":len(prev),
            "new_native_q":newq,
            "old_q_supported_CI_signs":prior_q_sign,
            "source_niche_position_same_nominal_supported_direction":bool(position_agree),
            "new_native_daylength_slope":newday,
            "old_daylength_supported_CI_signs":prior_d_sign,
            "new_slope_supported_opposite_to_at_least_one_old_source":bool(opposite),
            "new_slope_supported_same_as_at_least_one_old_source":bool(same_slope),
            "prior_same_source_tag_population_independence_proven":False,
            "prior_original_GPS_exposure":row["known_previous_original_GPS_exposure_status"],
            "is_genuinely_untouched_original_GPS_external_replication":False,
        })
    if len(out)!=6 or sum(x["source_niche_position_same_nominal_supported_direction"] for x in out)!=3:
        raise ValueError("Original first q cross-source candidate positions changed")
    if sum(x["new_slope_supported_opposite_to_at_least_one_old_source"] for x in out)!=2:
        raise ValueError("Two newly observed opposite seasonal slope sources changed")
    return {
        "version":VERSION,
        "stage":"POST_FIRST_FOLLOWUP_SOURCE_OUTCOME_ECOLOGICAL_AXIS_DIAGNOSTIC",
        "earlier_source_DOI":63,
        "previous_source_native_taxon_rows":88,
        "new_original_source_study_IDs":6,
        "new_original_source_q_respondents":3,
        "new_original_source_daylength_respondents":3,
        "new_q_same_all_prior_supported_direction":3,
        "new_q_opposite_to_prior_supported_direction":0,
        "new_q_HOLD_or_incomplete":3,
        "new_supported_daylength_opposite_to_any_prior_supported_source":2,
        "new_supported_daylength_same_as_any_prior_supported_source":0,
        "original_new_source_axis_rows_including_HOLD":out,
        "prior_GPS_source_independence_not_verified":True,
        "original_tagged_animal_nominal95_intervals_not_crosssource_tests":True,
        "post_outcome_and_species_target_adaptive":True,
        "original_raw_GPS_reopened":False,
        "old_48_source_external_solar_direction_not_replicated":True,
        "old_12_source_ERA5_raw_SHA_HOLD":True,
        "limits":"These counts are descriptive of selected archived DOI/tag-bootstrap intervals, not rates over world populations. A same-sign q interval and an opposite slope interval in two source studies are not a calibrated study-difference CI or causal seasonal response. Source-level tag overlap and telemetry/calendar domain may differ."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    j=analyze()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(j,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_POSITION_SEASON_SLOPE_POSTOUTCOME",json.dumps({
        "source_n":j["new_original_source_study_IDs"],
        "supported_q_same":j["new_q_same_all_prior_supported_direction"],
        "opposite_daylight_slope":j["new_supported_daylength_opposite_to_any_prior_supported_source"],
        "nonresponding_q":j["new_q_HOLD_or_incomplete"],
        "cases":[{"species":r["scientific_name"],"same_position":r["source_niche_position_same_nominal_supported_direction"],
                  "opposite_slope":r["new_slope_supported_opposite_to_at_least_one_old_source"]}
                  for r in j["original_new_source_axis_rows_including_HOLD"]]
    },sort_keys=True))

if __name__=="__main__":
    main()
