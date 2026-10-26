"""Summarize new RSDS same-tagged-animal q and solar-direction source results.

Pure post-first ecological analysis from public aggregate receipts: no
additional GPS/ERA5 provider access and no animal/coordinate/time disclosure.
Never substitutes missing source effects with zero or reclassifies earlier
negative external radiation replication.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=ROOT/"research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
PREVIOUS=ROOT/"research/RSDS_CROSS_AXIS_GENERALITY_FIRST_RESULT_20261026.json"
RADIATION=ROOT/"research/SOLAR_RADIATION_TWO_24_SOURCE_NONREPLICATION_DIAGNOSTIC_FIRST_RESULT_20261025.json"
OUTPUT=ROOT/"research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_ECOLOGICAL_SYNTHESIS_20261026.json"
FIXED={"Utah":["Antilocapra americana","Cervus elaphus","Odocoileus hemionus","Ovis canadensis"],
       "Wyoming":["Odocoileus hemionus","Ovis canadensis"]}
NAMES=("sun_cos_per100W","sun_sin_per100W","q_daypart_per100W")
FROZEN_NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def signed_source(summary):
    if summary is None or summary.get("source_mean") is None:
        return "UNAVAILABLE"
    ci=summary.get("source_animal_bootstrap95")
    if not (isinstance(ci,list) and len(ci)==2 and all(isinstance(x,(float,int)) and math.isfinite(x) for x in ci)):
        return "UNAVAILABLE"
    if ci[0]>ci[1] or not ci[0]<=summary["source_mean"]<=ci[1]:
        raise ValueError("Invalid original source-animal descriptive bootstrap")
    if ci[0]>0:return "POSITIVE_SOURCE_BOOTSTRAP"
    if ci[1]<0:return "NEGATIVE_SOURCE_BOOTSTRAP"
    return "UNRESOLVED_SIGN"


def check_complete_original_sources(data,old):
    if data["version"]!="rsds-original-two-studies-same-animal-solar-weather-v0":
        raise ValueError("Not a source-frozen first-run same-animal receipt")
    if data["fixed_original_study_count"]!=2 or data["total_source_native_taxon_rows"]!=6:
        raise ValueError("Two original studies and six native taxon rows must remain")
    if data["frozen_external_radiation_nonreplication_UNCHANGED"]!=FROZEN_NEGATIVE:
        raise ValueError("Previously unreplicated negative Sun direction effect was reclassified")
    src={x["region"]:x for x in data["sources"]}
    if set(src)!=set(FIXED):
        raise ValueError("Original regional studies not retained")
    for name, taxa in FIXED.items():
        x=src[name]
        if x["original_native_taxa"]!=taxa:
            raise ValueError("Original taxon roster changed in source-first result")
        if x["status"] not in ("SAME_ORIGINAL_ANIMAL_WEATHER_FRAME_COMPLETE","ORIGINAL_SOURCE_HOLD"):
            raise ValueError("Unknown original source first state")
        observed=[t["taxon"] for t in x["taxon_results"]]
        if observed!=taxa:
            raise ValueError("Original source first report lacks native taxon row")
        if x["source_original_replayed_before_NEW_weather"]:
            if x["status"]!="SAME_ORIGINAL_ANIMAL_WEATHER_FRAME_COMPLETE":
                raise ValueError("Source replay cannot imply unverified new weather success")
            if x["original_provider_integrity"]["provider_checksum_verified"] is not True:
                raise ValueError("Original provider MD5 not validated")
        elif x["status"]!="ORIGINAL_SOURCE_HOLD":
            raise ValueError("Missing original q/daylight replay must be explicit SOURCE_HOLD")
    if old["sunlight_weather_temporal_observation_overlap"]["preregistered_heldout_decision_UNCHANGED"]!=FROZEN_NEGATIVE:
        raise ValueError("Previously frozen older RSDS radiation verdict changed")
    return src


def analyze(original, previous, rad):
    src=check_complete_original_sources(original,previous)
    if rad["prior_frozen_external_replication_decision_UNCHANGED"]!=FROZEN_NEGATIVE:
        raise ValueError("External negative-direction nonreplication was silently changed")
    taxa=[]
    for region, study in src.items():
        for source_taxon in study["taxon_results"]:
            pair=source_taxon.get("source_animal_summary") or {}
            summaries={name:pair.get(name) for name in NAMES}
            cells=source_taxon.get("selected_contexts") or []
            context_flips={}
            for key,inner in (("sun_cos_per100W","direction_sun_cos"),
                              ("sun_sin_per100W","direction_sun_sin"),
                              ("q_daypart_per100W","daypart_q_vs_irradiance")):
                signed=[signed_source(c.get(inner)) for c in cells]
                context_flips[key]={
                    "eligible_distinct_original_contexts":len(cells),
                    "positive_excluding_zero":signed.count("POSITIVE_SOURCE_BOOTSTRAP"),
                    "negative_excluding_zero":signed.count("NEGATIVE_SOURCE_BOOTSTRAP"),
                    "both_directions_present_in_selected_contexts":(
                        "POSITIVE_SOURCE_BOOTSTRAP" in signed and
                        "NEGATIVE_SOURCE_BOOTSTRAP" in signed),
                }
            row={
                "region":region,"source_doi":study["original_source_doi"],
                "taxon":source_taxon["taxon"],
                "source_original_q_replayed":source_taxon.get("original_source_replay_success",False),
                "source_animal_days_selected_before_weather":sum(c.get("original_qualified_q_days",0) for c in cells),
                "weather_verified_contexts":sum(c.get("ERA5_original_provider_response_sha256") is not None for c in cells),
                "paired_source_individuals":pair.get("paired_source_animals",0),
                "measured_source_axis_coefficients":{
                    key:({"source_mean":s.get("source_mean"),
                           "animal_bootstrap95":s.get("source_animal_bootstrap95"),
                           "eligible_original_animals":s.get("n_source_animals"),
                           "sign_support":signed_source(s)} if s is not None else
                         {"source_mean":None,"animal_bootstrap95":None,
                          "eligible_original_animals":0,"sign_support":"UNAVAILABLE"})
                    for key,s in summaries.items()},
                "within_source_animal_rank_correlation_q_vs_solar_cos":pair.get("within_source_animal_rho_q_vs_sun_cos"),
                "within_source_animal_rank_correlation_q_vs_solar_sin":pair.get("within_source_animal_rho_q_vs_sun_sin"),
                "context_sign_sensitivity":context_flips,
                "taxon_source_status":source_taxon.get("status","UNKNOWN"),
                "original_total_source_studies_represented_by_taxon_row":1,
            }
            taxa.append(row)
    comparator=[]
    for species in ("Odocoileus hemionus","Ovis canadensis"):
        rows={x["region"]:x for x in taxa if x["taxon"]==species}
        if set(rows)!={"Utah","Wyoming"}:
            raise ValueError("Both preselected original native taxa missing")
        diff={}
        for name in NAMES:
            ut=rows["Utah"]["measured_source_axis_coefficients"][name]
            wy=rows["Wyoming"]["measured_source_axis_coefficients"][name]
            diff[name]={
                "Utah":ut,"Wyoming":wy,
                "signed_point_difference_Utah_minus_Wyoming":(
                    ut["source_mean"]-wy["source_mean"] if
                    ut["source_mean"] is not None and wy["source_mean"] is not None
                    else None),
                "both_groups_signed_in_same_direction":(
                    ut["source_mean"]*wy["source_mean"]>0
                    if ut["source_mean"] is not None and wy["source_mean"] is not None
                    else None),
            }
        comparator.append({
            "taxon":species,
            "two_original_source_studies":[rows["Utah"]["source_doi"],rows["Wyoming"]["source_doi"]],
            "original_study_unit_replications":2,
            "weather_time_direction":diff,
            "interpretation":"Two source climate/geography/protocol sampling frames, not two independently random world biomes. Coefficient differences are descriptive; bootstrap source animal intervals not independent causal tests."
        })
    return {
        "version":"rsds-same-tagged-animal-solar-weather-ecological-synthesis-v0",
        "date":"2026-10-26",
        "stage":"POST_FIRST_ORIGINAL_SOURCE_REPLAY_DESCRIPTIVE",
        "original_studies":2,
        "native_taxon_study_rows":6,
        "first_original_full_source_state":{k:v["status"] for k,v in src.items()},
        "frozen_external_radiation_nonreplication_UNCHANGED":FROZEN_NEGATIVE,
        "all_taxon_study_results":taxa,
        "two_study_native_taxon_comparisons":comparator,
        "boundaries":[
            "Only studied original tagged animals eligible for original native day/night source dates and selected source gridyears; missingness is not an ecological zero.",
            "Directional and daypart weather partial associations reuse same source animals and eligible days but target different GPS daytime/nighttime intervals.",
            "Multiple original source gridyear contexts share original animals and years, not independent ecosystems.",
            "Original context-specific weather values may reflect provider version differences on rerun; source JSON SHA preserved.",
            "ER​A5 sun-adjusted radiation anomalies are neither measured body irradiance nor randomized shortwave exposure.",
            "The source-level bootstrap represents selected source animals only; two regions not enough for worldwide macrobiome generality."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=OUTPUT)
    opts=p.parse_args()
    data=analyze(
        json.loads(ORIGINAL.read_text(encoding="utf-8")),
        json.loads(PREVIOUS.read_text(encoding="utf-8")),
        json.loads(RADIATION.read_text(encoding="utf-8"))
    )
    opts.output.parent.mkdir(parents=True,exist_ok=True)
    opts.output.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_ORIGINAL_TWO_REGION_SAME_ANIMAL_SYNTHESIS",
          json.dumps({
             "all_studies":data["original_studies"],
             "all_taxa":data["native_taxon_study_rows"],
             "source_states":data["first_original_full_source_state"],
             "paired_available":{x["region"]+":"+x["taxon"]:x["paired_source_individuals"] for x in data["all_taxon_study_results"]}
          },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
