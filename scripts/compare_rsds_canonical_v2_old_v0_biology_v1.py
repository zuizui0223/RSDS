"""Compare original v0 shortwave coefficients with NEW canonical ERA5 source model.

Original 2026-10-25/26 model effects were frozen, but old ERA5 physical
hourly ARRAYS were not archived, so a byte-hash mismatch cannot establish
weather-value change. This compares *previously saved fitted effects* with
new independently source-physical-SHA-verified fitted effects and documents
same-original-tag two-context descriptive contrasts. NOT an independent
heldout field-study replication, exposure treatment or climate causality.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import HERE

ROOT=HERE
FILES={
"research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json":
    "d8d84ca3583cc26ff373fd929f3f98156ca9b0cb",
"research/RSDS_ERA5_ORIGINAL_SOURCES_CANONICAL_FIRST_RESULT_20261028.json":
    "d39a1ee209ca81a895df6115110759ee4d0acd84",
"research/RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_FIRST_RESULT_20261028.json":
    "ac87bdab689567ed145aa91ed8d24bab8cdce4f1",
}
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
PER_CONTEXT={
"sun_cos":("direction_sun_cos","sun_cos"),
"sun_sin":("direction_sun_sin","sun_sin"),
"daypart_q":("daypart_q_vs_irradiance","daypart_q"),
}
PER_SOURCE={
"sun_cos":("sun_cos_per100W","sun_cos_vs_irradiance"),
"sun_sin":("sun_sin_per100W","sun_sin_vs_irradiance"),
"daypart_q":("q_daypart_per100W","daypart_q_vs_irradiance"),
}
OUTPUT=HERE/"research/RSDS_CANONICAL_ERA5_V2_OLD_V0_EFFECT_REPLAY_AND_TAG_PAIRS_20261028.json"


def load_pinned():
    x={}
    for path,expected in FILES.items():
        raw=(ROOT/path).read_bytes()
        got=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
        if got!=expected:
            raise ValueError("Prior original data result changed: "+path)
        x[path]=json.loads(raw.decode("utf-8"))
    return list(x.values())


def sign_interval(ci):
    if ci is None:
        return "NOT_ESTIMABLE"
    a,b=ci
    if a>0:return "NOMINAL_POSITIVE"
    if b<0:return "NOMINAL_NEGATIVE"
    return "INCLUDES_ZERO"


def run():
    old,physical,new=load_pinned()
    if (old["fixed_original_study_count"]!=2 or old["total_source_native_taxon_rows"]!=6
            or physical["current_provider_contexts_two_response_comparisons"]!=12
            or physical["within_current_session_canonical_UTC_hourly_values_equal"]!=12
            or new["original_study_count"]!=2 or new["original_source_taxon_rows"]!=6
            or new["selected_contexts_reported"]!=12
            or new["current_canonical_weather_contexts_exact_source_SHA_replayed"]!=12
            or not new["historical_12_ERA5_raw_JSON_SHA_HOLDs_UNCHANGED"]
            or new["historic_48_source_negative_radiation_replication_UNCHANGED"]!=NEG):
        raise ValueError("Frozen original GPS, weather provider or scientific source frame changed")
    oldrows={(s["region"],t["taxon"]):t for s in old["sources"]
             for t in s["taxon_results"]}
    if len(oldrows)!=6:
        raise ValueError("Original first native source taxon keys duplicate")
    deltas_context=[]
    deltas_source=[]
    rows=[]
    source_paired={}
    for rec in new["all_fixed_original_study_taxon_results"]:
        key=(rec["region"],rec["taxon"])
        parent=oldrows.pop(key,None)
        if parent is None:
            raise ValueError("Unexpected original study native taxon label")
        if len(rec["original_contexts"])!=2 or len(parent["selected_contexts"])!=2:
            raise ValueError("Original q source context count changed")
        comparisons=[]
        for i,(fresh,prior) in enumerate(zip(rec["original_contexts"],
                                             parent["selected_contexts"]),1):
            if (fresh["original_q_coverage_context_rank"]!=
                    prior["selected_gridyear_original_coverage_rank"] or
                fresh["old_source_q_qualified_days"]!=prior["original_qualified_q_days"] or
                fresh["old_source_q_animals"]!=prior["original_q_qualified_animals"] or
                fresh["old_sun_moving_steps"]!=prior["original_moving_daylight_steps"] or
                fresh["new_provider_physical_content_status"]!="CANONICAL_PHYSICAL_EXACT_SHA_MATCH"):
                raise ValueError("Old original GPS source context changed before reanalysis")
            axes={}
            for axis,(oldkey,newkey) in PER_CONTEXT.items():
                y0=(prior.get(oldkey) or {}).get("source_mean")
                y1=(fresh.get("new_source_effects") or {}).get(newkey,{}).get("source_mean")
                if y0 is None or y1 is None:
                    axes[axis]={"status":"UNESTIMABLE_CONTEXT","old":y0,"new":y1,"delta":None}
                    continue
                d=y1-y0
                deltas_context.append(abs(d))
                axes[axis]={"status":"NUMERICALLY_REPRODUCED_1E-12" if abs(d)<=1e-12 else "DIFFERENT_FIRST_ESTIMATE",
                            "old_source_mean":y0,
                            "new_canonical_source_mean":y1,
                            "new_minus_old_fitted_coefficient":d}
            comparisons.append({"rank_by_original_q_coverage":i,"axes":axes})
        group={}
        for axis,(oldkey,newkey) in PER_SOURCE.items():
            y0=parent["source_animal_summary"][oldkey]["source_mean"]
            y1=rec["source_animal_response_per100W"][newkey]["mean"]
            if y0 is None or y1 is None:
                group[axis]={"status":"UNESTIMABLE","delta":None}
                continue
            d=y1-y0
            deltas_source.append(abs(d))
            group[axis]={"status":"NUMERICALLY_REPRODUCED_1E-12" if abs(d)<=1e-12 else "DIFFERENT_FIRST_ESTIMATE",
                         "old_source_mean":y0,"new_canonical_source_mean":y1,
                         "new_minus_old_fitted_coefficient":d,
                         "new_animal_source_bootstrap95":rec["source_animal_response_per100W"][newkey]["animal_bootstrap95"]}
        pair=rec["paired_tagged_context_response_RANK2_MINUS_RANK1"]
        paired={}
        if pair:
            for axis in ("q","sun_cos","sun_sin"):
                a=pair["per_axis"][axis]
                paired[axis]={
                    "tagged_animals_in_both_contexts_with_axis":a["paired_source_animals"],
                    "rank2_minus_rank1_original_tag_mean":a["rank2_minus_rank1_source_tag_mean"],
                    "descriptive_paired_animal_bootstrap95":a["paired_animal_bootstrap95"],
                    "nominal_source_interval_sign":sign_interval(a["paired_animal_bootstrap95"])
                }
        if len(paired)==0:
            paired={k:{"nominal_source_interval_sign":"NOT_ESTIMABLE"} for k in ("q","sun_cos","sun_sin")}
        for axis,stats in paired.items():
            typ=source_paired.setdefault(axis,{"NOMINAL_POSITIVE":0,"NOMINAL_NEGATIVE":0,"INCLUDES_ZERO":0,"NOT_ESTIMABLE":0})
            typ[stats["nominal_source_interval_sign"]]+=1
        rows.append({
            "region":rec["region"],"native_taxon":rec["taxon"],
            "two_original_q_ranked_contexts":comparisons,
            "old_v0_vs_new_canonical_source_animal_effects":group,
            "same_original_tagged_animal_between_contexts":rec["same_original_tagged_animals_in_two_q_coverage_ranked_contexts"],
            "paired_context_effect_status":rec["paired_tagged_animal_contrast_status"],
            "paired_same_tag_delta_axes":paired
        })
    if oldrows or len(rows)!=6:
        raise ValueError("One original source native taxon missing")
    return {
        "version":"RSDS_V0_V2_ORIGINAL_ERA5_COEFFICIENT_REPLAY_WITH_TAG_PAIRS",
        "date":"2026-10-28",
        "analysis_stage":"POST_NEW_SOURCE_EFFECTS_DESCRIPTIVE_EFFECT_REPLAY_NOT_OLD_PHYSICAL_ARRAY_REPLAY",
        "original_field_studies":2,
        "original_source_taxon_groups":6,
        "original_source_weather_contexts":12,
        "native_original_context_and_three_axis_comparisons":len(deltas_context),
        "max_abs_new_vs_old_original_context_fitted_coefficient_difference":max(deltas_context),
        "native_original_source_animal_mean_comparisons":len(deltas_source),
        "max_abs_new_vs_old_original_source_tag_means_difference":max(deltas_source),
        "all_original_36_context_fitted_coefficients_reproduced_to_1e-12":
            len(deltas_context)==36 and max(deltas_context)<=1e-12,
        "all_original_18_source_animal_mean_coefficients_reproduced_to_1e-12":
            len(deltas_source)==18 and max(deltas_source)<=1e-12,
        "paired_original_same_tag_context_nominal_descriptive_bootstrap_sign_counts":source_paired,
        "all_six_source_taxon_comparisons":rows,
        "current_provider_numeric_physical_sha_matched_contexts":12,
        "historic_12_raw_JSON_SHA_source_HOLDS_UNCHANGED":True,
        "old_v0_weather_physical_array_identity":"NOT_IDENTIFIABLE_FROM_UNARCHIVED_HISTORIC_HOURLY_ARRAYS",
        "historical_external_negative_radiation_nonreplication_UNCHANGED":NEG,
        "interpretation_ceiling":[
            "36 original context and 18 tagged-source fitted effect values matching numerically proves repeatability of those saved statistical summaries under the new snapshot, NOT byte identity or equality of all unarchived old hourly weather arrays.",
            "This is a same-source same-sample post-result reanalysis, NOT external biological field study replication, true climate model version equivalence or causal solar navigation.",
            "Only q-coverage-ranked context2 minus context1 within same tagged animals is contrasted; not higher/lower irradiance or common ecological gradient.",
            "Paired animal nominal bootstrap 95% intervals are exploratory and NOT corrected for multiple testing. Zero-spanning does not prove no effect.",
            "The old standalone strict raw-JSON equality ERA5_VERSION_HOLD from the earlier V1 route is retained as a historic source-identity outcome, even when biological fitted coefficients are numerically reproduced.",
            "Utah/Wyoming source field studies are nonrandom and field/tag/season/route-dependent, so this is not a world/global hierarchical ecology estimate."
        ],
        "original_blob_SHA_inputs":FILES
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,default=OUTPUT)
    a=cli.parse_args()
    d=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                     encoding="utf-8")
    print("RSDS_CANONICAL_EFFECT_REPLAY",json.dumps({
        "sources":d["original_source_taxon_groups"],
        "context_comparisons":d["native_original_context_and_three_axis_comparisons"],
        "max_context_abs_diff":d["max_abs_new_vs_old_original_context_fitted_coefficient_difference"],
        "source_means":d["native_original_source_animal_mean_comparisons"],
        "max_source_abs_diff":d["max_abs_new_vs_old_original_source_tag_means_difference"],
        "paired":d["paired_original_same_tag_context_nominal_descriptive_bootstrap_sign_counts"],
        "old_version_hold":d["historic_12_raw_JSON_SHA_source_HOLDS_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
