"""Original two-study, six taxon new CANONICAL ERA5 biological source estimates.

A new numeric source version is a post-result reanalysis, NOT replication of
the older unarchived ERA5 physical arrays or the earlier external 48-source
negative radiation finding. Six fixed native taxa, 12 frozen source contexts
and all weak support / version mismatch HOLDs remain in the result.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.run_rsds_canonical_era5_same_tag_ecology_v2 import (
    frozen,NEG,CONTRACT
)
from scripts.movebank_native_ecology_20261019 import HERE

ROOT=HERE
RESULT=HERE/"research/RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_FIRST_RESULT_20261028.json"
AXES=("daypart_q_vs_irradiance","sun_cos_vs_irradiance","sun_sin_vs_irradiance")
SIGN=("NOMINAL_POSITIVE_ANIMAL_BOOTSTRAP","NOMINAL_NEGATIVE_ANIMAL_BOOTSTRAP",
      "INCLUDES_ZERO","UNSUPPORTED")


def classify(axis):
    if axis is None:
        return {"status":"UNSUPPORTED","mean":None,"n_source_animals":0,"animal_bootstrap95":None}
    if not isinstance(axis,dict):
        raise ValueError("Not a source native tagged original biological effect")
    n=axis["n_source_animals"]
    ci=axis["source_animal_bootstrap95"]
    y=axis["source_mean"]
    if n<5:
        if y is not None or ci is not None:
            raise ValueError("Unsupported original source group has a fake climate effect")
        return {"status":"UNSUPPORTED","mean":None,"n_source_animals":n,"animal_bootstrap95":None}
    if y is None or not isinstance(ci,list) or len(ci)!=2 or not ci[0]<=y<=ci[1]:
        raise ValueError("Invalid original animal source bootstrap")
    status=("NOMINAL_POSITIVE_ANIMAL_BOOTSTRAP" if ci[0]>0 else
            "NOMINAL_NEGATIVE_ANIMAL_BOOTSTRAP" if ci[1]<0 else "INCLUDES_ZERO")
    return {"status":status,"mean":y,"n_source_animals":n,"animal_bootstrap95":ci}


def run(utah,wyoming):
    spec,old,physical=frozen()
    sources={"Utah":utah,"Wyoming":wyoming}
    samples={"Utah":spec["study_coverage"]["utah"],
             "Wyoming":spec["study_coverage"]["wyoming"]}
    rows=[]
    statuses=Counter()
    axis_support={k:Counter() for k in AXES}
    paired=Counter()
    ctx=0
    phys_ok=0
    for region in ("Utah","Wyoming"):
        result=sources[region]
        target=samples[region]
        if (result["region"]!=region or
                result["source_study_DOI"]!=target["doi"] or
                result["source_native_taxa"]!=target["native_taxa"] or
                result["version"]!="RSDS_CANONICAL_ERA5_SAME_TAG_ECOLOGY_V2_ORIGINAL_SOURCE_RESULT"
                or not result["previous_ERA5_raw_SHA_12_context_HOLDS_UNCHANGED"]
                or result["prior_48_source_radiation_external_NONREPLICATED_UNCHANGED"]!=NEG
                or not result["new_ERA5_shortwave_effect_is_NOT_replication_of_old_weather_dataset"]):
            raise ValueError("Original source study/outcome scope changed")
        if (result["status"]=="ORIGINAL_SOURCE_REANALYSIS_COMPLETED" and
                (not result["original_GPS_provider_MD5_replayed"]
                 or not result["original_q_and_daylength_replayed_before_weather"])):
            raise ValueError("New weather fitted before source MD5/q/daylength was verified")
        source_taxa=result["taxa"]
        if len(source_taxa)!=len(target["native_taxa"]) or \
           {x["taxon"] for x in source_taxa}!=set(target["native_taxa"]):
            raise ValueError("Original six source taxon frame incomplete or duplicate")
        for t in source_taxa:
            contexts=t.get("context_results",[])
            if contexts and len(contexts)!=2:
                raise ValueError("New current-source weather source context number differs")
            ctx+=len(contexts)
            phys_ok+=sum(x["new_provider_physical_content_status"]==
                         "CANONICAL_PHYSICAL_EXACT_SHA_MATCH" for x in contexts)
            a={name:classify(t.get("source_animal_aggregate",{}).get(name)
                             if isinstance(t.get("source_animal_aggregate"),dict) else None)
               for name in AXES}
            for name,x in a.items():
                axis_support[name][x["status"]]+=1
            status=t["status"]
            statuses[status]+=1
            pair=t.get("within_same_tag_rank2_minus_rank1")
            pstatus="PAIRED_UNAVAILABLE"
            paired_n=0
            if pair is not None:
                paired_n=pair["identical_tag_ids_in_both_original_contexts"]
                if paired_n>=5 and any(x["status"]=="PAIRED_ESTIMABLE"
                                       for x in pair["per_axis"].values()):
                    pstatus="PAIRED_2_CONTEXT_EFFECT_ESTIMABLE"
                elif paired_n<5:
                    pstatus="INSUFFICIENT_PAIRED_TAGS"
            paired[pstatus]+=1
            rows.append({
                "region":region,"source_doi":result["source_study_DOI"],
                "taxon":t["taxon"],
                "new_CANONICAL_weather_state":status,
                "old_historic_ERA5_version_HOLD_remains":True,
                "original_contexts":contexts,
                "source_animal_response_per100W":a,
                "same_original_tagged_animals_in_two_q_coverage_ranked_contexts":paired_n,
                "paired_tagged_animal_contrast_status":pstatus,
                "paired_tagged_context_response_RANK2_MINUS_RANK1":pair,
            })
    if len(rows)!=6:
        raise ValueError("Original two source studies not represented by all six native taxa")
    return {
        "version":"RSDS_CANONICAL_ERA5_SOURCE_ORIGINAL_SAME_TAG_BIOLOGY_V2_FIRST_RESULT",
        "date":"2026-10-28",
        "stage":"POST_OUTCOME_NEW_CANONICAL_PHYSICAL_WEATHER_VERSION_REANALYSIS",
        "original_study_count":2,
        "original_source_taxon_rows":6,
        "selected_original_GPS_gridyear_context_denominator":12,
        "selected_contexts_reported":ctx,
        "current_canonical_weather_contexts_exact_source_SHA_replayed":phys_ok,
        "native_taxon_estimator_states":dict(statuses),
        "source_animal_effect_nominal_sign_support":{k:dict(v) for k,v in axis_support.items()},
        "same_tag_two_context_effect_support":dict(paired),
        "all_fixed_original_study_taxon_results":rows,
        "previous_first_biological_effects_unmodified":True,
        "historical_12_ERA5_raw_JSON_SHA_HOLDs_UNCHANGED":True,
        "old_ERA5_numeric_version_equivalence":"NOT_IDENTIFIABLE_UNARCHIVED_OLD_ARRAYS",
        "historic_48_source_negative_radiation_replication_UNCHANGED":NEG,
        "causal_or_global_claims":[
            "This is a NEW version reanalysis of original two already observed GPS studies; it is NOT an independent heldout source replication.",
            "The frozen current physical hourly ERA5 numeric SHA is checked before effect fitting. Old raw SHA and old physical arrays cannot be retroactively upgraded or compared.",
            "The sign of sun cos or sin response to shortwave anomalies may reflect migration route, habitat, microclimate, timing, wind, temperature and nonrandom tagged individual source selection.",
            "Original source tag/daypart q is a GPS displacement proxy, not latent daily activity time.",
            "Rank2-minus-rank1 context contrast does not assign lower/higher radiation, latitude, phenology, migration or environmental stress.",
            "Source animal bootstrap nominal 95% is not a multiple-testing adjusted world population confidence interval, source ecosystem unit count remains two.",
            "All unavailable weather numeric versions, q/clock support holds, original native taxa and old externally negative radiation replication are retained.",
            "Original animal identifiers, exact GPS coordinates/dates and individual source weather effects not published."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input",type=Path,required=True)
    cli.add_argument("--out",type=Path,default=RESULT)
    a=cli.parse_args()
    d={}
    for reg in ("utah","wyoming"):
        p=a.input/("source_"+reg+"_new_physical_biology_v2.json")
        if not p.is_file():
            raise ValueError("A fixed original source is missing "+reg)
        d[reg]=json.loads(p.read_text(encoding="utf-8"))
    r=run(d["utah"],d["wyoming"])
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,sort_keys=True,indent=2,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    print("RSDS_CANONICAL_ERA5_ORIGINAL_TAG_ECOLOGY",json.dumps({
        "studies":r["original_study_count"],"taxa":r["original_source_taxon_rows"],
        "contexts":r["selected_contexts_reported"],
        "canonical_sha_matched":r["current_canonical_weather_contexts_exact_source_SHA_replayed"],
        "source_signs":r["source_animal_effect_nominal_sign_support"],
        "paired":r["same_tag_two_context_effect_support"],
        "historical_HOLD":r["historical_12_ERA5_raw_JSON_SHA_HOLDs_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
