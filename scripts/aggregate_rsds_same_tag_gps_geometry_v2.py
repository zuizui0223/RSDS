"""Aggregate new weather-independent same-tagged-animal GPS time+direction results."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from scripts.run_rsds_same_tag_gps_geometry_v2 import CONTRACT, ERA5_FIRST_HOLD, INITIAL_NEGATIVE

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_FIRST_RESULT_20261026.json"


def aggregate(sources):
    contract=json.loads(CONTRACT.read_text(encoding="utf-8"))
    old=json.loads(ERA5_FIRST_HOLD.read_text(encoding="utf-8"))
    if old["heldout_negative_radiation_nonreplication_UNCHANGED"] != INITIAL_NEGATIVE:
        raise ValueError("Original radiation external decision has changed")
    expected={v["region"]:v for v in contract["source_original_frames"]}
    if set(sources)!=set(expected) or len(expected)!=2:
        raise ValueError("Not two pre-fixed original GPS source studies")
    rows=[]
    states={}
    for region in ("Utah","Wyoming"):
        s=sources[region]
        if (s.get("region")!=region or s.get("version")!=contract["version"]
                or s.get("source_doi")!=expected[region]["doi"]):
            raise ValueError("Original source identity does not match contract")
        taxa=s["source_taxa"]
        if {r["taxon"] for r in taxa}!=set(expected[region]["taxa"]):
            raise ValueError("Incomplete original native taxon frame")
        states[region]=s["status"]
        for x in taxa:
            paired=x.get("paired_context_effect")
            if paired is not None:
                n=paired["same_original_tagged_animals_all_axes"]
                if n < 0:
                    raise ValueError("Negative paired original animal count")
                for axis,item in paired["source_effects"].items():
                    if item["common_original_tagged_animals"] != n:
                        raise ValueError("Not same original tagged animals across axes")
                    if item["status"] == "PAIRED_TAG_ESTIMABLE" and n < 5:
                        raise ValueError("Source effect overclaims unsupported original tags")
                    if item["status"] == "PAIRED_TAG_UNAVAILABLE" and item["rank2_minus_rank1_tag_mean"] is not None:
                        raise ValueError("Undefined source response was recoded numeric")
            rows.append({
                "region":region,"source_doi":s["source_doi"],"taxon":x["taxon"],
                "source_original_q_daylight_exact_replay":x.get("original_source_replay_success",False),
                "new_context_status":x.get("geometry_source_status","SOURCE_HOLD"),
                "all_contexts_selected_by_original_q_coverage":x.get("selection_is_coverage_only",False),
                "original_context_coverage":{
                    k:x.get(k) for k in (
                        "all_original_q_gridyears","eligible_q_coverage_gridyears",
                        "eligible_context_pairs","original_q_tags_in_both_selected_contexts",
                        "selected_context_qualified_q_days",
                        "selected_context_qualified_q_animals")
                },
                "source_geometry_paired":paired
            })
    if len(rows)!=6:
        raise ValueError("Unrecognized taxon/source denominator")
    repeated=[]
    for taxon in ("Odocoileus hemionus","Ovis canadensis"):
        pair=[r for r in rows if r["taxon"]==taxon]
        if len(pair)!=2 or {r["region"] for r in pair}!={"Utah","Wyoming"}:
            raise ValueError("Original two-region common taxon unavailable")
        contrast={}
        for axis in ("q","sun_cos","sun_sin","geographic_marginal","joint_circular"):
            vals=[]
            for r in pair:
                p=r["source_geometry_paired"]
                a=p["source_effects"][axis] if p else None
                vals.append(a["rank2_minus_rank1_tag_mean"] if a and a["status"]=="PAIRED_TAG_ESTIMABLE" else None)
            contrast[axis]={
                "Utah":vals[0] if pair[0]["region"]=="Utah" else vals[1],
                "Wyoming":vals[1] if pair[1]["region"]=="Wyoming" else vals[0],
                "between_study_difference":(vals[0]-vals[1]) if None not in vals else None,
                "status":"TWO_SOURCE_DESCRIPTIVE" if None not in vals else "NOT_ESTIMABLE_TWO_SOURCES"
            }
        repeated.append({"taxon":taxon,"original_studies":2,"axes":contrast})
    return {
        "version":"RSDS_SAME_TAG_GPS_GEOMETRY_V2_FIRST_RESULT",
        "date":"2026-10-26",
        "original_source_count":2,
        "source_taxon_rows":len(rows),
        "original_source_states":states,
        "source_taxa_with_paired_tag_estimates":sum(r["new_context_status"]=="PAIRED_TAG_ESTIMABLE" for r in rows),
        "all_original_source_taxon_results":rows,
        "two_region_same_native_taxon_repeatability":repeated,
        "no_ERA5_retrieved_or_reinterpreted":True,
        "previous_raw_ERA5_SHA_HOLD_unchanged":True,
        "external_48_original_source_shortwave_nonreplication_UNCHANGED":INITIAL_NEGATIVE,
        "nonclaims":[
            "Original q-coverage selected gridyear ranks are not treatment groups sorted by irradiance, latitude, solar season or movement sign.",
            "Same tagged animal contrasts are observational differences between selected recorded GPS dates/sites, not causal Sun guidance or temporal niche plasticity.",
            "Raw Sun-relative movement cosine may change purely from fixed compass bearing and different Sun positions. Geographic marginal and joint circular components are reported separately.",
            "Annual/source q from tracked GPS speed is not true latent activity duration. Individual tags from a field study do not represent independent global ecosystems.",
            "Old ERA5 JSON version holds and old external radiation nonreplication remain immutable. No radiation response inferred in this GPS-only route."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input",type=Path,required=True)
    cli.add_argument("--out",type=Path,default=OUTPUT)
    args=cli.parse_args()
    sources={}
    for region in ("Utah","Wyoming"):
        p=args.input/("source_"+region.lower()+"_geometry_v2.json")
        if not p.is_file():
            raise ValueError("Missing fixed source original GPS result: "+region)
        sources[region]=json.loads(p.read_text(encoding="utf-8"))
    d=aggregate(sources)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_SAME_TAG_GPS_GEOMETRY_V2",json.dumps({
        "original_studies":d["original_source_count"],
        "study_taxa":d["source_taxon_rows"],
        "paired_source_taxa":d["source_taxa_with_paired_tag_estimates"],
        "source_states":d["original_source_states"],
        "old_ERA5_HOLD":d["previous_raw_ERA5_SHA_HOLD_unchanged"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
