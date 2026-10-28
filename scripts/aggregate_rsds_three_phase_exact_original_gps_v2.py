"""Aggregate eight frozen Movebank sources' SAME-ANIMAL SAME-DATE GPS UTC triple.

New original GPS route intersects original/EVEN/ODD, after discovering that
two separately paired samples can show opposite descriptive bias signs.
One original DOI is one field study, taxa from Utah/Wyoming NOT independent
ecosystems. Every unsupported animal/day condition and historic source HOLD
stays explicit; q is observed displacement not activity time.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean

from scripts.run_rsds_gps_three_phase_exact_dates_v2 import (
    HERE,AXES,NEGATIVE,frozen
)

OUT=HERE/"research/RSDS_GPS_THREE_PHASE_EXACT_TAG_DATE_V2_FIRST_RESULT_20261028.json"


def aggregate(sources):
    _,first=frozen()
    if len(sources)!=8 or set(sources)!=set(range(8)):
        raise ValueError("One of the original eight field studies missing")
    original=first["all_original_source_study_and_native_taxa"]
    studies=[]
    taxa=[]
    support=Counter()
    estimates={axis:Counter() for axis in AXES}
    phase_status={}
    for i,old in enumerate(original):
        n=sources[i]
        if (n.get("original_study_index")!=i or
                n.get("original_study_doi")!=old["original_study_doi"] or
                n.get("version")!="RSDS_EXACT_TRIPLE_MATCHED_ORIGINAL_GPS_STUDY_V2" or
                n.get("historical_12_ERA5_JSON_SOURCE_HOLDS_UNCHANGED") is not True or
                n.get("historical_signed_shortwave_external_nonreplication_UNCHANGED")!=NEGATIVE):
            raise ValueError("Wrong original public source study identity or previous negative verdict")
        required={s["taxon"] for s in old["source_native_taxon_results"]}
        found=n.get("native_taxon_results",[])
        if (set(n.get("original_native_taxa",[]))!=required or
                len(found)!=len(required) or {t["taxon"] for t in found}!=required):
            raise ValueError("Original native taxon denominator incomplete")
        status=n["source_status"]
        if old["source_status"]=="PREDECLARED_FILE_METADATA_HOLD" and status!="PREDECLARED_FILE_METADATA_HOLD":
            raise ValueError("Historical size-cap source HOLD was reclassified")
        if status=="TRIPLE_SOURCE_SUPPORT_EVALUATED" and (
            n.get("original_source_provider_checksum_verified") is not True or
            n.get("native_q_daylength_baseline_replayed_before_new_phase_comparison") is not True):
            raise ValueError("Triple GPS response without provider MD5/original biological replay")
        support[status]+=1
        studyrow={"original_study_index":i,
                  "original_study_doi":old["original_study_doi"],
                  "original_native_taxa":sorted(required),
                  "source_status":status,
                  "source_MD5":n.get("original_source_provider_checksum_verified"),
                  "original_q_replayed":n.get("native_q_daylength_baseline_replayed_before_new_phase_comparison"),
                  "source_hold_error_type":n.get("source_error_type"),
                  "taxa":[]}
        for item in found:
            oldnative=next(x for x in old["source_native_taxon_results"] if x["taxon"]==item["taxon"])
            ctx=item.get("triple_source_effects")
            native={
                "taxon":item["taxon"],"source_status":item.get("source_status","SOURCE_HOLD"),
                "triple_original_animal_day_support":item.get("original_same_tag_date_triple_source_support"),
                "previous_two_separate_parity_pair_support":item.get("original_separate_phase_q_paired_tag_upper_bounds"),
                "new_strict_triple_original_tag_effects":ctx,
            }
            if status=="TRIPLE_SOURCE_SUPPORT_EVALUATED":
                if not isinstance(ctx,dict) or set(ctx)!=set(AXES):
                    raise ValueError("Qualified original native group lacks all six predefined triple axes")
                current=item["original_same_tag_date_triple_source_support"]
                count=current["source_original_tags_with_at_least_ten_triple_days"]
                firstPhase=oldnative["even_odd_cadence_effects"]
                if count>min(
                    firstPhase["EVEN_UTC_HOURS"]["matched_q_velocity_proxy_difference"][
                      "n_original_tagged_animals_paired_between_cadences"],
                    firstPhase["ODD_UTC_HOURS"]["matched_q_velocity_proxy_difference"][
                      "n_original_tagged_animals_paired_between_cadences"]
                ):
                    raise ValueError("Triple original tag sample exceeds first phase sample")
                for axis in AXES:
                    val=ctx[axis]
                    if val["source_original_same_tags_and_triple_dates"]>count:
                        raise ValueError("An axis was estimated from non-triple original tags")
                    st=val["status"]
                    if st not in ("UNAVAILABLE","NOMINAL_POSITIVE","NOMINAL_NEGATIVE","INCLUDES_ZERO"):
                        raise ValueError("Unknown strict triple source response status")
                    if st=="UNAVAILABLE" and (val["source_original_same_tags_and_triple_dates"]>=5 or
                        val["group_source_mean_delta"] is not None or val["descriptive_source_tag_bootstrap95"] is not None):
                        raise ValueError("Unsupported strict same-date original response has a fake nonzero estimate")
                    if st!="UNAVAILABLE":
                        ntag=val["source_original_same_tags_and_triple_dates"]
                        ci=val["descriptive_source_tag_bootstrap95"]
                        m=val["group_source_mean_delta"]
                        if (ntag<5 or not isinstance(ci,list) or len(ci)!=2 or
                           not all(math.isfinite(x) for x in (m,*ci)) or not ci[0]<=m<=ci[1]):
                            raise ValueError("Source original tagged-animal triple effect/uncertainty invalid")
                    estimates[axis][st]+=1
            phase_status[native["source_status"]]=phase_status.get(native["source_status"],0)+1
            taxa.append({"source_study_doi":old["original_study_doi"],**native})
            studyrow["taxa"].append(native)
        studies.append(studyrow)
    if len(studies)!=8 or len(taxa)!=12:
        raise ValueError("Wrong native original eight-study twelve-taxon study cohort")
    return {
        "version":"RSDS_EIGHT_SOURCE_EXACT_ORIGINAL_TAG_DATE_GPS_PARITY_V2_FIRST_RESULT",
        "date":"2026-10-28",
        "research_stage":"POST_FIRST_PARITY_OUTCOMES_STRICT_SHARED_ORIGINAL_TAG_DATE_REANALYSIS",
        "fixed_original_source_studies":8,
        "native_source_study_taxon_rows":12,
        "source_study_states":dict(support),
        "source_native_taxon_states":phase_status,
        "strict_three_phase_axis_nominal_tag_bootstrap_sign_counts":{
            axis:dict(c) for axis,c in estimates.items()
        },
        "all_original_source_studies_and_native_tag_triple_q":studies,
        "historical_original_UTC_fix_phase_alias_first_result_git_blob":
            "670199d6a40f3271f2fffa896c1e9971c66a737c",
        "old_12_ERA5_source_RAW_SHA_HOLDS_UNCHANGED":True,
        "old_48_original_source_negative_shortwave_direction_external_replication_UNCHANGED":NEGATIVE,
        "claim_ceiling":[
            "All 3 q models use the identical original animal, calendar year and original astronomy-qualified dates. Removing different source-calendar support does not turn the observation schedules into ecological treatments.",
            "The original-hourbin GPS series can be native 2-hour duty-cycle and does NOT represent a true 1-hour FIX acquisition at all source studies.",
            "Even UTC-hour versus odd UTC-hour phase is an instrument collection phase, not ecological diurnal phase. Recording bias, GPS missingness, moving interval chord length and true travel path curvature remain.",
            "Under five triple-identical source tags, effects are UNAVAILABLE, not zero or failure to replicate.",
            "Native taxa within original Utah/Wyoming studies are correlated field cohorts; eight original DOI units are purposively sampled not random representative ecosystems.",
            "All original native source MD5/q/daylength first outcomes and prior eight-source phase statistics are protected by sha and replay tests.",
            "Earlier 12 original ERA5 source raw JSON version HOLDS, 48-study external nonreplication and egret original GPS 200 MB HOLD stay unchanged.",
            "No independent solar compass, latent daily time budget, causal photoperiod or shortwave radiation response is established.",
            "Never publish original individual GPS tags, precise locations, UTC event/calendar dates, exact grid years or per-tag scientific estimates."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input",type=Path,required=True)
    cli.add_argument("--out",type=Path,default=OUT)
    a=cli.parse_args()
    full={}
    for i in range(8):
        p=a.input/("original_source_"+str(i)+"_triple_UTC_v2.json")
        if not p.is_file():
            raise ValueError("Missing one of eight frozen original GPS source results: "+str(i))
        full[i]=json.loads(p.read_text(encoding="utf-8"))
    report=aggregate(full)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                     encoding="utf-8")
    print("RSDS_TRIPLE_ORIGINAL_GPS_SOURCE",json.dumps({
        "studies":report["fixed_original_source_studies"],
        "native_taxon":report["native_source_study_taxon_rows"],
        "source_states":report["source_study_states"],
        "tag_states":report["source_native_taxon_states"],
        "q_even_odd":report["strict_three_phase_axis_nominal_tag_bootstrap_sign_counts"]["q_even_minus_odd"],
        "old_failure":report["old_48_original_source_negative_shortwave_direction_external_replication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
