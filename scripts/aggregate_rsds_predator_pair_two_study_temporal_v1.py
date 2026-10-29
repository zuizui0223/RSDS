"""Source-complete two-study coyote/bobcat diel q and photoperiod atlas.

Includes BOTH pre-registered source studies and EVERY native source taxon,
including the Illinois white-tailed deer. A source first solar-azimuth study
being already viewed is disclosed; temporal q was not selected by effect sign.
No source replacement on MD5/collar schedule/missing night-day support.

Interpretation is source-conditional descriptive, not two independent
species-level global parameters, causal seasonal light or predation behavior.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"research/RSDS_COYOTE_BOBCAT_TWO_STUDY_TIME_NICHE_V1_CONTRACT.json"
OUTPUT=ROOT/"research/RSDS_COYOTE_BOBCAT_TWO_STUDY_SHARED_DATE_TEMPORAL_V1_FIRST_RESULT_20261029.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
PRIMARY=("Canis latrans","Lynx rufus")


def check_axis(x,kind):
    if x is None:
        return {"status":"UNAVAILABLE","mean":None,"ci":None,"source_tags":0}
    if not isinstance(x,dict):
        raise ValueError("Original native source first temporal effect malformed")
    n=x["original_tagged_animals"] if kind=="shared" else x["n_individuals"]
    mean=x["source_animal_mean"] if kind=="shared" else x["mean"]
    ci=x["source_tag_bootstrap95"] if kind=="shared" else x["animal_bootstrap_95_percentile"]
    sign=x["nominal_sign"] if kind=="shared" else None
    if n<5:
        if mean is not None or ci is not None:
            raise ValueError("Original native source effect has fake mean on inadequate tag support")
        return {"status":"UNAVAILABLE","mean":None,"ci":None,"source_tags":n}
    if not isinstance(ci,list) or len(ci)!=2 or mean is None or ci[0]>mean or mean>ci[1]:
        raise ValueError("Source effect missing valid original animal bootstrap")
    now="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
    if sign is not None and sign!=now:
        raise ValueError("Native original source CI sign was rewritten")
    return {"status":now,"mean":mean,"ci":ci,"source_tags":n}


def aggregate(sources):
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (spec["version"]!="RSDS_COYOTE_BOBCAT_TWO_STUDY_TEMPORAL_NICHE_V1"
            or set(sources)!={"Washington","Illinois"}
            or len(spec["source_studies"])!=2):
        raise ValueError("Original fixed 2-source coyote/bobcat study census incomplete")
    full=[]
    status={}
    for src in spec["source_studies"]:
        name=src["region"]
        row=sources[name]
        if (row["source_region"]!=name or row["original_study_DOI"]!=src["doi"] or
            row["native_taxa"]!=src["native_taxa"] or
            row["version"]!="RSDS_COYOTE_BOBCAT_SHARED_DATE_TEMPORAL_ORIGINAL_STUDY_V1" or
            row["old_12_ERA5_RAW_SHA_HOLDS_UNCHANGED"] is not True or
            row["old_48_source_sun_direction_external_nonreplication_UNCHANGED"]!=NEGATIVE):
            raise ValueError("An original source DOI/native taxon was changed")
        if row["status"]=="ORIGINAL_SOURCE_SHARED_CALENDAR_EVALUATED" and (
            row["provider_MD5_verified"] is not True or
            row["first_solar_source_native_hour_fix_counts_replayed"] is not True
        ):
            raise ValueError("Source native temporal effect computed without MD5 and prior source replay")
        if {x["taxon"] for x in row["original_native_source_full_temporal_q_and_daylength"]}!=set(src["native_taxa"]):
            raise ValueError("Original native source taxon (including Illinois deer) omitted")
        native={}
        for x in row["original_native_source_full_temporal_q_and_daylength"]:
            t=x["taxon"]
            native[t]={
                "original_hourbin_source_tag_count":x.get("original_tag_count"),
                "original_eligible_q_tag_count":x.get("source_q_original_tag_support"),
                "q":check_axis(x.get("q"),"source"),
                "q_daylength":check_axis(x.get("q_daylength"),"source")
            }
        shared={}
        for stage,label in [("shared_original_calendar_comparison","ORIGINAL_ONLY"),
                            ("shared_triple_GPS_clock_comparison","THREE_SCHEDULES")]:
            x=row.get(stage)
            if x is None:
                shared[label]={"status":"UNAVAILABLE","source_calendar_days":0,"native_source_pair":None}
                continue
            qs=x["native_pair"]["coyote_minus_bobcat_q"]
            gs=x["native_pair"]["coyote_minus_bobcat_daylength"]
            sch=["original"] if label=="ORIGINAL_ONLY" else ["original","even","odd"]
            if set(qs)!=set(sch) or set(gs)!=set(sch):
                raise ValueError("Source predeclared clock-schedule native pair missing")
            contrasts={}
            for s in sch:
                a=qs[s]
                b=gs[s]
                contrasts[s]={
                    "q_coyote_minus_bobcat":a,
                    "daylength_coyote_minus_bobcat":b,
                    "nominal_q_direction":a["nominal_sign"],
                    "nominal_daylength_direction":b["nominal_sign"]
                }
            shared[label]={
                "status":x["status"],
                "source_calendar_days":x["source_calendar_dates_with_both_taxa"],
                "source_year_month_blocks":x["source_month_year_blocks"],
                "native_original_tags_on_shared_source_dates":x["source_native_taxon_tag_support_on_shared_calendar"],
                "source_native_pair":contrasts,
                "calendar_year_month_block_pair":x["date_level_pair_difference_vs_daylength"],
            }
        full.append({"source_region":name,"source_original_study_DOI":src["doi"],
                     "original_source_status":row["status"],
                     "prior_original_solar_analysis_was_opened":True,
                     "source_provider_MD5_replayed":row["provider_MD5_verified"],
                     "source_native_taxa":native,
                     "shared_calendar":shared,
                     "no_tag_ID_or_GPS_event_exported":True})
        status[name]=row["status"]
    w,i=full
    recur={}
    for part in ("ORIGINAL_ONLY","THREE_SCHEDULES"):
        phases=["original"] if part=="ORIGINAL_ONLY" else ["original","even","odd"]
        v={}
        for phase in phases:
            data={}
            for s in (w,i):
                case=s["shared_calendar"][part]
                z=case["native_source_pair"]
                data[s["source_region"]]={
                   axis:{"nominal_sign":z[phase]["nominal_"+axis+"_direction"],
                         "effect":z[phase][f"{axis}_coyote_minus_bobcat"]["coyote_minus_bobcat_mean"]}
                        if z is not None else {"nominal_sign":"UNAVAILABLE","effect":None}
                   for axis in ("q","daylength")
                }
            matched={}
            for axis in ("q","daylength"):
                signs={data[k][axis]["nominal_sign"] for k in data}
                matched[axis]=(
                    "SAME_NOMINAL_SIGN_TWO_SOURCE_STUDIES" if len(signs)==1 and
                    signs.pop() in ("POSITIVE","NEGATIVE") else
                    "OPPOSITE_NOMINAL_SIGNS_TWO_STUDIES" if signs=={"POSITIVE","NEGATIVE"}
                    else "INCONCLUSIVE_OR_SOURCE_HOLD")
            v[phase]={"source_study_native_pair":data,"two_source_same_sign_classification":matched}
        recur[part]=v
    return {
        "version":"RSDS_TWO_ORIGINAL_SOURCE_CARNIVORE_SHARED_DAY_TEMPORAL_V1_FIRST_RESULT",
        "date":"2026-10-29",
        "original_study_DOI_count":2,
        "native_study_taxon_rows":5,
        "preselected_primary_native_pair":list(PRIMARY),
        "source_original_outcome_states":status,
        "all_original_source_study_native_taxa":full,
        "two_study_native_pair_repeatability":recur,
        "prior_GPS_solar_azimuth_exploratory_source_opened":True,
        "new_temporal_effects_evaluated_after_specific_source_selection":True,
        "previous_12_ERA5_source_JSON_HOLDS_UNCHANGED":True,
        "previous_48_source_negative_direction_shortwave_external_replication_UNCHANGED":NEGATIVE,
        "scientific_ceilings":[
            "Two published source DOI studies were already exposed for solar azimuth analysis in 2026-10-25. This is a new temporally defined original-GPS analysis, not a blinded new-source empirical replication.",
            "The primary pair is sympatric coyote and bobcat recorded in Illinois versus northern Washington. The Illinois white-tailed deer remains a native source record, not a substitute for unsupported predator tags.",
            "A source-level shared calendar q comparison requires at least 5 tagged animals per native group and 10 valid dates/year; insufficient source calendar support is UNAVAILABLE not zero.",
            "Original one-fix-per-UTC-hourbin does not prove true hourly data capture, and UTC-even/odd schedules can be unsupported or biased.",
            "q is NIGHT minus DAY observed straight-line GPS speed, not latent biological hours of activity, predation risk or foraging time.",
            "An astronomical 43N standardized seasonal daylength is a common mathematical reference, not experienced local photoperiod or causal solar exposure.",
            "Study-local species contrasts in two published field studies cannot identify globally conserved guild niche, solar-thermal mechanism, species adaptation or direct predator-prey interactions.",
            "All source native/collar/calendar missingness and previous radiation nonreplication/ERA5 weather holds remain unchanged."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,default=OUTPUT)
    a=p.parse_args()
    files={}
    for s in ("Washington","Illinois"):
        file=a.input/("source_"+s.lower()+"_predator_time.json")
        if not file.is_file():
            raise ValueError("Original pinned source DOI missing, cannot substitute study "+s)
        files[s]=json.loads(file.read_text(encoding="utf-8"))
    d=aggregate(files)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_COYOTE_BOBCAT_TWO_STUDY_TIME",json.dumps({
        "studies":d["original_study_DOI_count"],
        "native_taxa":d["native_study_taxon_rows"],
        "status":d["source_original_outcome_states"],
        "original_q_contrast":d["two_study_native_pair_repeatability"]["ORIGINAL_ONLY"]["original"]["two_source_same_sign_classification"],
        "old_direction":d["previous_48_source_negative_direction_shortwave_external_replication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
