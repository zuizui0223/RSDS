"""Exact SAME original animal/year/astronomical date across THREE GPS clock phases.

Post-result audit of originally observed two parity-specific measurement biases:
a previous comparison matched ORIGINAL separately to EVEN and ODD 2h source
samples, but the two phase comparisons could have different original animals
and dates. This eliminates that source selection difference by intersecting
the exact qualifying DATE set for all three original fixed schedules.

No weather is retrieved or new biological selection made. No individual ID,
GPS position or observation timestamp leaves private CI runner.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from statistics import mean

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import slope
from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    verify_frame, source_reference, source_baseline_check,
    per_tag_year_days, PAIRS, NEGATIVE, git_blob
)

CONTRACT=HERE/"research/RSDS_GPS_TRIPLE_MATCHED_UTC_PARITY_V2_CONTRACT.json"
FIRST=HERE/"research/RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_FIRST_RESULT_20261028.json"
ALIAS=HERE/"research/RSDS_EIGHT_STUDY_GPS_CADENCE_PHASE_ALIAS_DIAGNOSTIC_20261028.json"
SEED=20261028
NBOOT=2000
AXES=("q_even_minus_original", "q_odd_minus_original", "q_even_minus_odd",
      "daylength_even_minus_original","daylength_odd_minus_original",
      "daylength_even_minus_odd")


def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c["version"]!="RSDS_EXACT_TRIPLE_TAG_DATE_UTC_PHASE_GPS_V2":
        raise ValueError("Original triple-match source inferential contract changed")
    def pin(p,expected):
        payload=p.read_bytes()
        if git_blob(payload)!=expected:
            raise ValueError("Prior ORIGINAL source phase result changed")
        return json.loads(payload.decode("utf-8"))
    first=pin(FIRST,c["first_result_git_blob_SHA"])
    diag=pin(ALIAS,c["phase_alias_git_blob_SHA"])
    if (first["original_field_studies"]!=8 or first["original_source_taxon_rows"]!=12
            or diag["original_field_studies"]!=8 or diag["source_taxon_rows"]!=12
            or diag["source_cohorts_with_native_UTC_phase_fix_schedule_imbalance_below_20pct"]!=4
            or first["historical_48_source_negative_radiation_replication_UNCHANGED"]!=NEGATIVE
            or not first["historic_ERA5_old_12_source_raw_SHA_HOLDS_UNCHANGED"]):
        raise ValueError("Original eight source study census, native taxon and holds changed")
    return c,first


def triple_from_private_original_source(full,even,odd):
    """Full, even and odd are original PRIVATE tag->year->date:(q,L) tables.

    Returns only private per-tag deltas and aggregate support. Each q delta is
    formed on the same triple-intersection dates/year, then equally weighted
    across eligible original animal calendar years.
    """
    keys=sorted(set(full)&set(even)&set(odd))
    effects={}
    summary={
        "original_source_tags_with_at_least_one_day_in_all_three_schedules":len(keys),
        "triply_qualified_tag_years_with_at_least_ten_identical_dates":0,
        "triply_qualified_original_animal_days":0,
        "source_original_tags_with_at_least_ten_triple_days":0,
        "source_original_tags_with_daylength_slope_triple_support":0,
    }
    for tag in keys:
        per_axis={name:[] for name in AXES}
        for year in sorted(set(full[tag])&set(even[tag])&set(odd[tag])):
            shared=sorted(set(full[tag][year])&set(even[tag][year])&set(odd[tag][year]))
            if len(shared)<10:
                continue
            summary["triply_qualified_tag_years_with_at_least_ten_identical_dates"]+=1
            summary["triply_qualified_original_animal_days"]+=len(shared)
            x=[full[tag][year][date][1] for date in shared]
            q0=[full[tag][year][date][0] for date in shared]
            qe=[even[tag][year][date][0] for date in shared]
            qo=[odd[tag][year][date][0] for date in shared]
            per_axis["q_even_minus_original"].append(mean(qe)-mean(q0))
            per_axis["q_odd_minus_original"].append(mean(qo)-mean(q0))
            per_axis["q_even_minus_odd"].append(mean(qe)-mean(qo))
            grads=[slope(list(zip(x,y))) for y in (q0,qe,qo)]
            if all(v is not None for v in grads):
                per_axis["daylength_even_minus_original"].append(grads[1]-grads[0])
                per_axis["daylength_odd_minus_original"].append(grads[2]-grads[0])
                per_axis["daylength_even_minus_odd"].append(grads[1]-grads[2])
        if per_axis["q_even_minus_original"]:
            effects[tag]={key:mean(vals) for key,vals in per_axis.items() if vals}
    summary["source_original_tags_with_at_least_ten_triple_days"]=len(effects)
    summary["source_original_tags_with_daylength_slope_triple_support"]=sum(
        "daylength_even_minus_original" in x for x in effects.values()
    )
    return effects,summary


def bootstrap(values,seed):
    if len(values)<5:return None
    rng=random.Random(seed)
    n=len(values)
    ordered=sorted(mean(values[rng.randrange(n)] for _ in range(n))
                   for _ in range(NBOOT))
    return [ordered[49],ordered[1949]]


def aggregate_private_triple(effects):
    """Do not leak PRIVATE keys or individual coefficients."""
    result={}
    for k,axis in enumerate(AXES):
        vals=[d[axis] for d in effects.values() if axis in d]
        n=len(vals)
        ci=bootstrap(vals,SEED+23*k)
        if ci is None:
            sign="UNAVAILABLE"
            m=None
        else:
            m=mean(vals)
            sign=("NOMINAL_POSITIVE" if ci[0]>0 else
                  "NOMINAL_NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO")
        result[axis]={
            "source_original_same_tags_and_triple_dates":n,
            "group_source_mean_delta":m,
            "descriptive_source_tag_bootstrap95":ci,
            "status":sign,
            "positive_source_tag_changes":sum(v>0 for v in vals),
            "negative_source_tag_changes":sum(v<0 for v in vals),
        }
    qa=result["q_even_minus_original"]["group_source_mean_delta"]
    qb=result["q_odd_minus_original"]["group_source_mean_delta"]
    qc=result["q_even_minus_odd"]["group_source_mean_delta"]
    if qa is not None and abs(qa-qb-qc)>1e-12:
        raise ValueError("Exact same-date original q parity mean difference not algebraically consistent")
    a=result["daylength_even_minus_original"]["group_source_mean_delta"]
    b=result["daylength_odd_minus_original"]["group_source_mean_delta"]
    c=result["daylength_even_minus_odd"]["group_source_mean_delta"]
    if a is not None and abs(a-b-c)>1e-12:
        raise ValueError("Exact same-date q/daylength slope parity decomposition not consistent")
    return result


def sample_native_taxon(taxon,tracks,archived):
    reference=next(x for x in archived["source_native_taxon_results"] if x["taxon"]==taxon)
    first=reference.get("even_odd_cadence_effects")
    if not isinstance(first,dict) or set(first)!=set(PAIRS):
        raise ValueError("Original source native taxon phase first result missing")
    full,n0=per_tag_year_days(tracks,"ORIGINAL_1H")
    even,ne=per_tag_year_days(tracks,PAIRS[0])
    odd,no=per_tag_year_days(tracks,PAIRS[1])
    for n in (n0,ne,no):
        if n["original_selected_1hour_source_fixes"]!=reference[
                "original_hourly_GPS_data_support"]["original_selected_1hour_source_fixes"]:
            raise ValueError("Original source first GPS bin count changed")
    if n0["tags_with_some_q_observable"]!=reference[
            "original_hourly_GPS_data_support"]["tags_with_some_q_observable"]:
        raise ValueError("Original native source q qualified-tag coverage changed")
    for phase,now in zip(PAIRS,(ne,no)):
        prev=first[phase]["data_support"]
        for key in ("UTC_phase_retained_GPS_fixes","qualified_same_solar_period_GPS_intervals",
                    "animal_calendar_days_with_ge3_day_and_ge3_night_segments"):
            if now[key]!=prev[key]:
                raise ValueError("Original source phase GPS interval/day support changed: "+key)
    paired,denom=triple_from_private_original_source(full,even,odd)
    max_even=first[PAIRS[0]]["matched_q_velocity_proxy_difference"][
        "n_original_tagged_animals_paired_between_cadences"]
    max_odd=first[PAIRS[1]]["matched_q_velocity_proxy_difference"][
        "n_original_tagged_animals_paired_between_cadences"]
    if denom["source_original_tags_with_at_least_ten_triple_days"]>min(max_even,max_odd):
        raise ValueError("Triple original same-tag support exceeds previously frozen phase-specific support")
    effect=aggregate_private_triple(paired)
    return {
        "taxon":taxon,
        "original_same_tag_date_triple_source_support":denom,
        "original_separate_phase_q_paired_tag_upper_bounds":{
            "even":max_even,"odd":max_odd
        },
        "triple_source_effects":effect,
        "source_status":"TRIPLE_ORIGINAL_TAG_UTC_PARITY_ESTIMABLE" if len(paired)>=5
                        else "INSUFFICIENT_TRIPLE_TAG_DATE_SUPPORT",
        "not_sunlight_or_biological_phase_causality":True
    }


def run_source(index,out):
    c,prior=frozen()
    old,new,fixed,meta=verify_frame()
    if index not in range(8):
        raise ValueError("Original 8-study roster cannot be changed")
    context=("Utah" if index==0 else "Wyoming" if index==1 else str(index))
    src,taxa,ref,provider=source_reference(context,old,new)
    first=prior["all_original_source_study_and_native_taxa"][index]
    if (first["fixed_original_study_index"]!=index or
            first["original_study_doi"]!=src["doi"] or
            {r["taxon"] for r in first["source_native_taxon_results"]}!=set(taxa)):
        raise ValueError("Original cadence first study and taxon identities not frozen")
    result={
        "version":"RSDS_EXACT_TRIPLE_MATCHED_ORIGINAL_GPS_STUDY_V2",
        "date":"2026-10-28","original_study_index":index,
        "original_study_doi":src["doi"],"original_native_taxa":list(taxa),
        "source_status":"ORIGINAL_SOURCE_HOLD",
        "original_source_provider_checksum_verified":False,
        "native_q_daylength_baseline_replayed_before_new_phase_comparison":False,
        "native_taxon_results":[],
        "historical_12_ERA5_JSON_SOURCE_HOLDS_UNCHANGED":True,
        "historical_signed_shortwave_external_nonreplication_UNCHANGED":NEGATIVE
    }
    if first["source_status"]=="PREDECLARED_FILE_METADATA_HOLD":
        result["source_status"]="PREDECLARED_FILE_METADATA_HOLD"
        result["source_hold_reason"]="ORIGINAL_SOURCE_FILE_OVER_PREDECLARED_SIZE_CAP"
        result["native_taxon_results"]=[{
            "taxon":taxon,"source_status":"SOURCE_HOLD",
            "triple_source_effects":None
        } for taxon in taxa]
    else:
        is_old=index<2
        pkg=src["package_uuid"] if is_old else src["id"]
        ident=src["original_bitstream_uuid"] if is_old else provider["original_file"]["bitstream_uuid"]
        size=src["original_bytes"] if is_old else provider["original_file"]["bytes"]
        md5=src["original_md5"] if is_old else provider["original_file"]["checksum"]["value"]
        original=provider["original_file"]["filename"].lower() if provider else ""
        suffix=".zip" if index==0 or original.endswith(".zip") else ".csv.gz" if original.endswith(".csv.gz") else ".csv"
        file=out/("private_original_triple_"+str(index)+suffix)
        try:
            options=[f for f in original_files(pkg) if f["bitstream_uuid"]==ident]
            if (len(options)!=1 or options[0]["bytes"]!=size or
                    options[0]["checksum"]["value"].lower()!=md5.lower()):
                raise ValueError("Frozen original source provider UUID/file bytes/checksum changed")
            verify=checksum_and_download(options[0],file)
            if not verify["provider_checksum_verified"]:
                raise ValueError("Frozen original original-file provider MD5 not reproduced")
            result["original_source_provider_checksum_verified"]=True
            if index==0:
                animals,counts,other=extract_time_only_original_source(file,taxa)
            elif index==1:
                animals,counts,other=streaming_time_selected_csv(file,taxa)
            else:
                animals,counts,other=load_source_taxa(file,provider["original_file"]["filename"],taxa)
            source_baseline_check(animals,ref)
            result["native_q_daylength_baseline_replayed_before_new_phase_comparison"]=True
            print("RSDS_TRIPLE_UTC_NATIVE_ORIGINAL_Q_REPLAY",index,len(taxa),flush=True)
            for taxon in taxa:
                result["native_taxon_results"].append(
                    sample_native_taxon(taxon,animals[taxon],first)
                )
                animals[taxon].clear()
            result["source_status"]="TRIPLE_SOURCE_SUPPORT_EVALUATED"
        except Exception as e:
            result["source_status"]="ORIGINAL_SOURCE_OR_TRIPLE_MATCH_HOLD"
            result["source_error_type"]=type(e).__name__
        finally:
            file.unlink(missing_ok=True)
        if not result["native_taxon_results"]:
            result["native_taxon_results"]=[{
                "taxon":taxon,"source_status":"SOURCE_HOLD",
                "triple_source_effects":None
            } for taxon in taxa]
    out.mkdir(parents=True,exist_ok=True)
    filename=out/("original_source_"+str(index)+"_triple_UTC_v2.json")
    filename.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                        encoding="utf-8")
    print("RSDS_TRIPLE_ORIGINAL_STUDY_RESULT",json.dumps({
        "study_index":index,"status":result["source_status"],
        "native_taxa":len(result["native_taxon_results"]),
        "q_replayed":result["native_q_daylength_baseline_replayed_before_new_phase_comparison"],
        "native_taxa_with_5_triple_tags":sum(
            x.get("source_status")=="TRIPLE_ORIGINAL_TAG_UTC_PARITY_ESTIMABLE"
            for x in result["native_taxon_results"]
        )
    },sort_keys=True),flush=True)
    return result


def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--source-index",type=int,required=True)
    arg.add_argument("--out",type=Path,required=True)
    a=arg.parse_args()
    if a.source_index not in range(8):
        arg.error("Exactly previously frozen eight original source studies, indices 0..7")
    a.out.mkdir(parents=True,exist_ok=True)
    run_source(a.source_index,a.out)


if __name__=="__main__":
    main()
