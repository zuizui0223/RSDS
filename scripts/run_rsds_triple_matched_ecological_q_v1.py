"""How much of original temporal ecology survives EXACT tag/year/day GPS cadence matching?

Original eight published Movebank studies / twelve native source-taxon rows.
Replays provider original GPS UUID/size/MD5 and archived q/daylength first
source effects before constructing new absolute q and q~astronomical daylength
on identical tagged animal/year/solar dates under ORIGINAL/EVEN/ODD schedules.
The per-tag triply-shared original date set and source delta must EXACTLY
reproduce the independently archived original three-phase measurement study.

All original GPS positions, timestamps, tags and individual metrics are
ephemeral. Only source-level summaries and original support counts published.
No weather is requested; these are sampled DISPLACEMENT, not activity duration.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from pathlib import Path
from statistics import mean
from collections import Counter

from scripts.movebank_native_ecology_20261019 import HERE, original_files, checksum_and_download
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import slope
from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    verify_frame, source_reference, source_baseline_check,
    per_tag_year_days, PAIRS, NEGATIVE, git_blob
)
from scripts.run_rsds_gps_three_phase_exact_dates_v2 import frozen as previous_frozen

CONTRACT=HERE/"research/RSDS_TRIPLE_DATE_ECOLOGICAL_Q_DAYLENGTH_INVARIANCE_V1_CONTRACT.json"
PREVIOUS=HERE/"research/RSDS_GPS_THREE_PHASE_EXACT_TAG_DATE_V2_FIRST_RESULT_20261028.json"
SCHEDULES=("ORIGINAL_UTC_HOURBIN","EVEN_UTC_2H","ODD_UTC_2H")
Q_KEYS=("q_original","q_even","q_odd")
SLOPE_KEYS=("daylength_original","daylength_even","daylength_odd")
NBOOT=2000
SEED=20261028


def freeze():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if spec["version"]!="RSDS_SAME_TAG_TRIPLE_DATE_ECOLOGICAL_Q_AND_DAYLENGTH_INVARIANCE_V1":
        raise ValueError("Original ecological source sign stability scope changed")
    prior_spec,first=previous_frozen()
    raw=PREVIOUS.read_bytes()
    pinned={r["path"]:r["git_blob_sha"] for r in spec["pinned_prior_source_group_results"]}
    if git_blob(raw)!=pinned[str(PREVIOUS.relative_to(HERE))]:
        raise ValueError("Original first triply-matched GPS source group result changed")
    if first["fixed_original_source_studies"]!=8 or first["native_source_study_taxon_rows"]!=12:
        raise ValueError("Original eight study/twelve taxon denominator changed")
    if (first["old_48_original_source_negative_shortwave_direction_external_replication_UNCHANGED"]!=NEGATIVE
            or first["old_12_ERA5_source_RAW_SHA_HOLDS_UNCHANGED"] is not True):
        raise ValueError("Original prior shortwave or weather HOLD reclassified")
    return spec,first


def private_tag_q_and_daylength(full,even,odd):
    """Private per-original-tag means over exactly identical qualified years/days.

    Output is PRIVATE and must NEVER be serialized. For each tag all three
    q estimands share the same years/dates. For q~daylength all three slopes
    also share the same subset of slope-identifiable years.
    """
    q_records={}
    slopes={}
    denominator={
        "original_tags_with_any_qualified_day_in_all_three_schedules":0,
        "original_tags_with_at_least_ten_triple_qualified_days":0,
        "source_tag_years_with_10_or_more_identical_dates":0,
        "matched_original_tagged_animal_calendar_days":0,
        "original_tags_with_common_photoperiod_slope_support":0,
    }
    for tag in sorted(set(full)&set(even)&set(odd)):
        denominator["original_tags_with_any_qualified_day_in_all_three_schedules"]+=1
        tag_q=[]
        tag_s=[]
        for year in sorted(set(full[tag])&set(even[tag])&set(odd[tag])):
            shared=sorted(set(full[tag][year])&set(even[tag][year])&set(odd[tag][year]))
            if len(shared)<10:
                continue
            denominator["source_tag_years_with_10_or_more_identical_dates"]+=1
            denominator["matched_original_tagged_animal_calendar_days"]+=len(shared)
            x=[full[tag][year][d][1] for d in shared]
            response=[
                [full[tag][year][d][0] for d in shared],
                [even[tag][year][d][0] for d in shared],
                [odd[tag][year][d][0] for d in shared],
            ]
            tag_q.append([mean(y) for y in response])
            gradient=[slope(list(zip(x,y))) for y in response]
            if all(v is not None and math.isfinite(v) for v in gradient):
                tag_s.append(gradient)
        if tag_q:
            q_records[tag]={k:mean(v[i] for v in tag_q)
                            for i,k in enumerate(Q_KEYS)}
        if tag_s:
            slopes[tag]={k:mean(v[i] for v in tag_s)
                         for i,k in enumerate(SLOPE_KEYS)}
    denominator["original_tags_with_at_least_ten_triple_qualified_days"]=len(q_records)
    denominator["original_tags_with_common_photoperiod_slope_support"]=len(slopes)
    return q_records,slopes,denominator


def paired_bootstrap(triples, keys, seed=SEED):
    """One resampled tagged-animal index vector reused across all 3 schedule axes."""
    if not triples:
        return {k:{"status":"UNAVAILABLE","n_source_tags":0,"group_original_tag_mean":None,
                   "nominal_source_tag_bootstrap95":None} for k in keys}
    names=sorted(triples)
    n=len(names)
    ci={}
    draws={k:[] for k in keys}
    if n>=5:
        rng=random.Random(seed)
        for i in range(NBOOT):
            sample=[names[rng.randrange(n)] for _ in range(n)]
            for k in keys:
                draws[k].append(mean(triples[tag][k] for tag in sample))
    for k in keys:
        values=[triples[tag][k] for tag in names]
        if n>=5:
            ordered=sorted(draws[k])
            bounds=[ordered[49],ordered[1949]]
            m=mean(values)
            assert bounds[0]-1e-12<=m<=bounds[1]+1e-12
            label=("POSITIVE" if bounds[0]>0 else
                   "NEGATIVE" if bounds[1]<0 else "INCLUDES_ZERO")
        else:
            m=None
            bounds=None
            label="UNAVAILABLE"
        ci[k]={
            "status":label,"n_source_tags":n,
            "group_original_tag_mean":m,
            "nominal_source_tag_bootstrap95":bounds,
            "original_tags_positive":sum(x>0 for x in values),
            "original_tags_negative":sum(x<0 for x in values),
        }
    return ci


def classify_stability(group,keys):
    statuses=[group[k]["status"] for k in keys]
    if all(s=="POSITIVE" for s in statuses): return "THREE_SCHEDULES_NOMINAL_POSITIVE"
    if all(s=="NEGATIVE" for s in statuses): return "THREE_SCHEDULES_NOMINAL_NEGATIVE"
    if any(s=="UNAVAILABLE" for s in statuses): return "UNAVAILABLE"
    if "POSITIVE" in statuses and "NEGATIVE" in statuses:
        return "OPPOSITE_NOMINAL_SOURCE_TAG_SIGN"
    return "UNRESOLVED_AT_LEAST_ONE_ZERO_SPANNING_INTERVAL"


def compare_previous_source_delta(q,grad,prior):
    """Current strict common-tag means must reproduce prior PRIVATE delta result."""
    old=prior["new_strict_triple_original_tag_effects"]
    if old is None or set(old)!={
        "q_even_minus_original","q_odd_minus_original","q_even_minus_odd",
        "daylength_even_minus_original","daylength_odd_minus_original","daylength_even_minus_odd"}:
        raise ValueError("Previous original same-tag/date group effect missing")
    mapped=(
        (q,("q_even","q_original"),"q_even_minus_original"),
        (q,("q_odd","q_original"),"q_odd_minus_original"),
        (q,("q_even","q_odd"),"q_even_minus_odd"),
        (grad,("daylength_even","daylength_original"),"daylength_even_minus_original"),
        (grad,("daylength_odd","daylength_original"),"daylength_odd_minus_original"),
        (grad,("daylength_even","daylength_odd"),"daylength_even_minus_odd"),
    )
    for current,(a,b),name in mapped:
        n=current[a]["n_source_tags"]
        original=old[name]
        if n!=original["source_original_same_tags_and_triple_dates"]:
            raise ValueError("Triply matched original tags denominator differs from first source")
        m=current[a]["group_original_tag_mean"]
        prior_m=original["group_source_mean_delta"]
        if n<5:
            if m is not None or prior_m is not None:
                raise ValueError("Missing original source effect was converted to numeric")
        elif not math.isclose(m-current[b]["group_original_tag_mean"],prior_m,
                              rel_tol=0,abs_tol=1e-12):
            raise ValueError("Absolute common-date q/photoperiod source levels fail to reconstruct original paired delta")


def analyze_taxon(taxon,original,prior):
    full,_=per_tag_year_days(original,"ORIGINAL_1H")
    even,_=per_tag_year_days(original,PAIRS[0])
    odd,_=per_tag_year_days(original,PAIRS[1])
    levels,slopes,den=private_tag_q_and_daylength(full,even,odd)
    previous=prior["triple_original_animal_day_support"]
    if previous is None:
        raise ValueError("An original source-verified native taxon lacks triple source eligibility")
    for now,old in (
        ("original_tags_with_any_qualified_day_in_all_three_schedules","original_source_tags_with_at_least_one_day_in_all_three_schedules"),
        ("original_tags_with_at_least_ten_triple_qualified_days","source_original_tags_with_at_least_ten_triple_days"),
        ("source_tag_years_with_10_or_more_identical_dates","triply_qualified_tag_years_with_at_least_ten_identical_dates"),
        ("matched_original_tagged_animal_calendar_days","triply_qualified_original_animal_days"),
        ("original_tags_with_common_photoperiod_slope_support","source_original_tags_with_daylength_slope_triple_support"),
    ):
        if den[now]!=previous[old]:
            raise ValueError("Original first source year/day/tag matched support changed: "+now)
    q=paired_bootstrap(levels,Q_KEYS,seed=SEED)
    light=paired_bootstrap(slopes,SLOPE_KEYS,seed=SEED+1000)
    compare_previous_source_delta(q,light,prior)
    return {
        "taxon":taxon,"original_three_schedule_source_support":den,
        "q_original_tag_means_on_exact_same_dates":q,
        "q_daylength_per_hour_on_exact_same_dates":light,
        "source_q_sign_recurrence":classify_stability(q,Q_KEYS),
        "source_daylength_slope_sign_recurrence":classify_stability(light,SLOPE_KEYS),
        "methodological_warning":"Three original timestamp schedules measure different displacement chords and observation opportunities, NOT experimental biological times. The original time-bin series does not require 1h actual fix spacing.",
        "source_status":"THREE_SCHEDULE_ECOLOGICAL_Q_ESTIMABLE" if len(levels)>=5 else "INSUFFICIENT_ORIGINAL_TRIPLE_TAG_DATE_SUPPORT"
    }


def run_source(index,out):
    spec,first=freeze()
    old,new,fixed,metadata=verify_frame()
    if index not in range(8):
        raise ValueError("Only eight frozen original study IDs")
    region="Utah" if index==0 else "Wyoming" if index==1 else str(index)
    src,taxa,reference,provider=source_reference(region,old,new)
    previous=first["all_original_source_studies_and_native_tag_triple_q"][index]
    if (previous["original_study_index"]!=index or
            previous["original_study_doi"]!=src["doi"] or
            {r["taxon"] for r in previous["taxa"]}!=set(taxa)):
        raise ValueError("Previously verified native original study/taxon frame changed")
    result={
        "version":"RSDS_8_STUDY_TRIPLE_DATE_ECOLOGICAL_Q_AND_DAYLENGTH_V1_SOURCE",
        "index":index, "original_study_doi":src["doi"],
        "native_taxa":list(taxa),"status":"ORIGINAL_SOURCE_HOLD",
        "original_provider_bitstream_MD5_confirmed":False,
        "original_source_q_and_daylength_replayed_before_new_q_levels":False,
        "native_taxon_results":[],
        "previous_12_ERA5_raw_version_HOLDS_UNCHANGED":True,
        "previous_48_source_signed_shortwave_external_NONREPLICATED_UNCHANGED":NEGATIVE,
    }
    if previous["source_status"]=="PREDECLARED_FILE_METADATA_HOLD":
        result["status"]="PREDECLARED_FILE_METADATA_HOLD"
        result["native_taxon_results"]=[
            {"taxon":t,"source_status":"SOURCE_HOLD","q_original_tag_means_on_exact_same_dates":None,
             "q_daylength_per_hour_on_exact_same_dates":None} for t in taxa
        ]
    else:
        is_old=index<2
        pkg=src["package_uuid"] if is_old else src["id"]
        bit=src["original_bitstream_uuid"] if is_old else provider["original_file"]["bitstream_uuid"]
        nbytes=src["original_bytes"] if is_old else provider["original_file"]["bytes"]
        digest=src["original_md5"] if is_old else provider["original_file"]["checksum"]["value"]
        original_name=provider["original_file"]["filename"].lower() if provider else ""
        ext=".zip" if index==0 or original_name.endswith(".zip") else (
            ".csv.gz" if original_name.endswith(".csv.gz") else ".csv"
        )
        local=out/("private_source_original_q_ecology_"+str(index)+ext)
        try:
            candidates=[x for x in original_files(pkg) if x["bitstream_uuid"]==bit]
            if (len(candidates)!=1 or candidates[0]["bytes"]!=nbytes or
                    candidates[0]["checksum"]["value"].lower()!=digest.lower()):
                raise ValueError("Provider original source UUID/file bytes/checksum changed")
            verified=checksum_and_download(candidates[0],local)
            if not verified["provider_checksum_verified"]:
                raise ValueError("Original source GPS MD5 failed")
            result["original_provider_bitstream_MD5_confirmed"]=True
            if index==0:
                tracks,counts,unknown=extract_time_only_original_source(local,taxa)
            elif index==1:
                tracks,counts,unknown=streaming_time_selected_csv(local,taxa)
            else:
                tracks,counts,unknown=load_source_taxa(local,provider["original_file"]["filename"],taxa)
            source_baseline_check(tracks,reference)
            result["original_source_q_and_daylength_replayed_before_new_q_levels"]=True
            print("RSDS_TRIPLE_Q_BASELINE_SOURCE_REPLAY",index,len(taxa),flush=True)
            for taxon in taxa:
                prior=next(x for x in previous["taxa"] if x["taxon"]==taxon)
                result["native_taxon_results"].append(analyze_taxon(taxon,tracks[taxon],prior))
                tracks[taxon].clear()
            result["status"]="SOURCE_TRIPLE_ECOLOGICAL_Q_COMPLETED"
        except Exception as exc:
            result["status"]="ORIGINAL_PROVIDER_OR_SUPPORT_HOLD"
            result["source_error_type"]=type(exc).__name__
        finally:
            local.unlink(missing_ok=True)
        if not result["native_taxon_results"]:
            result["native_taxon_results"]=[{"taxon":t,"source_status":"SOURCE_HOLD",
                                              "q_original_tag_means_on_exact_same_dates":None,
                                              "q_daylength_per_hour_on_exact_same_dates":None}
                                              for t in taxa]
    out.mkdir(parents=True,exist_ok=True)
    (out/("source_"+str(index)+"_triple_q_ecology.json")).write_text(
        json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8"
    )
    print("RSDS_TRIPLE_Q_ECOLOGICAL_SOURCE",json.dumps({
        "study_index":index,"status":result["status"],
        "native_taxa":len(result["native_taxon_results"]),
        "exact_original_1h_q_daylength_replay":result["original_source_q_and_daylength_replayed_before_new_q_levels"],
        "source_q_supported":sum(t.get("source_status")=="THREE_SCHEDULE_ECOLOGICAL_Q_ESTIMABLE"
                                 for t in result["native_taxon_results"])
    },sort_keys=True),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index",type=int,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    if a.index not in range(8):
        parser.error("Only the eight original study DOI indices 0...7 are valid")
    a.out.mkdir(parents=True,exist_ok=True)
    run_source(a.index,a.out)


if __name__=="__main__":
    main()
