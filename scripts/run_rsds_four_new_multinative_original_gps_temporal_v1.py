"""Independent multi-native original GPS temporal ecology at shared source days.

Four ORIGINAL public Movebank field studies were fixed by study DOI/source ID/
registered taxa, then exact original primary GPS file name/UUID/bytes/MD5
was sealed BEFORE viewing NEW source animal q outcomes. This program verifies
that seal and reads original animal GPS ONLY in a temporary private CI runner.

For each native taxon it first reports the source-independent individual-level
NIGHT−DAY observed GPS step-speed proxy q and its astronomical daylength
association. It then intersects astronomy-qualified calendar days over ALL
registered taxa in that SAME published study, retaining no favorable subset.
This is co-observation of calendar days, NOT proof of spatial sympatry or
actual animal encounters. All null, cadence, species and study HOLDS remain.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean,median

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download, daylength
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon,slope
from scripts.run_rsds_eight_study_gps_cadence_v1 import per_tag_year_days,git_blob
from scripts.rsds_within_month_photoperiod_shared_calendar_v1 import month_conditional_slope

CONTRACT=HERE/"research/RSDS_FOUR_NEW_MULTINATIVE_SOURCE_GPS_V1_CONTRACT.json"
SOURCE_ONLY=HERE/"research/RSDS_FOUR_NEW_MULTINATIVE_SOURCE_METADATA_FIRST_V1_20261029.json"
PIN=HERE/"research/RSDS_FOUR_NEW_MULTINATIVE_SOURCE_PRE_EVENT_SHA_LOCK_V1.json"
OLD=HERE/"research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
SOURCE_KEY="RSDS_FOUR_NEW_MULTI_NATIVE_GPS_FIRST_ORIGINAL_TEMPORAL_SOURCE_V1"
SEED=20261029
NBOOT=2000


def fixed_sources():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    seal=json.loads(PIN.read_text(encoding="utf-8"))
    raw=SOURCE_ONLY.read_bytes()
    manifest=json.loads(raw.decode("utf-8"))
    old=json.loads(OLD.read_text(encoding="utf-8"))
    if (spec["version"]!="RSDS_FOUR_NEW_MULTINATIVE_ORIGINAL_GPS_STUDY_TEMPORAL_V1"
            or seal["version"]!="RSDS_FOUR_NEW_MULTI_NATIVE_GPS_PRE_EVENT_SOURCE_BLOB_PIN_V1"
            or git_blob(raw)!=seal["original_source_manifest_git_blob_SHA"]
            or manifest["original_published_study_count"]!=4
            or manifest["registered_native_taxon_rows"]!=10
            or manifest["animal_GPS_events_read_by_source_only_preflight"]!=0
            or old["original_study_denominator"]!=8
            or old["old_48_radiation_study_external_negative_direction_test_UNCHANGED"]!=NEG
            or not old["old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED"]):
        raise ValueError("Original source-only pre-event metadata or historical science boundary changed")
    entries=manifest["source_file_entries"]
    assert len(entries)==len(spec["sources"])==4
    for e,s in zip(entries,spec["sources"]):
        if (e["original_repository_package_uuid"]!=s["package_uuid"]
                or e["original_study_id"]!=s["study_id"]
                or e["published_DOI"]!=s["doi"]
                or e["registered_native_taxa"]!=s["registered_taxa"]):
            raise ValueError("Source registry/species or source pre-event metadata modified")
    return spec,entries


def source_reference_latitude(tracks):
    """Study reference lat from native original GPS positions, BEFORE q outcomes.

    Per original tagged animal mean GPS lat; study reference is median of those
    means across ALL registered native taxa. NEVER export this coordinate.
    """
    individuals=[]
    for taxon in tracks:
        for _,fixes in tracks[taxon].items():
            if fixes:
                individuals.append(mean(v[1] for v in fixes.values()))
    if not individuals:
        return None
    center=median(individuals)
    if not -80<=center<=80 or not math.isfinite(center):
        raise ValueError("Source study anatomical solar daylength reference unsupported")
    return center


def common_original_calendar(native):
    """Fixed-point intersection of admissible native tag-years across ALL taxa.

    Selection depends only on source day support, not native q magnitude or
    sign. A tag-year must have ≥10 source q-qualified dates. Study q date
    intersection must be common across every registered native species.
    """
    if len(native)<2:
        raise ValueError("Multispecies original source requires at least two taxa")
    active={
        sp:{tag:{yr:dict(ds) for yr,ds in years.items() if len(ds)>=10}
             for tag,years in original.items() if any(len(ds)>=10 for ds in years.values())}
        for sp,original in native.items()
    }
    rounds=0
    while True:
        rounds+=1
        if rounds>50:raise ValueError("Fixed-point original multispecies calendar did not converge")
        sets=[
            {(yr,d) for years in active[sp].values() for yr,days in years.items() for d in days}
            for sp in native
        ]
        shared=set.intersection(*sets) if sets else set()
        next_group={}
        for sp,tag_data in active.items():
            group={}
            for tag,years in tag_data.items():
                new_years={}
                for yr,daydata in years.items():
                    kept={d:q for d,q in daydata.items() if (yr,d) in shared}
                    if len(kept)>=10:
                        new_years[yr]=kept
                if new_years:
                    group[tag]=new_years
            next_group[sp]=group
        if next_group==active:
            break
        active=next_group
    date_by={
        sp:{(yr,d) for years in active[sp].values() for yr,days in years.items() for d in days}
        for sp in native
    }
    if any(v!=shared for v in date_by.values()):
        raise ValueError("Native source calendar date sets differ across fixed original taxa")
    return active,shared,rounds


def tag_level_q_and_season(group,latitude):
    q_tag={};slope_tag={}
    for tag,years in group.items():
        levels=[];coef=[]
        for yr,days in years.items():
            order=sorted(days)
            if len(order)<10: continue
            v=[days[d][0] for d in order]
            levels.append(mean(v))
            x=[daylength(latitude,d) for d in order]
            if max(x)-min(x)<.5:
                continue
            beta=slope(list(zip(x,v)))
            if beta is not None and math.isfinite(beta):
                coef.append(beta)
        if levels:q_tag[tag]=mean(levels)
        if coef:slope_tag[tag]=mean(coef)
    return q_tag,slope_tag


def bootstrap(values,seed=SEED):
    v=list(values)
    if len(v)<5:return None
    rng=random.Random(seed)
    n=len(v)
    draws=sorted(mean(v[rng.randrange(n)] for _ in range(n))
                 for k in range(NBOOT))
    return [draws[49],draws[1949]]


def effect(values,seed):
    arr=list(values.values()) if isinstance(values,dict) else list(values)
    n=len(arr)
    ci=bootstrap(arr,seed)
    return {
        "qualified_original_tagged_animals":n,
        "mean":mean(arr) if ci else None,
        "source_original_tag_bootstrap95":ci,
        "nominal_interval_sign":(
            "POSITIVE" if ci and ci[0]>0 else
            "NEGATIVE" if ci and ci[1]<0 else
            "INCLUDES_ZERO" if ci else "UNAVAILABLE"
        ),
        "positive_original_tagged_animals":sum(x>0 for x in arr),
        "negative_original_tagged_animals":sum(x<0 for x in arr)
    }


def pair_native_q(a,b,seed):
    if len(a)<5 or len(b)<5:
        return {
            "native_tag_support":[len(a),len(b)],
            "first_minus_second_q":None,"original_independent_tag_bootstrap95":None,
            "nominal_interval_sign":"UNAVAILABLE"
        }
    first=sorted(a);second=sorted(b)
    rng=random.Random(seed)
    n=len(first);m=len(second)
    diff=mean(a[x] for x in first)-mean(b[x] for x in second)
    vals=sorted(
        mean(a[first[rng.randrange(n)]] for _ in range(n))-
        mean(b[second[rng.randrange(m)]] for _ in range(m))
        for j in range(NBOOT)
    )
    ci=[vals[49],vals[1949]]
    return {
        "native_tag_support":[n,m],
        "first_minus_second_q":diff,
        "original_independent_tag_bootstrap95":ci,
        "nominal_interval_sign":"POSITIVE" if ci[0]>0 else
                                "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
    }


def daily_source_q(group,dates):
    v=defaultdict(list)
    for _,years in group.items():
        for yr,by_date in years.items():
            for d,z in by_date.items():
                v[(yr,d)].append(z[0])
    if set(v)!=dates:
        raise ValueError("Not all common original dates represented by this registered taxon")
    return {k:mean(x) for k,x in v.items()}


def pair_season(source_a,source_b,dates,lat,seed):
    contrast={k:source_a[k]-source_b[k] for k in dates}
    x=[daylength(lat,d) for yr,d in dates]
    block=len({(yr,d.month) for yr,d in dates})
    acceptable=(len(dates)>=40 and block>=5 and max(x)-min(x)>=.5) if dates else False
    full=slope([(daylength(lat,d),contrast[(yr,d)]) for yr,d in dates]) if acceptable else None
    within=month_conditional_slope(contrast,lat,daylength,seed=seed)
    return {
        "source_shared_calendar_days":len(dates),
        "source_year_month_blocks":block,
        "study_reference_astronomical_daylength_hours_range":
            (max(x)-min(x)) if x else None,
        "all_yearmonth_reference_daylength_slope":full,
        "all_yearmonth_slope_status":
            "DESCRIPTIVE_ESTIMABLE" if full is not None else "UNAVAILABLE",
        "within_month_reference_daylength":within,
        "interpretation":"Reference astronomical daylength is determined by date and study GPS median latitude; a within-month slope is not biological photoperiod causality."
    }


def one_source(index,out):
    spec,manifest=fixed_sources()
    if index not in range(4):
        raise ValueError("Only four pre-event-registered original studies")
    study=spec["sources"][index]
    meta=manifest[index]
    taxa=study["registered_taxa"]
    result={
        "version":SOURCE_KEY,"date":"2026-10-29",
        "source_study_index":index,
        "original_study_DOI":study["doi"],
        "registered_native_taxa":taxa,
        "original_MD5_verified":False,
        "source_file_metadata_checked_again":False,
        "status":"ORIGINAL_SOURCE_HOLD",
        "original_native_taxon_results":[{"taxon":t,"status":"SOURCE_HOLD"}
                                       for t in taxa],
        "shared_source_calendar_support":None,
        "all_native_study_pair_shared_calendar_results":[],
        "previous_12_ERA5_source_version_HOLDS_UNCHANGED":True,
        "previous_48_source_negative_directional_radiation_replication_UNCHANGED":NEG,
        "raw_original_animal_identifiers_or_coordinates_exported":False
    }
    out.mkdir(parents=True,exist_ok=True)
    if meta["source_status"]!="SOURCE_ONLY_ORIGINAL_FILE_ADMISSIBLE":
        result["status"]="PRE_EVENT_SOURCE_METADATA_HOLD"
        result["source_hold_reason"]=meta.get("hold_reason")
    else:
        receipt=meta["source_original_file"]
        original=out/("private_original_source_"+str(index)+(
            ".zip" if receipt["filename"].lower().endswith(".zip") else
            ".csv.gz" if receipt["filename"].lower().endswith(".csv.gz") else ".csv"
        ))
        try:
            files=[x for x in original_files(study["package_uuid"])
                   if x["bitstream_uuid"]==receipt["bitstream_uuid"]]
            if (len(files)!=1 or files[0]["bytes"]!=receipt["bytes"]
                    or files[0]["checksum"]!=receipt["checksum"]):
                raise ValueError("Pre-event frozen provider original file identity not reproduced")
            result["source_file_metadata_checked_again"]=True
            if not checksum_and_download(files[0],original)["provider_checksum_verified"]:
                raise ValueError("Original public source checksum was not verified")
            result["original_MD5_verified"]=True
            tracks,coverage,other=load_source_taxa(original,receipt["filename"],taxa)
            status_by_taxon={}
            original_days={}
            for k,taxon in enumerate(taxa):
                first,_=one_taxon(taxon,tracks[taxon])
                original_days[taxon],daystatus=per_tag_year_days(tracks[taxon],"ORIGINAL_1H")
                status_by_taxon[taxon]=daystatus
                result["original_native_taxon_results"][k]={
                    "taxon":taxon,"status":"ORIGINAL_NATIVE_Q_EVALUATED",
                    "original_source_metadata_taxonomy_coverage":coverage[taxon],
                    "original_qualified_day_source_coverage":daystatus,
                    "original_q_full_source_tag_effect":first["day_night_speed_contrast"],
                    "original_q_daylength_full_source_tag_effect":first["daylength_association"],
                    "original_source_animals_with_diel_q":first["admitted_animals_diel"],
                    "original_source_animals_with_q_daylength":first["admitted_animals_photoperiod_slope"],
                }
            # Reference latitude is determined by telemetry positions, not
            # response q or species sign, and NEVER published as coordinate.
            latitude=source_reference_latitude(tracks)
            if latitude is None:
                result["status"]="SOURCE_NO_VALID_GPS_FIXES"
            else:
                source_calendar,dates,loops=common_original_calendar(original_days)
                source_monthblocks=len({(yr,d.month) for yr,d in dates})
                q_by_taxon={};grad_by_taxon={}
                for k,taxon in enumerate(taxa):
                    q,grad=tag_level_q_and_season(source_calendar[taxon],latitude)
                    q_by_taxon[taxon]=q
                    grad_by_taxon[taxon]=grad
                    t=result["original_native_taxon_results"][k]
                    t["same_study_all_native_common_calendar_q"]=effect(q,SEED+10*k+index)
                    t["same_study_all_native_common_calendar_daylength_slope"]=effect(grad,SEED+100+10*k+index)
                    t["native_original_source_tag_years_on_shared_calendar"]=sum(
                        len(y) for y in source_calendar[taxon].values()
                    )
                result["shared_source_calendar_support"]={
                    "source_calendar_dates_shared_by_all_registered_taxa":len(dates),
                    "source_yearmonth_blocks":source_monthblocks,
                    "registered_source_taxa":len(taxa),
                    "source_registered_native_q_tag_counts":
                        {t:len(q_by_taxon[t]) for t in taxa},
                    "fixed_point_shared_calendar_support_iterations":loops,
                    "source_shared_calendar_dates_and_GPS_coordinates_not_exported":True,
                    "reference_astronomical_latitude_exact_value_not_exported":True,
                }
                date_values=( {t:daily_source_q(source_calendar[t],dates)
                               for t in taxa} if dates else None )
                pairs=[]
                for j,(a,b) in enumerate(itertools.combinations(taxa,2)):
                    q=pair_native_q(q_by_taxon[a],q_by_taxon[b],
                                    SEED+1000*index+15*j)
                    weather=pair_season(date_values[a],date_values[b],dates,
                                        latitude,SEED+3000*index+15*j
                                        ) if dates else None
                    pairs.append({
                        "native_first":a,"native_second":b,
                        "same_source_calendar_dates_across_ALL_taxa":len(dates),
                        "source_original_animal_q_pair":q,
                        "source_common_day_q_difference_daylength":weather,
                        "not_proven_spatial_sympatry":True
                    })
                result["all_native_study_pair_shared_calendar_results"]=pairs
                result["status"]="ORIGINAL_MULTINATIVE_SHARED_CALENDAR_COMPLETED"
        except Exception as exc:
            result["status"]="ORIGINAL_SOURCE_OR_NATIVE_SHARED_CALENDAR_HOLD"
            result["source_hold_error_type"]=type(exc).__name__
        finally:
            original.unlink(missing_ok=True)
    target=out/("study_"+str(index)+"_original_multinative_q_v1.json")
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",
                      encoding="utf-8")
    print("RSDS_FOUR_ORIGINAL_MULTINATIVE_STUDY",json.dumps({
        "index":index,"doi":result["original_study_DOI"],
        "status":result["status"],"MD5":result["original_MD5_verified"],
        "native_taxa":len(result["original_native_taxon_results"]),
        "qualified_source_taxa_q":sum(z.get("original_q_full_source_tag_effect") is not None
                                     for z in result["original_native_taxon_results"]),
        "common_source_days":(result["shared_source_calendar_support"] or {}).get(
            "source_calendar_dates_shared_by_all_registered_taxa")
    },sort_keys=True),flush=True)
    return result


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--index",type=int,required=True)
    cli.add_argument("--out",type=Path,required=True)
    args=cli.parse_args()
    if args.index not in range(4):cli.error("Only four source-only original studies")
    args.out.mkdir(parents=True,exist_ok=True)
    one_source(args.index,args.out)


if __name__=="__main__":
    main()
