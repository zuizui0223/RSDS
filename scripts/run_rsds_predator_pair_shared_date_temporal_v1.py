"""Two original Movebank studies of sympatric coyotes and bobcats:
    relative diel GPS displacement and seasonal daylength association.

Washington and Illinois sources were previously frozen and already analyzed
for SOLAR HEADING in the original 12-source archive. This is an explicitly
post-source-selection NEW TEMPORAL q endpoint, not untouched external data.
No previously observed temporal q/sign or q~daylength source outcome is used
to choose which study, taxon, calendar dates, or GPS schedule to retain.

Primary: ORIGINAL one-record-per-UTC-hourbin q on study-SHARED astronomical
calendar dates of coyote (Canis latrans) vs bobcat (Lynx rufus).
Secondary: same ORIGINAL+EVEN+ODD supported dates and tagged animal-years,
only if three GPS duty-cycle schedules independently have enough segments.

Native Illinois white-tailed deer retained in the FULL source q/daylength
registry and is not silently dropped from the source census. All released
results are DOI/native taxon/group aggregate; original animal IDs, precise
GPS events, locations, calendar years, per-tag coefficients remain private.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download, daylength
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon, slope
from scripts.run_rsds_eight_study_gps_cadence_v1 import per_tag_year_days, PAIRS
from scripts.run_rsds_triple_matched_ecological_q_v1 import git_blob

CONTRACT=HERE/"research/RSDS_COYOTE_BOBCAT_TWO_STUDY_TIME_NICHE_V1_CONTRACT.json"
MANIFEST=HERE/"research/GLOBAL_SOLAR_AZIMUTH_12_NEW_SOURCE_FROZEN_ORIGINAL_FILE_MANIFEST_20261025.json"
SOLAR=HERE/"research/GLOBAL_SOLAR_AZIMUTH_12_NEW_SOURCE_FIRST_ORIGINAL_GPS_RESULT_20261025.json"
PRIOR=HERE/"research/RSDS_UTAH_WYOMING_TWO_STUDY_SHARED_DATE_PHOTOPERIOD_V1_RESULT_20261029.json"
TAXA=("Canis latrans","Lynx rufus")
SCHEDULES=("original","even","odd")
NBOOT=2000
SEED=20261029
REFERENCE=43.0
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c["version"]!="RSDS_COYOTE_BOBCAT_TWO_STUDY_TEMPORAL_NICHE_V1" or \
       len(c["source_studies"])!=2 or \
       any(x["primary_two"]!=list(TAXA) for x in c["source_studies"]):
        raise ValueError("Frozen source studies and new native two-taxon time ecology changed")
    docs={}
    for item in c["pinned_previous"]:
        target=HERE/item["path"]
        raw=target.read_bytes()
        if git_blob(raw)!=item["git_blob_SHA"]:
            raise ValueError("Original pre-existing Movebank source provenance changed")
        docs[item["path"]]=json.loads(raw.decode("utf-8"))
    m=docs[str(MANIFEST.relative_to(HERE))]
    old=docs[str(SOLAR.relative_to(HERE))]
    oldcmp=docs[str(PRIOR.relative_to(HERE))]
    if (len(m["sample_original_file_source_entries"])!=12 or
            old["registered_source_archives"]!=12 or
            oldcmp["fixed_original_study_DOIs"]!=2 or
            oldcmp["previous_48_source_shortwave_direction_external_nonreplication_UNCHANGED"]!=NEGATIVE):
        raise ValueError("First archived original 12-source GPS/previous ecological source frame changed")
    prior_file={x["package_uuid"]:x for x in m["sample_original_file_source_entries"]}
    prior_effect={x["package_uuid"]:x for x in old["all_registered_source_rows"]}
    if len({s["package_uuid"] for s in c["source_studies"]})!=2:
        raise ValueError("The original two native predator source DOI are not independent IDs")
    for s in c["source_studies"]:
        p,e=prior_file.get(s["package_uuid"]),prior_effect.get(s["package_uuid"])
        if p is None or e is None:
            raise ValueError("New native temporal source not in October 25 original source archive")
        orig=p.get("original_file") or {}
        if (p["status"]!="ORIGINAL_FILE_METADATA_ADMISSIBLE"
            or e["source_state"]!="original_MD5_verified_source_analyzed"
            or p["source_study_ids"]!=[s["study_id"]]
            or p["original_species_labels"]!=s["native_taxa"]
            or orig["bitstream_uuid"]!=s["bitstream_uuid"]
            or orig["bytes"]!=s["bytes"] or orig["checksum"]["value"].lower()!=s["md5"]):
            raise ValueError("Previous original source MD5/taxonomy/source-study frozen identity changed")
        if {x["taxon"] for x in e["source_taxa"]}!=set(s["native_taxa"]):
            raise ValueError("Original source taxa not preserved")
    return c,prior_effect


def check_original_source_coverage(tracks,source,reference,total):
    """No temporal q calculation can begin until original solar pilot support replay."""
    meta=reference["source_original_coverage"]
    if total["original_total_rows"]!=meta["original_total_rows"] or \
       total["source_rows_outside_registered_taxa"]!=meta["source_rows_outside_registered_taxa"]:
        raise ValueError("Original source parsed native GPS row count changed")
    if set(tracks)!=set(source["native_taxa"]):
        raise ValueError("Native species missing from original GPS archive")
    for taxon in source["native_taxa"]:
        expected=next(x for x in reference["source_taxa"] if x["taxon"]==taxon)
        cov=expected["source_original_coverage"]
        tag=tracks[taxon]
        hourly=sum(len(v) for v in tag.values())
        if hourly!=cov["hour_selected_fixes"] or len(tag)!=cov["source_individuals"]:
            raise ValueError("Prior original source GPS timestamp-selected fix/tag support changed")
    return True


def private_tag_calendar(tracks,mode):
    """Return private tag/year/calendar solar date -> schedule:q mappings."""
    labels=("original",) if mode=="original" else SCHEDULES
    all_raw={}
    for s in labels:
        gpsmode="ORIGINAL_1H" if s=="original" else PAIRS[0] if s=="even" else PAIRS[1]
        all_raw[s],_=per_tag_year_days(tracks,gpsmode)
    result={}
    common=set.intersection(*(set(a) for a in all_raw.values()))
    for tag in sorted(common):
        years=set.intersection(*(set(t[tag]) for t in all_raw.values()))
        selected={}
        for year in sorted(years):
            shared=set.intersection(*(set(t[tag][year]) for t in all_raw.values()))
            if len(shared)<10:
                continue
            dates=sorted(shared)
            x=[daylength(REFERENCE,d) for d in dates]
            if max(x)-min(x)<.5:
                continue
            selected[year]={d:{s:all_raw[s][tag][year][d][0] for s in labels}
                            for d in dates}
        if selected:
            result[tag]=selected
    return result


def same_native_calendar(first,second,labels):
    """Prune via SOURCE date eligibility only; never inspect q values to choose."""
    valid={TAXA[0]:first,TAXA[1]:second}
    rounds=0
    while True:
        rounds+=1
        if rounds>30:raise ValueError("Source common astronomy calendar selection did not converge")
        dates={}
        for species in TAXA:
            dates[species]={(yr,d) for tag in valid[species].values()
                            for yr,obs in tag.items() for d in obs}
        common=dates[TAXA[0]]&dates[TAXA[1]]
        changed={}
        for species in TAXA:
            dest={}
            for tag,years in valid[species].items():
                ys={}
                for yr,records in years.items():
                    remaining={d:x for d,x in records.items() if (yr,d) in common}
                    if len(remaining)>=10:
                        dl=[daylength(REFERENCE,d) for d in remaining]
                        if max(dl)-min(dl)>=.5:
                            ys[yr]=remaining
                if ys:dest[tag]=ys
            changed[species]=dest
        if changed==valid:
            break
        valid=changed
    result={s:{(yr,d) for ys in valid[s].values() for yr,ds in ys.items() for d in ds}
            for s in TAXA}
    if result[TAXA[0]]!=result[TAXA[1]]:
        raise ValueError("Source calendar dates inconsistent after fixed-point filtering")
    return valid,result[TAXA[0]],rounds


def tag_effects(pool,labels):
    level={}
    grads={}
    for tag,years in pool.items():
        qyear=[]
        gyear=[]
        for year,day in years.items():
            dates=sorted(day)
            if len(dates)<10:continue
            x=[daylength(REFERENCE,d) for d in dates]
            if max(x)-min(x)<.5:continue
            response={s:[day[d][s] for d in dates] for s in labels}
            qyear.append({s:mean(response[s]) for s in labels})
            slope_y={s:slope(list(zip(x,response[s]))) for s in labels}
            if all(k is not None for k in slope_y.values()):
                gyear.append(slope_y)
        if qyear:level[tag]={s:mean(z[s] for z in qyear) for s in labels}
        if gyear:grads[tag]={s:mean(z[s] for z in gyear) for s in labels}
    return level,grads


def summarize_tags(tags,labels,seed):
    names=sorted(tags)
    n=len(names)
    draws={s:[] for s in labels}
    if n>=5:
        rng=random.Random(seed)
        for _ in range(NBOOT):
            sample=[names[rng.randrange(n)] for k in range(n)]
            for s in labels:
                draws[s].append(mean(tags[t][s] for t in sample))
    result={}
    for s in labels:
        vals=[tags[t][s] for t in names]
        if n>=5:
            m=mean(vals)
            a=sorted(draws[s]);ci=[a[49],a[1949]]
            sign="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
        else:
            m=ci=None;sign="UNAVAILABLE"
        result[s]={"source_animal_mean":m,"source_tag_bootstrap95":ci,
                    "nominal_sign":sign,"original_tagged_animals":n}
    return result


def pair_contrast(left,right,labels,seed):
    a,b=sorted(left),sorted(right)
    n,m=len(a),len(b)
    rng=random.Random(seed)
    simulations={s:[] for s in labels}
    if n>=5 and m>=5:
        for _ in range(NBOOT):
            aa=[a[rng.randrange(n)] for _ in range(n)]
            bb=[b[rng.randrange(m)] for _ in range(m)]
            for s in labels:
                simulations[s].append(mean(left[t][s] for t in aa)-mean(right[t][s] for t in bb))
    response={}
    for s in labels:
        if n>=5 and m>=5:
            y=mean(left[t][s] for t in a)-mean(right[t][s] for t in b)
            z=sorted(simulations[s]);ci=[z[49],z[1949]]
            sign="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
        else:
            y=ci=None;sign="UNAVAILABLE"
        response[s]={"coyote_minus_bobcat_mean":y,"source_tag_bootstrap95":ci,
                     "nominal_sign":sign,"native_tags":{TAXA[0]:n,TAXA[1]:m}}
    return response


def calendar_date_slope(pool,dates,labels,seed):
    """Equal-date coyote-minus-bobcat q slope, year×month block bootstrap."""
    source={}
    for sp in TAXA:
        a=defaultdict(lambda:{s:[] for s in labels})
        for tag in pool[sp].values():
            for yr,obs in tag.items():
                for d,z in obs.items():
                    for s in labels:a[(yr,d)][s].append(z[s])
        source[sp]={key:{s:mean(vals[s]) for s in labels}
                    for key,vals in a.items()}
    if set(source[TAXA[0]])!=set(source[TAXA[1]]) or set(source[TAXA[0]])!=dates:
        raise ValueError("Source native date support differs between taxa")
    daydiff={x:{s:source[TAXA[0]][x][s]-source[TAXA[1]][x][s] for s in labels} for x in dates}
    blocks=defaultdict(list)
    for (yr,d),vals in daydiff.items():
        blocks[(yr,d.month)].append((daylength(REFERENCE,d),vals))
    keys=sorted(blocks)
    supported=len(dates)>=20 and len(keys)>=5 and (
        max(daylength(REFERENCE,d) for yr,d in dates)-
        min(daylength(REFERENCE,d) for yr,d in dates)>=.5)
    rng=random.Random(seed)
    effects={}
    for s in labels:
        value=slope([(daylength(REFERENCE,d),z[s]) for (yr,d),z in daydiff.items()]) if supported else None
        sampled=[]
        if supported and value is not None:
            for k in range(NBOOT):
                picks=[keys[rng.randrange(len(keys))] for v in keys]
                dataset=[rec for b in picks for rec in blocks[b]]
                val=slope([(v,y[s]) for v,y in dataset])
                if val is not None:sampled.append(val)
        if len(sampled)>=1600 and value is not None:
            ordered=sorted(sampled)
            ci=[ordered[int(.025*len(ordered))],ordered[int(.975*len(ordered))-1]]
            sign="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
        else:
            ci=None;sign="UNAVAILABLE"
        effects[s]={"same_calendar_coyote_minus_bobcat_q_daylength_slope":value,
                    "calendar_month_block_bootstrap95":ci,"nominal_sign":sign}
    return {"coobserved_source_calendar_dates":len(dates),
            "coobserved_source_year_month_blocks":len(keys),"slopes":effects}


def compare_study_pair(tracks,mode,seed):
    labels=("original",) if mode=="ORIGINAL_ONLY" else SCHEDULES
    pools={s:private_tag_calendar(tracks[s],("original" if mode=="ORIGINAL_ONLY" else "all_three"))
           for s in TAXA}
    matched,dates,iters=same_native_calendar(pools[TAXA[0]],pools[TAXA[1]],labels)
    result={"source_native_taxon_tag_support_on_shared_calendar":{s:len(matched[s]) for s in TAXA},
            "source_calendar_dates_with_both_taxa":len(dates),
            "source_month_year_blocks":len({(y,d.month) for y,d in dates}),
            "source_fixed_point_support_iterations":iters,
            "exact_source_UTC_schedule_labels":list(labels),
            "reference_astronomical_daylength_latitude":REFERENCE,
            "native":{},"native_pair":{}}
    records={}
    for i,sp in enumerate(TAXA):
        q,g=tag_effects(matched[sp],labels)
        records[sp]=(q,g)
        result["native"][sp]={"n_source_q_tags":len(q),"n_source_daylength_tags":len(g),
             "q":summarize_tags(q,labels,seed+i*100),
             "q_daylength":summarize_tags(g,labels,seed+i*100+3)}
    a,b=(records[sp] for sp in TAXA)
    result["native_pair"]={
        "coyote_minus_bobcat_q":pair_contrast(a[0],b[0],labels,seed+501),
        "coyote_minus_bobcat_daylength":pair_contrast(a[1],b[1],labels,seed+701),
    }
    result["date_level_pair_difference_vs_daylength"]=(
        calendar_date_slope(matched,dates,labels,seed+900) if dates else {
            "coobserved_source_calendar_dates":0,
            "coobserved_source_year_month_blocks":0,
            "slopes":{s:{"nominal_sign":"UNAVAILABLE",
                         "same_calendar_coyote_minus_bobcat_q_daylength_slope":None,
                         "calendar_month_block_bootstrap95":None} for s in labels}
        })
    result["status"]="BOTH_NATIVE_TAXA_SOURCE_COOCCURRENT_ESTIMABLE" if min(
        len(records[s][0]) for s in TAXA)>=5 else "INSUFFICIENT_COOCCURRENT_SOURCE_NATIVE_TAGS"
    return result


def run_source(region,out):
    spec,prior=frozen()
    s=next(x for x in spec["source_studies"] if x["region"]==region)
    old=prior[s["package_uuid"]]
    result={"version":"RSDS_COYOTE_BOBCAT_SHARED_DATE_TEMPORAL_ORIGINAL_STUDY_V1",
            "source_region":region,"original_study_DOI":s["doi"],
            "original_source_ID":s["study_id"],
            "native_taxa":s["native_taxa"],
            "status":"ORIGINAL_SOURCE_HOLD","provider_MD5_verified":False,
            "first_solar_source_native_hour_fix_counts_replayed":False,
            "original_native_source_full_temporal_q_and_daylength":[],
            "shared_original_calendar_comparison":None,
            "shared_triple_GPS_clock_comparison":None,
            "old_12_ERA5_RAW_SHA_HOLDS_UNCHANGED":True,
            "old_48_source_sun_direction_external_nonreplication_UNCHANGED":NEGATIVE}
    local=out/("private_rsds_predator_source_"+region.lower()+".csv")
    try:
        files=[t for t in original_files(s["package_uuid"]) if t["bitstream_uuid"]==s["bitstream_uuid"]]
        if (len(files)!=1 or files[0]["bytes"]!=s["bytes"] or
                files[0]["checksum"]["value"].lower()!=s["md5"].lower()):
            raise ValueError("Original pre-event source-file UUID MD5/size changed")
        recorded=checksum_and_download(files[0],local)
        if not recorded["provider_checksum_verified"]:
            raise ValueError("Original GPS MD5 provider source file did not match")
        result["provider_MD5_verified"]=True
        source,coverage,total=load_source_taxa(local,s["file_name"],s["native_taxa"])
        check_original_source_coverage(source,s,old,total)
        result["first_solar_source_native_hour_fix_counts_replayed"]=True
        print("RSDS_COYOTE_BOBCAT_PROVIDER_GPS_REPLAY",region,len(source),flush=True)
        for species in s["native_taxa"]:
            q,detail=one_taxon(species,source[species])
            result["original_native_source_full_temporal_q_and_daylength"].append({
                "taxon":species,
                "original_hour_selected_fix_count":coverage[species]["hour_selected_fixes"],
                "original_tag_count":coverage[species]["source_individuals"],
                "q":q["day_night_speed_contrast"],
                "q_daylength":q["daylength_association"],
                "source_q_original_tag_support":q["admitted_animals_diel"],
                "source_daylength_original_tag_support":q["admitted_animals_photoperiod_slope"],
            })
        pair={k:source[k] for k in TAXA}
        result["shared_original_calendar_comparison"]=compare_study_pair(pair,"ORIGINAL_ONLY",SEED)
        result["shared_triple_GPS_clock_comparison"]=compare_study_pair(pair,"THREE_SCHEDULES",SEED+10000)
        result["status"]="ORIGINAL_SOURCE_SHARED_CALENDAR_EVALUATED"
    except Exception as e:
        result["status"]="ORIGINAL_SOURCE_OR_SHARED_CALENDAR_HOLD"
        result["hold_error_type"]=type(e).__name__
    finally:
        local.unlink(missing_ok=True)
    if not result["original_native_source_full_temporal_q_and_daylength"]:
        result["original_native_source_full_temporal_q_and_daylength"]=[
            {"taxon":sp,"q":None,"q_daylength":None} for sp in s["native_taxa"]]
    out.mkdir(parents=True,exist_ok=True)
    (out/("source_"+region.lower()+"_predator_time.json")).write_text(
        json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_COYOTE_BOBCAT_TEMPORAL_SOURCE",json.dumps({
        "region":region,"status":result["status"],
        "MD5":result["provider_MD5_verified"],
        "native_taxa":len(result["native_taxa"]),
        "full_source_q_tags":{
            x["taxon"]:x.get("source_q_original_tag_support")
            for x in result["original_native_source_full_temporal_q_and_daylength"]},
        "calendar_original_support":
            (result["shared_original_calendar_comparison"] or {}).get("source_native_taxon_tag_support_on_shared_calendar"),
        "triple_support":
            (result["shared_triple_GPS_clock_comparison"] or {}).get("source_native_taxon_tag_support_on_shared_calendar")
    },sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region",choices=("Washington","Illinois"),required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    run_source(a.region,a.out)


if __name__=="__main__":
    main()
