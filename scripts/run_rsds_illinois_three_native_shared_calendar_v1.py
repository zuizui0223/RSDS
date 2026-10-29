"""Illinois original Movebank source: coobserved coyote/bobcat/white-tailed deer q.

This is an OUTCOME-AWARE post-PR16 source support restriction (previous
native source q means already seen), not a new independent field study or
preregistered predator-prey behavior test. Use exact published source provider
UUID/MD5, replay all three earlier native q/daylength source means and
individual counts before source-common-date restriction.

No individual GPS ID/position/UTC event/precise calendar date leaves memory.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from scripts.movebank_native_ecology_20261019 import HERE,original_files,checksum_and_download,daylength
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon,slope
from scripts.run_rsds_predator_pair_shared_date_temporal_v1 import (
    frozen as previous_frozen, check_original_source_coverage,
    private_tag_calendar,tag_effects,summarize_tags
)
from scripts.run_rsds_triple_matched_ecological_q_v1 import git_blob

CONTRACT=HERE/"research/RSDS_ILLINOIS_THREE_NATIVE_SHARED_DATE_V1_CONTRACT.json"
PREVIOUS=HERE/"research/RSDS_COYOTE_BOBCAT_TWO_STUDY_SHARED_DATE_TEMPORAL_V1_FIRST_RESULT_20261029.json"
TAXA=("Canis latrans","Lynx rufus","Odocoileus virginianus")
PAIRS=((TAXA[0],TAXA[1]),(TAXA[0],TAXA[2]),(TAXA[1],TAXA[2]))
REFERENCE=43.0
SEED=20261029
NBOOT=2000
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def freeze():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (spec["version"]!="RSDS_ILLINOIS_THREE_NATIVE_SYMPATRIC_DATE_Q_V1"
            or spec["source_original_study"]["native_taxa"]!=list(TAXA)):
        raise ValueError("Precommitted full three native Illinois study source changed")
    older,source=previous_frozen()
    receipt=PREVIOUS.read_bytes()
    pins={x["path"]:x["git_blob_SHA"] for x in spec["pinned_previous"]}
    if git_blob(receipt)!=pins[str(PREVIOUS.relative_to(HERE))]:
        raise ValueError("First Illinois native q and two predator original source results altered")
    before=json.loads(receipt.decode("utf-8"))
    if (before["original_study_DOI_count"]!=2 or before["native_study_taxon_rows"]!=5 or
        before["previous_48_source_negative_direction_shortwave_external_replication_UNCHANGED"]!=NEGATIVE):
        raise ValueError("Previously accepted original source study/taxa and solar heldout verdict changed")
    one=next(x for x in before["all_original_source_study_native_taxa"]
             if x["source_region"]=="Illinois")
    study=spec["source_original_study"]
    if one["source_original_study_DOI"]!=study["doi"] or set(one["source_native_taxa"])!=set(TAXA):
        raise ValueError("Original native Illinois source taxa changed")
    source_old=source[study["package_uuid"]]
    if set(x["taxon"] for x in source_old["source_taxa"])!=set(TAXA):
        raise ValueError("Prior original solar source native taxonomy changed")
    return spec,one,source_old


def triple_common_calendar(private):
    """Fixed point intersection of source calendar date supports in ALL 3 taxa.

    Source original q VALUES NEVER enter the support selection criterion.
    """
    if set(private)!=set(TAXA):
        raise ValueError("Exactly three original co-tracked native taxa required")
    valid=private
    loops=0
    while True:
        loops+=1
        if loops>50:
            raise ValueError("Triply coobserved original source calendar did not converge")
        pools={sp:{(yr,d) for tagged in valid[sp].values()
                   for yr,obs in tagged.items() for d in obs}
               for sp in TAXA}
        same=set.intersection(*(pools[sp] for sp in TAXA))
        next_state={}
        for sp in TAXA:
            group={}
            for tag,years in valid[sp].items():
                kept={}
                for yr,record in years.items():
                    shared={d:v for d,v in record.items() if (yr,d) in same}
                    if len(shared)>=10:
                        dx=[daylength(REFERENCE,d) for d in shared]
                        if max(dx)-min(dx)>=.5:
                            kept[yr]=shared
                if kept:
                    group[tag]=kept
            next_state[sp]=group
        if next_state==valid:
            break
        valid=next_state
    dates={sp:{(yr,d) for tag in valid[sp].values()
               for yr,rec in tag.items() for d in rec}
           for sp in TAXA}
    if any(dates[sp]!=dates[TAXA[0]] for sp in TAXA[1:]):
        raise ValueError("Three native taxa do not share identical original source astronomy dates")
    return valid,dates[TAXA[0]],loops


def pair_delta(a,b,seed):
    aa,bb=sorted(a),sorted(b)
    n,m=len(aa),len(bb)
    if min(n,m)<5:
        return {"source_tag_counts":[n,m],"status":"UNAVAILABLE",
                "first_minus_second_mean":None,"original_tag_bootstrap95":None}
    center=mean(a[x]["original"] for x in aa)-mean(b[x]["original"] for x in bb)
    rng=random.Random(seed)
    draws=[]
    for _ in range(NBOOT):
        left=[aa[rng.randrange(n)] for _ in range(n)]
        right=[bb[rng.randrange(m)] for _ in range(m)]
        draws.append(mean(a[t]["original"] for t in left)-
                     mean(b[t]["original"] for t in right))
    draws.sort()
    ci=[draws[49],draws[1949]]
    return {"source_tag_counts":[n,m],"status":(
        "POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"),
        "first_minus_second_mean":center,
        "original_tag_bootstrap95":ci}


def datewise_pair_slope(a,b,dates,seed):
    if not dates:
        return {"source_calendar_days":0,"source_year_month_blocks":0,
                "slope":None,"date_block95":None,"status":"UNAVAILABLE"}
    dated={}
    for group in (a,b):
        means=defaultdict(list)
        for year_set in group.values():
            for year,days in year_set.items():
                for d,x in days.items(): means[(year,d)].append(x["original"])
        dated[id(group)]={key:mean(vals) for key,vals in means.items()}
    if any(set(dated[id(x)])!=dates for x in (a,b)):
        raise ValueError("Source-native datewise observation sets differ within shared date interval")
    delta={key:dated[id(a)][key]-dated[id(b)][key] for key in dates}
    block=defaultdict(list)
    for (year,d),x in delta.items():
        block[(year,d.month)].append((daylength(REFERENCE,d),x))
    keys=sorted(block)
    good=(len(dates)>=20 and len(keys)>=5 and
          max(daylength(REFERENCE,d) for yr,d in dates)-
          min(daylength(REFERENCE,d) for yr,d in dates)>=.5)
    coeff=slope([(daylength(REFERENCE,d),v) for (yr,d),v in delta.items()]) if good else None
    if coeff is None:
        return {"source_calendar_days":len(dates),"source_year_month_blocks":len(keys),
                "slope":None,"date_block95":None,"status":"UNAVAILABLE"}
    rng=random.Random(seed)
    vals=[]
    for k in range(NBOOT):
        sampled=[keys[rng.randrange(len(keys))] for item in keys]
        y=slope([x for part in sampled for x in block[part]])
        if y is not None:vals.append(y)
    if len(vals)<1600:
        return {"source_calendar_days":len(dates),"source_year_month_blocks":len(keys),
                "slope":coeff,"date_block95":None,"status":"BLOCK_BOOTSTRAP_UNSTABLE"}
    vals.sort();ci=[vals[int(.025*len(vals))],vals[int(.975*len(vals))-1]]
    return {"source_calendar_days":len(dates),"source_year_month_blocks":len(keys),
            "slope":coeff,"date_block95":ci,"status":(
                "POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO")}


def run(out):
    spec,first,old_solar=freeze()
    src=spec["source_original_study"]
    result={
        "version":"RSDS_ILLINOIS_ALL_THREE_ORIGINAL_NATIVE_SHARED_DATE_Q_V1_RESULT",
        "date":"2026-10-29","study_DOI":src["doi"],"native_taxa":list(TAXA),
        "status":"ORIGINAL_SOURCE_HOLD","MD5_verified":False,
        "first_all_native_q_daylength_and_solar_source_replayed":False,
        "source_all_native_q_levels_on_same_original_calendar":None,
        "source_all_three_native_support":None,
        "three_source_native_pair_q_and_daylength":None,
        "previous_washington_temporal_no_eligible_tags_HOLD_UNCHANGED":True,
        "old_12_ERA5_raw_SHA_HOLDS_UNCHANGED":True,
        "prior_48_negative_shortwave_direction_external_UNCHANGED":NEGATIVE
    }
    out.mkdir(parents=True,exist_ok=True)
    file=out/"private_illinois_three_native_original.csv"
    try:
        fs=[f for f in original_files(src["package_uuid"])
            if f["bitstream_uuid"]==src["bitstream_uuid"]]
        if len(fs)!=1 or fs[0]["bytes"]!=src["bytes"] or \
           fs[0]["checksum"]["value"].lower()!=src["md5"]:
            raise ValueError("Original Illinois provider source file UUID/size/MD5 changed")
        verification=checksum_and_download(fs[0],file)
        if not verification["provider_checksum_verified"]:
            raise ValueError("Original Illinois source provider MD5 invalid")
        result["MD5_verified"]=True
        tracks,coverage,total=load_source_taxa(file,src["file_name"],list(TAXA))
        check_original_source_coverage(tracks,src,old_solar,total)
        previous=first["source_native_taxa"]
        for sp in TAXA:
            native,_=one_taxon(sp,tracks[sp])
            a=previous[sp]
            for name,field in (("day_night_speed_contrast","q"),
                               ("daylength_association","q_daylength")):
                z=native[name]
                r=a[field]
                if z is None and r["mean"] is not None:
                    raise ValueError("Original source native q or q/daylength first result changed")
                if z is not None:
                    if z["n_individuals"]!=r["source_tags"] or abs(z["mean"]-r["mean"])>1e-12:
                        raise ValueError("Native Illinois source first q/daylength and tags not replayed")
        result["first_all_native_q_daylength_and_solar_source_replayed"]=True
        print("RSDS_ILLINOIS_3NATIVE_PRIOR_Q_AND_MD5_REPLAY",len(TAXA),flush=True)
        pools={sp:private_tag_calendar(tracks[sp],"original") for sp in TAXA}
        aligned,dates,loops=triple_common_calendar(pools)
        q={};season={};native={}
        for i,sp in enumerate(TAXA):
            q[sp],season[sp]=tag_effects(aligned[sp],("original",))
            native[sp]={
                "q_source_tagged_animals":len(q[sp]),
                "daylength_source_tagged_animals":len(season[sp]),
                "q_original":summarize_tags(q[sp],("original",),SEED+101*i)["original"],
                "daylength_original":summarize_tags(season[sp],("original",),SEED+201*i)["original"]
            }
        result["source_all_three_native_support"]={
            "native_original_tag_count_on_all_three_shared_source_dates":
                {sp:len(aligned[sp]) for sp in TAXA},
            "shared_source_calendar_dates_all_three_native_taxa":len(dates),
            "shared_source_year_month_blocks":len({(yr,d.month) for yr,d in dates}),
            "fixed_point_source_support_iterations":loops,
            "reference_astronomical_daylength_latitude_degrees":REFERENCE,
            "original_calendar_years_and_dates_not_exported":True
        }
        result["source_all_native_q_levels_on_same_original_calendar"]=native
        pair={}
        for i,(a,b) in enumerate(PAIRS):
            tag=f"{a} MINUS {b}"
            pair[tag]={
                "native_first":a,"native_second":b,
                "q":pair_delta(q[a],q[b],SEED+701+i),
                "q_daylength":pair_delta(season[a],season[b],SEED+801+i),
                "daily_mean_q_difference_vs_standard_daylength":
                    datewise_pair_slope(aligned[a],aligned[b],dates,SEED+901+i)
            }
        result["three_source_native_pair_q_and_daylength"]=pair
        result["status"]="THREE_NATIVE_SHARED_SOURCE_DATE_EVALUATED" if min(
            len(q[sp]) for sp in TAXA)>=5 else "INSUFFICIENT_THREE_NATIVE_COMMON_DATE_TAGS"
    except Exception as exc:
        result["status"]="ORIGINAL_ILLINOIS_SOURCE_OR_SHARED_DATE_HOLD"
        result["source_hold_error_type"]=type(exc).__name__
    finally:
        file.unlink(missing_ok=True)
    (out/"illinois_three_native_shared_date_result.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_ILLINOIS_3NATIVE_SHARED_SOURCE_DATE",json.dumps({
        "status":result["status"],"original_MD5":result["MD5_verified"],
        "prior_native_q_replay":result["first_all_native_q_daylength_and_solar_source_replayed"],
        "shared_dates":(result["source_all_three_native_support"] or {}).get(
          "shared_source_calendar_dates_all_three_native_taxa"),
        "shared_native_tags":(result["source_all_three_native_support"] or {}).get(
          "native_original_tag_count_on_all_three_shared_source_dates")
    },sort_keys=True),flush=True)
    return result


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,required=True)
    a=cli.parse_args()
    run(a.out)


if __name__=="__main__":
    main()
