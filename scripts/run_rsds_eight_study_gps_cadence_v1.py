"""Study-stratified GPS observation-cadence sensitivity on SAME tagged animals/dates.

Do not interpret as animals changing activity: observed straight-line speed
depends on temporal sample spacing, daypart boundary crossing and source-tag
detection effort. Original one-fix/UTC-hour sources and MD5 files are frozen.
For 2h duty-cycle alternate parity, use EVEN and ODD original UTC hour fix
bins with actual segment elapsed 1.5-2.5h, no interpolation. Pair q/slope
only on the same original tag/year/date records in both resolutions.

All original 8 source DOIs, 12 study-taxon rows, one fixed egret HOLD, negative
old external shortwave direction result, and 12 old ERA5 source-version HOLDS
remain immutable. Only aggregated source tagged-animal effects are released.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from statistics import mean, median

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download, daylength, solar_period
)
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon, slope
from scripts.movebank_solar_movement_20261019 import (
    valid_interval, great_circle_km
)
from scripts.movebank_utah_four_taxa_solar_20261019 import (
    extract_time_only_original_source
)
from scripts.movebank_wyoming_two_ungulates_20261022 import (
    streaming_time_selected_csv
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.run_rsds_same_animal_solar_weather_v0 import source_meta
from scripts.run_rsds_new_six_original_gps_ecology_v0 import load_frozen as six_sources

CONTRACT=HERE/"research/RSDS_ORIGINAL_EIGHT_STUDY_GPS_CADENCE_V1_CONTRACT.json"
PREV_2=HERE/"research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
PREV_6=HERE/"research/RSDS_NEW_SIX_ORIGINAL_GPS_ECOLOGY_FIRST_RESULT_20261026.json"
METADATA_PIN=HERE/"research/RSDS_NEW_SIX_ORIGINAL_GPS_SOURCE_METADATA_BLOB_V0.json"
ATLAS=HERE/"research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
SEED=20261028

PAIRS=("EVEN_UTC_HOURS","ODD_UTC_HOURS")


def git_blob(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()


def verify_frame():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if c["version"]!="RSDS_ORIGINAL_8_STUDY_UTC_CADENCE_MATCHED_DAY_Q_V1":
        raise ValueError("Frozen GPS original cadence contract version invalid")
    f=c["frozen_previous_blobs"]
    checks={
        PREV_2:f["original_two_study_source_first"],
        PREV_6:f["original_new_six_study_effects"],
        METADATA_PIN:f["original_new_six_primary_file_SHA_receipt"],
        ATLAS:f["original_eight_study_macro"],
    }
    data={}
    for path,sha in checks.items():
        raw=path.read_bytes()
        if git_blob(raw)!=sha:
            raise ValueError("Original published GPS source first result has changed")
        data[path]=json.loads(raw.decode("utf-8"))
    old, new, _, macro=(data[x] for x in (PREV_2,PREV_6,METADATA_PIN,ATLAS))
    if old["fixed_original_study_count"]!=2 or old["total_source_native_taxon_rows"]!=6:
        raise ValueError("Original two study/taxon frame changed")
    if new["frozen_source_original_study_count"]!=6 or macro["original_study_denominator"]!=8:
        raise ValueError("First original six and eight study cohorts changed")
    if (old["frozen_external_radiation_nonreplication_UNCHANGED"]!=NEGATIVE
            or new["previous_external_48_study_solar_radiation_nonreplication_UNCHANGED"]!=NEGATIVE
            or macro["old_48_radiation_study_external_negative_direction_test_UNCHANGED"]!=NEGATIVE
            or not macro["old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED"]):
        raise ValueError("Old null/heldout ecological results were reclassified")
    six,metadata=six_sources()
    if len(six)!=6 or len(metadata)!=6:
        raise ValueError("Missing one original source or source metadata HOLD")
    return old,new,six,metadata


def segment(a,b,mode):
    """Return native (astronomical day,part,km/h,lat) or no detected interval."""
    if mode=="ORIGINAL_1H":
        out=valid_interval(a[0],a[1],a[2],b[0],b[1],b[2])
        return (out[0],out[1],out[2],out[4]) if out is not None else None
    hours=(b[0]-a[0]).total_seconds()/3600
    if not 1.5<=hours<=2.5:return None
    midpoint=a[0]+(b[0]-a[0])/2
    lat=(a[1]+b[1])/2
    lon=(a[2]+b[2])/2
    before=solar_period(a[0],a[1],a[2])
    after=solar_period(b[0],b[1],b[2])
    center=solar_period(midpoint,lat,lon)
    if before!=after or before!=center:return None
    velocity=great_circle_km(a[1],a[2],b[1],b[2])/hours
    if not math.isfinite(velocity) or velocity>60:return None
    return center[0],center[1],velocity,lat


def per_tag_year_days(animals,mode):
    """Private source original tag-year-day tables, NEVER serialized."""
    if mode not in ("ORIGINAL_1H",)+PAIRS:
        raise ValueError("Unregistered UTC cadence parity")
    out={}
    accepted_segments=proposed=0
    qualified_days=0
    source_hourly_fixes=0
    for tag, fix_by_hour in animals.items():
        source_hourly_fixes += len(fix_by_hour)
        original=sorted(
            (value for utc,value in fix_by_hour.items()
             if mode=="ORIGINAL_1H" or
             utc.hour%2==(0 if mode=="EVEN_UTC_HOURS" else 1)),
            key=lambda x:x[0]
        )
        days=defaultdict(lambda:{"day":[],"night":[],"lat":[]})
        for a,b in zip(original,original[1:]):
            proposed+=1
            val=segment(a,b,mode)
            if val is None:continue
            accepted_segments+=1
            date,part,speed,lat=val
            z=days[date]
            z[part].append(speed)
            z["lat"].append(lat)
        year=defaultdict(dict)
        for d,entry in days.items():
            if len(entry["day"])<3 or len(entry["night"])<3:
                continue
            qualified_days+=1
            avg_lat=mean(entry["lat"])
            q=math.log1p(median(entry["night"])) - math.log1p(median(entry["day"]))
            year[d.year][d]=(q,daylength(avg_lat,d))
        if year:
            out[tag]=year
    return out,{
        "original_selected_1hour_source_fixes":source_hourly_fixes,
        "candidate_adjacent_GPS_intervals":proposed,
        "qualified_same_solar_period_GPS_intervals":accepted_segments,
        "animal_calendar_days_with_ge3_day_and_ge3_night_segments":qualified_days,
        "tags_with_some_q_observable":len(out)
    }


def source_baseline_check(animals,reference):
    """Fail BEFORE opening parity if base q/native slope is not exact."""
    qrows=[]
    lightrows=[]
    for taxon, animal_data in animals.items():
        observed,_=one_taxon(taxon,animal_data)
        expected=reference.get(taxon)
        if expected is None:
            raise ValueError("Native original study taxon not pinned in prior source")
        for item in ("day_night_speed_contrast","daylength_association"):
            val=observed[item]
            archived=expected[item]
            if (val is None)!=(archived is None):
                raise ValueError("Original source native biology effect eligibility changed: "+item)
            if val is not None and (val["n_individuals"]!=archived["n_individuals"]
                    or abs(val["mean"]-archived["mean"])>1e-12):
                raise ValueError("Original GPS native source q/daylength means differ from first result")
        qrows.append(observed["day_night_speed_contrast"])
        lightrows.append(observed["daylength_association"])
    return True


def matched_per_tag(full,alternate):
    """Pair by original source tag/year/solar date; values stay private."""
    differences_q={}
    differences_slope={}
    count_years=count_days=0
    source_phase_tag_years=0
    for tag in set(full)&set(alternate):
        qs=[]
        slopes=[]
        for year in sorted(set(full[tag])&set(alternate[tag])):
            same=sorted(set(full[tag][year]) & set(alternate[tag][year]))
            if len(same)<10:
                continue
            count_years+=1
            count_days+=len(same)
            qbase=[full[tag][year][d][0] for d in same]
            qthin=[alternate[tag][year][d][0] for d in same]
            qs.append(mean(qthin)-mean(qbase))
            daylen=[full[tag][year][d][1] for d in same]
            x=slope(list(zip(daylen,qbase)))
            y=slope(list(zip(daylen,qthin)))
            if x is not None and y is not None:
                slopes.append(y-x)
        if qs:
            differences_q[tag]=mean(qs)
            source_phase_tag_years+=len(qs)
        if slopes:
            differences_slope[tag]=mean(slopes)
    return differences_q,differences_slope,{
        "source_tag_years_with_ge10_matched_daypart_days":count_years,
        "original_tagged_days_identical_between_1h_and_2h":count_days,
        "source_tagged_animals_with_ge10_matching_days":len(differences_q),
        "source_tagged_animals_with_identifiable_matched_q_daylength_slope":len(differences_slope),
    }


def tag_bootstrap(values,seed):
    if len(values)<5:
        return None
    rng=random.Random(seed)
    n=len(values)
    b=sorted(mean(values[rng.randrange(n)] for _ in range(n))
             for _ in range(2000))
    return [b[49],b[1949]]


def summary(deltas,seed):
    vals=list(deltas.values())
    n=len(vals)
    ci=tag_bootstrap(vals,seed)
    return {
        "n_original_tagged_animals_paired_between_cadences":n,
        "status":"SOURCE_TAGGED_ANIMAL_DESCRIPTIVE" if n>=5 else "INSUFFICIENT_MATCHED_TAGS",
        "paired_tag_mean_2h_minus_same_day_1h":mean(vals) if n>=5 else None,
        "original_tag_resample95":ci,
        "n_positive_tag_deltas":sum(x>0 for x in vals),
        "n_negative_tag_deltas":sum(x<0 for x in vals),
        "nominal_interval_sign":(
            "POSITIVE" if ci and ci[0]>0 else
            "NEGATIVE" if ci and ci[1]<0 else
            "INCLUDES_ZERO" if ci else "UNAVAILABLE"
        )
    }


def analyze_taxon(taxon,animals):
    full,n0=per_tag_year_days(animals,"ORIGINAL_1H")
    ph={}
    for k,mode in enumerate(PAIRS):
        alt,attr=per_tag_year_days(animals,mode)
        dq,dl,common=matched_per_tag(full,alt)
        ph[mode]={
            "source_interval_gap_hours":"1.5_to_2.5",
            "animal_date_overlap_rule":"ge10_same_dates_per_calendar_year",
            "data_support":attr,
            "identical_original_tag_date_support":common,
            "matched_q_velocity_proxy_difference":summary(dq,SEED+17*k),
            "matched_q_daylength_slope_difference_per_hour":summary(dl,SEED+17*k+2),
        }
    return {"taxon":taxon,
            "source_original_hourly_eligible_GPS_cadence":"0.5_to_2.0",
            "source_original_full_hour_baseline_support":n0,
            "source_original_even_odd_2h_phase_sensitivity":ph,
            "not_latent_activity_duration_or_environmental_treatment":True}


def source_reference(region,old,new):
    if region in ("Utah","Wyoming"):
        parent=next(x for x in old["sources"] if x["region"]==region)
        src,taxa,baseline=source_meta(region)
        expected={}
        for taxon in taxa:
            n0,q,n1,light=baseline[taxon]
            # Original q/daylight source tagged animal CI fields appear only in
            # source first summary, not in this preflight scalar tuple.
            expected[taxon]={
                "day_night_speed_contrast":{"n_individuals":n0,"mean":q},
                "daylength_association":{"n_individuals":n1,"mean":light}
            }
        if src["doi"]!=parent["original_source_doi"] or len(taxa)!=(4 if region=="Utah" else 2):
            raise ValueError("Original Utah/Wyoming study DOI changed")
        return src,taxa,expected,None
    i=int(region)-2
    if not 0<=i<6:
        raise ValueError("Cadence study index outside previously frozen eight original study DOIs")
    fixed,meta=six_sources()
    src=fixed[i]
    parent=new["six_full_source_study_taxon_results"][i]
    if parent["study_doi"]!=src["doi"] or parent["native_taxon"]!=src["taxon"]:
        raise ValueError("Six original pre-event source-study identity changed")
    if meta[i]["source_metadata_status"]!="ORIGINAL_FILE_METADATA_ADMISSIBLE":
        return src,[src["taxon"]],{},meta[i]
    before={}
    for key,col in (
        ("day_night_speed_contrast","source_q_night_minus_day"),
        ("daylength_association","source_q_daylength_per_hour")
    ):
        x=parent[col]
        before[key]={"n_individuals":x["n_animals"],"mean":x["mean"]} \
            if x["mean"] is not None else None
    return src,[src["taxon"]],{src["taxon"]:before},meta[i]


def run_source(index,out_dir):
    old,new,fixed,meta=verify_frame()
    region=("Utah" if index==0 else "Wyoming" if index==1 else str(index))
    src,taxa,reference,preflight=source_reference(region,old,new)
    is_old=index<2
    doi=src["doi"] if is_old else src["doi"]
    output={
        "version":"RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_V1_SOURCE",
        "source_index":index,"original_study_doi":doi,
        "source_native_taxa":list(taxa),
        "source_status":"SOURCE_HOLD",
        "original_GPS_provider_MD5_verified":False,
        "baseline_1h_q_and_daylength_source_effects_replayed":False,
        "native_taxon_results":[],
        "unchanged_old_ERA5_12_RAW_SHA_HOLDS":True,
        "unchanged_48_source_signed_radiation_external_NONREPLICATION":NEGATIVE
    }
    if preflight and preflight["source_metadata_status"]!="ORIGINAL_FILE_METADATA_ADMISSIBLE":
        output.update(source_status="PREDECLARED_FILE_METADATA_HOLD",
                      hold_reason=preflight.get("metadata_hold_reason"))
    else:
        bit=(src["original_bitstream_uuid"] if is_old else preflight["original_file"]["bitstream_uuid"])
        pkg=(src["package_uuid"] if is_old else src["id"])
        length=(src["original_bytes"] if is_old else preflight["original_file"]["bytes"])
        md5=(src["original_md5"] if is_old else preflight["original_file"]["checksum"]["value"])
        suffix=(".zip" if (index==0 or
                (preflight is not None and preflight["original_file"]["filename"].lower().endswith(".zip")))
                else ".csv.gz" if preflight is not None and preflight["original_file"]["filename"].lower().endswith(".gz")
                else ".csv")
        local=out_dir/("private_source_"+str(index)+suffix)
        try:
            sources=[x for x in original_files(pkg) if x["bitstream_uuid"]==bit]
            if (len(sources)!=1 or sources[0]["bytes"]!=length or
                    sources[0]["checksum"]["value"].lower()!=md5.lower()):
                raise ValueError("Previously frozen original study source UUID/size/MD5 changed")
            recorded=checksum_and_download(sources[0],local)
            if not recorded["provider_checksum_verified"]:
                raise ValueError("Original GPS source MD5 could not be verified")
            output["original_GPS_provider_MD5_verified"]=True
            if index==0:
                animals,counts,other=extract_time_only_original_source(local,taxa)
            elif index==1:
                animals,counts,other=streaming_time_selected_csv(local,taxa)
            else:
                animals,counts,other=load_source_taxa(local,preflight["original_file"]["filename"],taxa)
            source_baseline_check(animals,reference)
            output["baseline_1h_q_and_daylength_source_effects_replayed"]=True
            print("RSDS_CADENCE_ORIGINAL_BASELINE_VERIFIED",json.dumps({
                "index":index,"n_native_taxa":len(taxa)},sort_keys=True),flush=True)
            for taxon in taxa:
                output["native_taxon_results"].append(analyze_taxon(taxon,animals[taxon]))
                animals[taxon].clear()
            output["source_status"]="SOURCE_CADENCE_SENSITIVITY_COMPLETED"
        except Exception as ex:
            output["source_status"]="ORIGINAL_SOURCE_OR_CADENCE_HOLD"
            output["source_hold_error_type"]=type(ex).__name__
        finally:
            local.unlink(missing_ok=True)
    if not output["native_taxon_results"]:
        output["native_taxon_results"]=[{"taxon":t,"source_status":"SOURCE_HOLD",
                                        "source_original_even_odd_2h_phase_sensitivity":None}
                                        for t in taxa]
    out_dir.mkdir(parents=True,exist_ok=True)
    dest=out_dir/("original_study_"+str(index)+"_cadence_v1.json")
    dest.write_text(json.dumps(output,sort_keys=True,indent=2,ensure_ascii=False)+"\n",
                    encoding="utf-8")
    print("RSDS_CADENCE_SOURCE_RESULT",json.dumps({
        "index":index,
        "source_status":output["source_status"],
        "original_q_replayed":output["baseline_1h_q_and_daylength_source_effects_replayed"],
        "native_taxa":len(output["native_taxon_results"]),
        "q_identifiable_by_phase":{
            phase:sum(t.get("source_original_even_odd_2h_phase_sensitivity") is not None
            and t["source_original_even_odd_2h_phase_sensitivity"][phase][
            "matched_q_velocity_proxy_difference"]["status"]=="SOURCE_TAGGED_ANIMAL_DESCRIPTIVE"
            for t in output["native_taxon_results"]) for phase in PAIRS
        }
    },sort_keys=True),flush=True)
    return output


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--study-index",type=int,required=True)
    cli.add_argument("--out",type=Path,required=True)
    args=cli.parse_args()
    if not 0<=args.study_index<8:
        cli.error("Exactly eight previously frozen original source studies, index 0..7")
    args.out.mkdir(parents=True,exist_ok=True)
    run_source(args.study_index,args.out)


if __name__=="__main__":
    main()
