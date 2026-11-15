"""Source-native within-tag-year-month q~daylength on five frozen Movebank sources.

This is a post-outcome seasonal confounding SENSITIVITY, not causal daylight
biology. All original provider bitstream UUID/MD5/bytes, earlier source-native
full-year tag means/n, and source identities are pinned before month analysis.
Raw GPS, tag IDs, dates, coordinates and individual results stay ephemeral.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import random
from statistics import mean
import tempfile

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned
from scripts.movebank_native_ecology_20261019 import (
    original_files, checksum_and_download, daylength
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.run_rsds_eight_study_gps_cadence_v1 import per_tag_year_days
from scripts.run_rsds_four_new_multinative_original_gps_temporal_v1 import source_reference_latitude

CONTRACT=HERE/"research/RSDS_FIVE_SOURCE_NATIVE_YEAR_MONTH_ADJUSTMENT_CONTRACT_20261115.json"
VERSION="RSDS_FIVE_SOURCE_NATIVE_WITHIN_TAG_YEAR_MONTH_FIRST_V1"
SOURCE_VERSION="RSDS_FIVE_SOURCE_NATIVE_YEAR_MONTH_CONDITIONED_PER_STUDY_V1"
FIRST_SEED=20261115
NBOOT=2000
MIN_MONTH_DATES=5
MIN_MONTH_DAYLENGTH_SPAN=.05
MIN_TAG_YEAR_MONTH_BLOCKS=2
MIN_TAG_YEAR_DATES=10
MIN_TAG_YEAR_SXX=.05
MIN_TAGS=5


def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (c.get("version")!=VERSION or
        c.get("registered_sources")!=5 or
        c.get("registered_species")!=["Anser albifrons","Odocoileus hemionus"] or
        c.get("max_gps_rows_per_source")!=2_000_000 or
        c.get("file_max_bytes")!=200_000_000 or
        c.get("baseline_replay_tolerance")!=1e-10 or
        c.get("legacy_48_sunlight_heading_external_replicated") is not False or
        c.get("legacy_12_ERA5_raw_SHA_HOLD") is not True or
        len(c.get("source_manifest") or [])!=5 or
        len(c.get("evidence_blob_pins") or [])!=6):
        raise ValueError("Month-sensitive original five-source contract was changed")
    for pinned_file in c["evidence_blob_pins"]:
        pinned(pinned_file["path"],pinned_file["git_blob_sha"])
    catalog=pinned(c["original_movebank_catalog"]["path"],
                   c["original_movebank_catalog"]["git_blob_sha"])
    if catalog.get("package_count")!=437 or len(catalog["packages"])!=437:
        raise ValueError("Original 437 package catalog changed")
    seen=set()
    for z in c["source_manifest"]:
        key=(z["doi"],str(z["study_id"]))
        if key in seen:raise ValueError("Original duplicate study/DOI")
        seen.add(key)
        candidates=[p for p in catalog["packages"] if p["package_uuid"]==z["package_uuid"]]
        if (len(candidates)!=1 or
            candidates[0]["repository_study_ids"]!=[str(z["study_id"])] or
            not any(d.strip().lower().replace("doi:","")==z["doi"].lower() for d in candidates[0]["dois"]) or
            candidates[0]["scientific_names"]!=z["registered_taxa"] or
            z["taxon"] not in z["registered_taxa"]):
            raise ValueError("Frozen source native title/ID/taxa differs from 437 catalog")
        md=z["original_file"]
        if (not 0<md["bytes"]<=c["file_max_bytes"] or
            len(md["md5"])!=32 or
            any(k not in "0123456789abcdef" for k in md["md5"].lower()) or
            not md["uuid"] or not md["filename"].lower().endswith((".csv",".csv.gz",".zip"))):
            raise ValueError("Frozen source original UUID MD5, size or filename invalid")
        for field in ("q","daylength"):
            e=z["baseline"][field]
            if (e["available"] is not True or
                not isinstance(e["original_source_tagged_animals"],int) or
                e["original_source_tagged_animals"]<5 or
                not isinstance(e["q"],(int,float)) or
                not math.isfinite(e["q"])):
                raise ValueError("Original source q/daylength baseline not prior-supported")
    if [z["source_key"] for z in c["source_manifest"]]!=[
        "goose_old_grazing","mule_old_utah","mule_old_wyoming",
        "goose_new_pannonic","mule_new_sublette"]:
        raise ValueError("Five fixed original source identities were reordered")
    return c


def group_month_stats(year_dates:dict, latitude:float):
    """Use ONLY the tag's original q-admitted calendar dates to screen blocks.

    Both local and fixed source latitude daylength variants share exactly the
    same admissible (tag,year,month,date) opportunity. Never filter on slope/q.
    """
    if len(year_dates)<10:
        return {"supported":False,"eligible_months":0,"eligible_days":0,
                "local_sxx":0.,"reference_sxx":0.,
                "local_slope":None,"reference_slope":None}
    grouped=defaultdict(list)
    for d,(q,Llocal) in year_dates.items():
        if not (math.isfinite(q) and math.isfinite(Llocal)):
            raise ValueError("Invalid original source native day q/daylength")
        Lfixed=daylength(latitude,d)
        if not math.isfinite(Lfixed):raise ValueError("Invalid fixed-lat original astronomy")
        grouped[d.month].append((float(q),float(Llocal),float(Lfixed)))
    usable=[]
    for month,records in sorted(grouped.items()):
        if len(records)<MIN_MONTH_DATES:continue
        if (max(x[1] for x in records)-min(x[1] for x in records)
             <MIN_MONTH_DAYLENGTH_SPAN or
            max(x[2] for x in records)-min(x[2] for x in records)
             <MIN_MONTH_DAYLENGTH_SPAN):
            continue
        bars=[mean(z[k] for z in records) for k in range(3)]
        sxxloc=sum((z[1]-bars[1])**2 for z in records)
        sxyloc=sum((z[1]-bars[1])*(z[0]-bars[0]) for z in records)
        sxxref=sum((z[2]-bars[2])**2 for z in records)
        sxyref=sum((z[2]-bars[2])*(z[0]-bars[0]) for z in records)
        usable.append((len(records),sxxloc,sxyloc,sxxref,sxyref))
    count=sum(z[0] for z in usable)
    lx=sum(z[1] for z in usable)
    rx=sum(z[3] for z in usable)
    ok=(len(usable)>=MIN_TAG_YEAR_MONTH_BLOCKS and
        count>=MIN_TAG_YEAR_DATES and
        lx>=MIN_TAG_YEAR_SXX and rx>=MIN_TAG_YEAR_SXX)
    return {"supported":ok,
        "eligible_months":len(usable),
        "eligible_days":count,
        "local_sxx":lx,
        "reference_sxx":rx,
        "local_slope":(sum(z[2] for z in usable)/lx) if ok else None,
        "reference_slope":(sum(z[4] for z in usable)/rx) if ok else None}


def same_source_month_effect(original_tag_year_data:dict,reference_lat:float,seed:int):
    local={};reference={}
    # Aggregate diagnostics count qualified months and tag-years, never dates.
    tag_total_support=0
    selected_years=0
    accepted_months=0
    accepted_days=0
    min_source_sxx=float("inf")
    for tag,years in original_tag_year_data.items():
        yr_local=[];yr_ref=[]
        for yr,ds in years.items():
            if len(ds)<10:continue
            stats=group_month_stats(ds,reference_lat)
            if stats["supported"]:
                yr_local.append(stats["local_slope"])
                yr_ref.append(stats["reference_slope"])
                selected_years+=1
                accepted_months+=stats["eligible_months"]
                accepted_days+=stats["eligible_days"]
                min_source_sxx=min(min_source_sxx,stats["local_sxx"],stats["reference_sxx"])
        if yr_local:
            local[tag]=mean(yr_local)
            reference[tag]=mean(yr_ref)
            tag_total_support+=1
    if set(local)!=set(reference):
        raise ValueError("Local/fixed-lat sample mismatch")
    def fitted(obj,seed):
        arr=list(obj.values())
        out={
            "status":"ESTIMABLE" if len(arr)>=MIN_TAGS else "INSUFFICIENT_NATIVE_TAG_SUPPORT",
            "tagged_animals_support":len(arr),
            "mean_slope_per_hour":None,
            "original_tag_bootstrap95":None,
            "nominal_interval_sign":"UNAVAILABLE",
        }
        if len(arr)<MIN_TAGS:return out
        rng=random.Random(seed)
        n=len(arr)
        draws=sorted(mean(arr[rng.randrange(n)] for i in range(n))
                     for j in range(NBOOT))
        ci=[draws[49],draws[1949]]
        out.update(mean_slope_per_hour=mean(arr),original_tag_bootstrap95=ci,
                   nominal_interval_sign=(
                       "POSITIVE" if ci[0]>0 else
                       "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"))
        return out
    return {
        "original_tagged_animals_with_year_month_conditions":tag_total_support,
        "original_tag_years_with_year_month_conditioned_slope":selected_years,
        "source_accepted_tag_month_blocks":accepted_months,
        "source_accepted_tag_solar_dates_count_not_independent":accepted_days,
        "minimum_admitted_tag_year_centered_daylength_Sxx_hours2":
            min_source_sxx if selected_years else None,
        "local_original_solar_daylength":fitted(local,seed),
        "fixed_source_reference_latitude_daylength":fitted(reference,seed+10000),
        "source_reference_latitude_not_exported":True,
        "original_animal_identifiers_and_daily_records_not_exported":True,
    }


def replay_original_baseline(source,tracks,epsilon):
    baseline=source["baseline"]
    output=one_taxon(source["taxon"],tracks[source["taxon"]])[0]
    for tag,key in (("q","day_night_speed_contrast"),("daylength","daylength_association")):
        expected=baseline[tag]
        observed=output[key]
        if observed is None:
            raise ValueError("Original full-source "+tag+" has become unavailable")
        if (expected["original_source_tagged_animals"]!=observed["n_individuals"] or
            not math.isclose(expected["q"],observed["mean"],rel_tol=0,abs_tol=epsilon)):
            raise ValueError("Original full-source "+tag+" baseline no longer replays")
    return {
        "q_native_tag_n":output["day_night_speed_contrast"]["n_individuals"],
        "q_native_original_mean":output["day_night_speed_contrast"]["mean"],
        "daylength_native_tag_n":output["daylength_association"]["n_individuals"],
        "daylength_native_original_mean":output["daylength_association"]["mean"],
        "original_q_and_daylength_baseline_numeric_replayed_before_month_conditioning":True,
    }


def process(index):
    spec=frozen()
    if type(index) is not int or not 0<=index<5:
        raise ValueError("Five original frozen source study slots only")
    s=spec["source_manifest"][index]
    output={
        "version":SOURCE_VERSION,
        "source_index":index,
        "original_study_DOI":s["doi"],
        "original_study_id":str(s["study_id"]),
        "original_native_target":s["taxon"],
        "all_original_registered_native_taxa":s["registered_taxa"],
        "old_original_source_GPS_access_state":s["source_file_prior_access"],
        "baseline_full_year_first_results_replayed":False,
        "original_file_MD5_bytes_verified":False,
        "status":"SOURCE_HOLD",
        "original_source_hold_stage":None,
        "original_source_hold_type":None,
        "source_original_animal_GPS_rows_read_only_ephemerally":0,
        "month_conditioned_native":None,
        "old_solar_direction_external_replication_still_negative":True,
        "old_original_ERA5_raw_SHA_hold":True,
        "no_raw_animal_GPS_tracks_ID_exact_coordinates_dates_or_tag_slopes_exported":True
    }
    ext=(".zip" if s["original_file"]["filename"].lower().endswith(".zip") else
         ".csv.gz" if s["original_file"]["filename"].lower().endswith(".csv.gz") else ".csv")
    stage="ORIGINAL_PROVIDER_SOURCE_UUID_BYTES_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_month_private_") as td:
            fpath=Path(td)/("original_source"+ext)
            candidates=[
                x for x in original_files(s["package_uuid"])
                if x["bitstream_uuid"]==s["original_file"]["uuid"]
            ]
            md=s["original_file"]
            if (len(candidates)!=1 or candidates[0]["bytes"]!=md["bytes"] or
                candidates[0]["checksum"]["checkSumAlgorithm"].upper()!="MD5" or
                candidates[0]["checksum"]["value"].lower()!=md["md5"].lower()):
                raise ValueError("Original file UUID/byte/MD5 changed")
            if not checksum_and_download(candidates[0],fpath)["provider_checksum_verified"]:
                raise ValueError("Original provider file checksum failure")
            output["original_file_MD5_bytes_verified"]=True
            stage="ORIGINAL_PRIVATE_NATIVE_GPS_AND_FULL_SOURCE_BASELINE_REPLAY"
            tracks,coverage,record=load_source_taxa(
                fpath,md["filename"],s["registered_taxa"]
            )
            output["source_original_animal_GPS_rows_read_only_ephemerally"]=record["original_total_rows"]
            output["original_registered_native_taxa_found"]={
                t:coverage[t]["original_taxon_rows"]>0 for t in s["registered_taxa"]
            }
            if not all(output["original_registered_native_taxa_found"].values()):
                raise ValueError("Original multi-native provider file taxonomy incomplete")
            old=replay_original_baseline(s,tracks,spec["baseline_replay_tolerance"])
            output["baseline_full_year_first_results_replayed"]=True
            output["baseline_original_native_first_effect"]=old
            stage="FROZEN_POSTHOC_WITHIN_TAG_YEAR_MONTH_CONDITIONING"
            original_tag_year_days,support=per_tag_year_days(
                tracks[s["taxon"]],"ORIGINAL_1H"
            )
            lat=source_reference_latitude(tracks)
            if lat is None:raise ValueError("Insufficient original source positions for fixed reference latitude")
            month=same_source_month_effect(
                original_tag_year_days,lat,seed=FIRST_SEED+1000*index
            )
            month["original_same_1hour_selected_source_qualified_days"]=support[
                "animal_calendar_days_with_ge3_day_and_ge3_night_segments"]
            month["original_selected_fix_count"]=support["original_selected_1hour_source_fixes"]
            output["month_conditioned_native"]=month
            output["status"]="MONTH_NATIVE_ESTIMABLE" if month["local_original_solar_daylength"]["status"]=="ESTIMABLE" else "MONTH_NATIVE_UNAVAILABLE_INSUFFICIENT_SOURCE_SUPPORT"
            tracks.clear()
    except Exception as e:
        output["status"]="SOURCE_HOLD"
        output["original_source_hold_stage"]=stage
        output["original_source_hold_type"]=type(e).__name__
        output["month_conditioned_native"]=None
    return output


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--index",type=int,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    r=process(a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_NATIVE_MONTH_ONLY",json.dumps({
        "source":r["source_index"],
        "status":r["status"],
        "q_replayed":r["baseline_full_year_first_results_replayed"],
        "month_native_tag_n":(r["month_conditioned_native"] or {}).get(
            "original_tagged_animals_with_year_month_conditions"),
        "month_slope_sign":(r["month_conditioned_native"] or {}).get(
            "local_original_solar_daylength",{}).get("nominal_interval_sign"),
        "hold_stage":r["original_source_hold_stage"],
        "hold_type":r["original_source_hold_type"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
