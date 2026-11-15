"""Paired source-native tag/date month-composition audit on five frozen studies.

Each source tag-year compares unadjusted q~L and year-month-FE q~L on the
EXACT SAME previously month-eligible original tag, month, and date sample.
Original source first full-year q and previous month-FE outputs are replayed
before the new paired (month-adjusted - unadjusted) contrast is computed.
Only source-group summaries can leave ephemeral Actions workers.
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

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned
from scripts.movebank_native_ecology_20261019 import (
    original_files, checksum_and_download, daylength
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.run_rsds_eight_study_gps_cadence_v1 import per_tag_year_days
from scripts.run_rsds_four_new_multinative_original_gps_temporal_v1 import source_reference_latitude
from scripts import rsds_five_source_native_tag_month_conditioned_v1 as old

CONTRACT=HERE/"research/RSDS_FIVE_SOURCE_MATCHED_DATE_MONTH_DECOMPOSITION_CONTRACT_20261115.json"
VERSION="RSDS_FIVE_SOURCE_NATIVE_MATCHED_DATE_MONTH_DECOMPOSITION_V1"
FIRST_RESULT_VERSION="RSDS_FIVE_SOURCE_NATIVE_SAME_TAG_DATE_MATCHED_UNADJUSTED_MONTH_V1"
NBOOT=2000
TOLERANCE=1e-10


def fixed():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    required={
        "version":VERSION,
        "previous_first_source_action_run_id":43722240065,
        "original_same_source_native_study_slots":5,
        "original_month_qualified_sources":5,
        "original_year_month_date_opportunity_min_days":old.MIN_MONTH_DATES,
        "original_year_month_date_opportunity_min_daylength_range_h":old.MIN_MONTH_DAYLENGTH_SPAN,
        "original_tag_year_min_qualified_months":old.MIN_TAG_YEAR_MONTH_BLOCKS,
        "original_tag_year_min_dates":old.MIN_TAG_YEAR_DATES,
        "original_tag_year_min_within_month_daylength_sxx_h2":old.MIN_TAG_YEAR_SXX,
        "min_original_independent_tags_per_source":old.MIN_TAGS,
        "original_bootstrap_replicates":NBOOT,
        "baseline_replay_required_before_new_month_support_effect":True,
        "original_source_MD5_required_before_GPS":True,
        "no_new_source_recruitment_or_identity_replacement":True,
        "historical_48_source_irradiance_heading_external_not_replicated":True,
        "historical_12_ERA5_original_raw_SHA_HOLD":True
    }
    if any(c.get(k)!=v for k,v in required.items()):
        raise ValueError("Frozen within-source matched-support contract modified")
    if (c.get("prior_five_source_month_script_git_blob") !=
        "82259963739e95d39f0f4b5f9ddc677b4a041586"):
        raise ValueError("Earlier source month method algorithm version changed")
    pinned(c["prior_five_source_month_script_path"],c["prior_five_source_month_script_git_blob"])
    pinned(c["prior_five_source_original_month_contract_path"],
           c["prior_five_source_original_month_contract_git_blob"])
    previous=pinned(c["prior_five_source_month_first_result_path"],
                    c["prior_five_source_month_first_result_git_blob"])
    original=old.frozen()
    if (previous.get("version")!="RSDS_FIVE_ORIGINAL_SOURCE_NATIVE_TAG_YEAR_MONTH_META_V1" or
        previous.get("source_status_count")!={"MONTH_NATIVE_ESTIMABLE":5} or
        len(previous["all_five_frozen_original_source_native_analyses_including_HOLD"])!=5 or
        len(original["source_manifest"])!=5):
        raise ValueError("Previous FIVE source first month results or original identities changed")
    for s,r in zip(original["source_manifest"],previous["all_five_frozen_original_source_native_analyses_including_HOLD"]):
        if (s["source_key"]!=r["source_key"] or s["doi"]!=r["original_DOI"] or
            s["taxon"]!=r["source_native_species"] or
            r["source_first_result_baseline_replayed"] is not True or
            r["source_original_MD5_verified"] is not True or
            r["source_gps_status"]!="MONTH_NATIVE_ESTIMABLE"):
            raise ValueError("Previous native month result/frozen source roster changed")
    return c,original,previous["all_five_frozen_original_source_native_analyses_including_HOLD"]


def calc_slope(x,y):
    a=mean(x);b=mean(y)
    xx=sum((t-a)**2 for t in x)
    if xx<=0:raise ValueError("Matched raw slope has no observational daylength range")
    return sum((t-a)*(z-b) for t,z in zip(x,y))/xx


def original_tag_year_selected_same_days(day_rows:dict,reference_lat:float):
    """Same date opportunity as source month FE in previous study, never q-sign selected."""
    prior=old.group_month_stats(day_rows,reference_lat)
    if not prior["supported"]:
        return None
    grouped=defaultdict(list)
    for d,(q,Llocal) in day_rows.items():
        grouped[d.month].append((float(q),float(Llocal),float(daylength(reference_lat,d))))
    valid=[]
    for month,records in sorted(grouped.items()):
        if len(records)<old.MIN_MONTH_DATES:continue
        if (max(z[1] for z in records)-min(z[1] for z in records)<old.MIN_MONTH_DAYLENGTH_SPAN
            or max(z[2] for z in records)-min(z[2] for z in records)<old.MIN_MONTH_DAYLENGTH_SPAN):
            continue
        valid.extend(records)
    if len(valid)!=prior["eligible_days"]:
        raise ValueError("Month-selected source native original calendar dates changed")
    if len(valid)<old.MIN_TAG_YEAR_DATES:
        raise ValueError("Qualified source tag year day support changed")
    q=[z[0] for z in valid]
    loc=[z[1] for z in valid]
    ref=[z[2] for z in valid]
    return {
        "same_original_dates":len(valid),
        "local_raw_same_dates":calc_slope(loc,q),
        "local_month_FE_same_dates":prior["local_slope"],
        "reference_raw_same_dates":calc_slope(ref,q),
        "reference_month_FE_same_dates":prior["reference_slope"],
    }


def source_matched_tag_effect(original_tag_year_data,reference_lat,seed):
    values={}
    selected_years=0
    tagged_days=0
    for tag,yrs in original_tag_year_data.items():
        yy=[]
        for yr,ds in yrs.items():
            if len(ds)<10:continue
            x=original_tag_year_selected_same_days(ds,reference_lat)
            if x is not None:
                yy.append(x)
                selected_years+=1
                tagged_days+=x["same_original_dates"]
        if yy:
            values[tag]={
                key:mean(y[key] for y in yy)
                for key in ("local_raw_same_dates","local_month_FE_same_dates",
                            "reference_raw_same_dates","reference_month_FE_same_dates")
            }
            values[tag]["local_month_minus_raw"]=(
                values[tag]["local_month_FE_same_dates"]-values[tag]["local_raw_same_dates"])
            values[tag]["reference_month_minus_raw"]=(
                values[tag]["reference_month_FE_same_dates"]-values[tag]["reference_raw_same_dates"])
    n=len(values)
    names=sorted(values)
    keys=("local_raw_same_dates","local_month_FE_same_dates","local_month_minus_raw",
          "reference_raw_same_dates","reference_month_FE_same_dates",
          "reference_month_minus_raw")
    result={
        "source_original_eligible_tagged_animals":n,
        "source_original_selected_tag_years":selected_years,
        "source_original_selected_tag_date_records_not_independent":tagged_days,
        "source_tag_bootstrap_replicates":NBOOT if n>=old.MIN_TAGS else 0,
        "all_contrasts_use_identical_original_tag_year_month_date_support":True,
        "source_q_daylength_effects":{}
    }
    if n<old.MIN_TAGS:
        for key in keys:
            result["source_q_daylength_effects"][key]={
                "tag_n":n,"status":"INSUFFICIENT_NATIVE_TAG_SUPPORT",
                "mean":None,"tag_bootstrap95":None,"nominal_sign":"UNAVAILABLE"}
        return result
    rng=random.Random(seed)
    samples=[[names[rng.randrange(n)] for _ in range(n)] for _ in range(NBOOT)]
    for key in keys:
        coeff=[values[tag][key] for tag in names]
        if not all(math.isfinite(x) for x in coeff):
            raise ValueError("Nonfinite tag matched source coefficient")
        draw=sorted(mean(values[t][key] for t in sample) for sample in samples)
        ci=[draw[49],draw[1949]]
        result["source_q_daylength_effects"][key]={
            "tag_n":n,"status":"ESTIMABLE","mean":mean(coeff),
            "tag_bootstrap95":ci,
            "nominal_sign":"POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"}
    return result


def source_process(index):
    c,original,previous=fixed()
    if type(index) is not int or not 0<=index<5:raise ValueError("Five source positions only")
    s=original["source_manifest"][index]
    p=previous[index]
    md=s["original_file"]
    result={
        "version":FIRST_RESULT_VERSION,
        "source_index":index,
        "original_source_key":s["source_key"],
        "original_DOI":s["doi"],
        "source_original_study_id":str(s["study_id"]),
        "source_native_taxon":s["taxon"],
        "all_source_registered_taxa":s["registered_taxa"],
        "old_q_and_daylength_original_baseline_replayed":False,
        "first_month_FE_identical_sample_effect_replayed":False,
        "original_file_UUID_MD5_bytes_verified":False,
        "original_ephemeral_GPS_source_rows_read":0,
        "source_status":"SOURCE_HOLD",
        "source_hold_stage":None,
        "source_hold_error_type":None,
        "matched_month_effect":None,
        "historical_48_solar_heading_negative_external_verdict_unchanged":True,
        "historical_12_ERA5_raw_SHA_HOLD_unchanged":True,
        "source_original_GPS_IDs_coordinates_exact_dates_and_individual_effects_never_exported":True,
    }
    ext=".zip" if md["filename"].lower().endswith(".zip") else ".csv.gz" if md["filename"].lower().endswith(".csv.gz") else ".csv"
    stage="ORIGINAL_PROVIDER_FROZEN_SOURCE_UUID_BYTES_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_match_month_private_") as td:
            path=Path(td)/("source"+ext)
            ff=[x for x in original_files(s["package_uuid"])
                if x["bitstream_uuid"]==md["uuid"]]
            if (len(ff)!=1 or ff[0]["bytes"]!=md["bytes"] or
                ff[0]["checksum"]["checkSumAlgorithm"].upper()!="MD5" or
                ff[0]["checksum"]["value"].lower()!=md["md5"]):
                raise ValueError("Original provider source MD5 differs from frozen file")
            if not checksum_and_download(ff[0],path)["provider_checksum_verified"]:
                raise ValueError("Provider MD5 not checked")
            result["original_file_UUID_MD5_bytes_verified"]=True
            stage="ORIGINAL_NATIVE_GPS_PARSE_AND_Q_REPLAY"
            if s["source_key"]=="mule_old_utah":
                tracks,cnts,unmatched=extract_time_only_original_source(path,s["registered_taxa"])
                original_rows=sum(x["rows"] for x in cnts.values())+unmatched
            else:
                tracks,coverage,report=load_source_taxa(path,md["filename"],s["registered_taxa"])
                original_rows=report["original_total_rows"]
            result["original_ephemeral_GPS_source_rows_read"]=original_rows
            old.replay_original_baseline(s,tracks,TOLERANCE)
            result["old_q_and_daylength_original_baseline_replayed"]=True
            stage="PREVIOUS_ORIGINAL_NATIVE_MONTH_FIRST_RESULT_REPLAY"
            data,_=per_tag_year_days(tracks[s["taxon"]],"ORIGINAL_1H")
            lat=source_reference_latitude(tracks)
            if lat is None:raise ValueError("No source median tagged-latitude for geometry")
            prior=old.same_source_month_effect(data,lat,seed=old.FIRST_SEED+1000*index)
            prev=p["month_native"]
            if (prior["original_tagged_animals_with_year_month_conditions"]!=
                prev["original_tagged_animals_with_year_month_conditions"]):
                raise ValueError("Previous month source tag support changed")
            for key in ("local_original_solar_daylength",
                        "fixed_source_reference_latitude_daylength"):
                for val in ("mean_slope_per_hour",):
                    if not math.isclose(prior[key][val],prev[key][val],rel_tol=0,abs_tol=1e-10):
                        raise ValueError("First source-native month-FE coefficient not replayed")
                if (prior[key]["nominal_interval_sign"]!=prev[key]["nominal_interval_sign"] or
                    prior[key]["tagged_animals_support"]!=prev[key]["tagged_animals_support"]):
                    raise ValueError("First source month bootstrap sign/tag count changed")
            result["first_month_FE_identical_sample_effect_replayed"]=True
            stage="NEW_SOURCE_NATIVE_MATCHED_TAG_DATE_MONTH_DECOMPOSITION"
            new=source_matched_tag_effect(data,lat,seed=20261115+4000*index)
            if new["source_original_eligible_tagged_animals"]!=prior["original_tagged_animals_with_year_month_conditions"]:
                raise ValueError("Matched tag count differs from already frozen month-FE source")
            for key,old_key in (
                ("local_month_FE_same_dates","local_original_solar_daylength"),
                ("reference_month_FE_same_dates","fixed_source_reference_latitude_daylength")
            ):
                if not math.isclose(
                    new["source_q_daylength_effects"][key]["mean"],
                    prior[old_key]["mean_slope_per_hour"],abs_tol=1e-10,rel_tol=0):
                    raise ValueError("Tag/date matched FE source slope differs from previous first source")
            result["matched_month_effect"]=new
            result["source_status"]="MATCHED_SUPPORT_PAIRED_SOURCE_ESTIMABLE" if new["source_original_eligible_tagged_animals"]>=5 else "MATCHED_SUPPORT_NO_NATIVE_TAGS"
            tracks.clear()
    except Exception as e:
        result["source_status"]="SOURCE_HOLD"
        result["source_hold_stage"]=stage
        result["source_hold_error_type"]=type(e).__name__
        result["matched_month_effect"]=None
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--index",required=True,type=int)
    parser.add_argument("--out",required=True,type=Path)
    a=parser.parse_args()
    d=source_process(a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_MATCHED_MONTH_Q",json.dumps({
        "index":d["source_index"],"status":d["source_status"],
        "q_replayed":d["old_q_and_daylength_original_baseline_replayed"],
        "prior_month_replayed":d["first_month_FE_identical_sample_effect_replayed"],
        "month_minus_same_sample_raw":((d["matched_month_effect"] or {}).get(
            "source_q_daylength_effects") or {}).get("local_month_minus_raw"),
        "hold_type":d["source_hold_error_type"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
