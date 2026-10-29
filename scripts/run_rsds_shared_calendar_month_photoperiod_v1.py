"""Two original Movebank coobserved-community photoperiod block-condition audit.

Original Illinois 3 taxa/169 shared source calendar dates and Wyoming 2 taxa/
193 shared dates. Before NEW slopes, verify original public provider source
bitstream UUID/size/MD5, prior original q and daylength source means/tag counts,
unchanged native taxon/date fixedpoint support and OLD unadjusted datewise q
pair slopes to 1e-12.

Outputs only DOI/taxon pair and aggregate numeric support/bootstrap. Original
dates, tagged animals, GPS, weather or per-date scientific effects NOT shared.
This is post-result diagnosis of seasonal month confounding; photoperiod is
deterministic with day-of-year so within-month association is NOT causal.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download, daylength
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon, slope
from scripts.run_rsds_predator_pair_shared_date_temporal_v1 import (
    check_original_source_coverage, private_tag_calendar, tag_effects, summarize_tags
)
from scripts.run_rsds_illinois_three_native_shared_calendar_v1 import (
    freeze as freeze_illinois, triple_common_calendar, TAXA as ILL_TAXA,
    PAIRS as ILL_PAIRS, REFERENCE as ILL_REFERENCE
)
from scripts.run_rsds_wyoming_coobserved_photoperiod_v1 import (
    frozen as freeze_wyoming, replay_original_triple,
    as_three_original_calendar, matched_two_taxon_calendar,
    same_calendar_date_difference,
    TAXA as WYO_TAXA, SCHEDULES as WYO_SCHEDULES,
    REFERENCE_LATITUDE as WYO_REFERENCE
)
from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    source_reference, source_baseline_check, verify_frame
)
from scripts.run_rsds_triple_matched_ecological_q_v1 import git_blob
from scripts.rsds_within_month_photoperiod_shared_calendar_v1 import (
    compare_frozen_datewise_prior, month_conditional_slope,
    interpret_month_adjustment
)

CONTRACT=HERE/"research/RSDS_SHARED_MONTH_ADJUSTED_PHOTOPERIOD_SOURCE_V1_CONTRACT.json"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
RESULT_NAMES={"Illinois":"illinois","Wyoming":"wyoming"}


def frozen(region):
    contract=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract["version"]!="RSDS_WITHIN_MONTH_SHARED_CALENDAR_PHOTOPERIOD_V1":
        raise ValueError("Source-original same-calendar within-month ecology was not frozen")
    refs={x["region"]:x for x in contract["fixed_original_sources"]}
    if set(refs)!=set(RESULT_NAMES):
        raise ValueError("Fixed both original field study DOI units absent")
    src=refs[region]
    raw=(HERE/src["source_result"]).read_bytes()
    if git_blob(raw)!=src["source_result_git_blob_sha"]:
        raise ValueError("Prior same-calendar original published native source result was modified")
    previous=json.loads(raw.decode("utf-8"))
    if src["region"]!=region or previous["status"] not in (
        "THREE_NATIVE_SHARED_SOURCE_DATE_EVALUATED",
        "STUDY_COMMON_DATE_EVALUATED"):
        raise ValueError("Previous native source has no valid shared calendar-day result")
    if (not previous["old_12_ERA5_raw_SHA_HOLDS_UNCHANGED"]
            if region=="Illinois" else not previous["historical_12_ERA5_raw_JSON_HOLDS_UNCHANGED"]):
        raise ValueError("Original ERA5 external weather version HOLD was modified")
    key=("prior_48_negative_shortwave_direction_external_UNCHANGED" if region=="Illinois"
         else "historical_48_source_signed_radiation_external_NONREPLICATED_UNCHANGED")
    if previous[key]!=NEG:
        raise ValueError("Historical independent 48 original source radiation failure changed")
    return contract,src,previous


def datewise_native_q(aligned,dates,taxa,schedule):
    """Only private source tag/year/date records; equal available source tags/day."""
    out={}
    for taxon in taxa:
        bydate=defaultdict(list)
        for tag,years in aligned[taxon].items():
            for year,observations in years.items():
                for day,samples in observations.items():
                    bydate[(year,day)].append(samples[schedule])
        if set(bydate)!=dates:
            raise ValueError("Shared source calendar dates not identical between taxa")
        out[taxon]={d:mean(x) for d,x in bydate.items()}
    return out


def one_pair(date_q,first,second,schedule,prior_slope,prior_sign,seed):
    keys=set(date_q[first])&set(date_q[second])
    if keys!=set(date_q[first]) or keys!=set(date_q[second]):
        raise ValueError("Source same-calendar species not truly co-observed")
    contrast={k:date_q[first][k]-date_q[second][k] for k in keys}
    original=slope([(daylength(43.0,day),v) for (year,day),v in contrast.items()])
    replay=compare_frozen_datewise_prior(original,prior_slope)
    report=month_conditional_slope(contrast,43.,daylength,seed=seed)
    return {
        "native_first":first,"native_second":second,"source_GPS_clock_schedule":schedule,
        "same_study_original_calendar_days":len(keys),
        "previous_all_months_source_daylength_contrast_coefficient":prior_slope,
        "previous_unadjusted_source_nominal_sign":prior_sign,
        "archived_first_unadjusted_original_coefficient_replayed":replay,
        "within_calendar_month_analysis":report,
        "between_vs_within_month_status":interpret_month_adjustment(
            prior_sign, report["nominal_interval_sign"]),
        "source_sample_is_one_original_field_study":True
    }


def illinois(out):
    contract,reg,first=frozen("Illinois")
    spec,prior_source,source_solar=freeze_illinois()
    original=spec["source_original_study"]
    if (original["doi"]!=reg["doi"] or reg["taxa"]!=list(ILL_TAXA)
            or first["source_all_three_native_support"]["shared_source_calendar_dates_all_three_native_taxa"]!=reg["shared_native_calendar_days"]
            or first["source_all_three_native_support"]["shared_source_year_month_blocks"]!=reg["year_month_blocks"]):
        raise ValueError("Illinois prior original study/source taxon calendar changed")
    result={
        "version":"RSDS_ILLINOIS_WITHIN_MONTH_SHARED_CALENDAR_PHOTOPERIOD_V1",
        "region":"Illinois","source_DOI":reg["doi"],
        "native_taxa":list(ILL_TAXA),"status":"ORIGINAL_SOURCE_HOLD",
        "provider_original_file_MD5_replayed":False,
        "prior_native_q_and_daylength_replayed_before_month_diagnosis":False,
        "same_original_calendar_days":None,"source_year_month_blocks":None,
        "source_native_pair_month_contrasts":[],
        "historic_12_ERA5_source_HOLDS_UNCHANGED":True,
        "historic_48_external_sunlight_signed_nonreplication_UNCHANGED":NEG
    }
    tmp=out/"private_month_illinois_original.csv"
    stage="SOURCE_MD5"
    try:
        fs=[f for f in original_files(original["package_uuid"])
            if f["bitstream_uuid"]==original["bitstream_uuid"]]
        if (len(fs)!=1 or fs[0]["bytes"]!=original["bytes"] or
                fs[0]["checksum"]["value"].lower()!=original["md5"]):
            raise ValueError("Original Illinois source MD5/bitstream identity mismatch")
        verified=checksum_and_download(fs[0],tmp)
        if not verified["provider_checksum_verified"]:
            raise ValueError("Original Illinois provider bytes are not MD5-verified")
        result["provider_original_file_MD5_replayed"]=True
        stage="ORIGINAL_NATIVE_Q_SOURCE_REPLAY"
        tracks,coverage,total=load_source_taxa(tmp,original["file_name"],list(ILL_TAXA))
        check_original_source_coverage(tracks,original,source_solar,total)
        for taxon in ILL_TAXA:
            from_original,_=one_taxon(taxon,tracks[taxon])
            prior=prior_source["source_native_taxa"][taxon]
            for result_name,first_name in (
                    ("day_night_speed_contrast","q"),
                    ("daylength_association","q_daylength")):
                d=from_original[result_name]
                ref=prior[first_name]
                if (d is None or ref["mean"] is None or
                        d["n_individuals"]!=ref["source_tags"] or
                        not math.isclose(d["mean"],ref["mean"],abs_tol=1e-12,rel_tol=0)):
                    raise ValueError("First Illinois original all-day q/daylength source not replayed")
        result["prior_native_q_and_daylength_replayed_before_month_diagnosis"]=True
        stage="ORIGINAL_THREE_NATIVE_SHARED_DATES_AND_Q_REPLAY"
        private={sp:private_tag_calendar(tracks[sp],"original") for sp in ILL_TAXA}
        aligned,dates,loops=triple_common_calendar(private)
        nblocks=len({(year,day.month) for year,day in dates})
        if len(dates)!=reg["shared_native_calendar_days"] or nblocks!=reg["year_month_blocks"]:
            raise ValueError("Original Illinois three-native calendar support changed")
        previous=first["source_all_native_q_levels_on_same_original_calendar"]
        for taxon in ILL_TAXA:
            q,light=tag_effects(aligned[taxon],("original",))
            if (len(q)!=previous[taxon]["q_source_tagged_animals"] or
                    len(light)!=previous[taxon]["daylength_source_tagged_animals"]):
                raise ValueError("Original shared-calendar tag population changed")
            if q:
                v=mean(x["original"] for x in q.values())
                if not math.isclose(v,previous[taxon]["q_original"]["source_animal_mean"],abs_tol=1e-12,rel_tol=0):
                    raise ValueError("Original Illinois shared-calendar native q changed")
            if light:
                v=mean(x["original"] for x in light.values())
                if not math.isclose(v,previous[taxon]["daylength_original"]["source_animal_mean"],abs_tol=1e-12,rel_tol=0):
                    raise ValueError("Original Illinois shared-calendar native photoperiod coefficient changed")
        date_q=datewise_native_q(aligned,dates,ILL_TAXA,"original")
        result["same_original_calendar_days"]=len(dates)
        result["source_year_month_blocks"]=nblocks
        stage="WITHIN_MONTH_SOURCE_NATIVE_CONTRAST"
        for i,(a,b) in enumerate(ILL_PAIRS):
            previous=first["three_source_native_pair_q_and_daylength"][f"{a} MINUS {b}"]["daily_mean_q_difference_vs_standard_daylength"]
            v=one_pair(date_q,a,b,"original",previous["slope"],previous["status"],20261029+i)
            result["source_native_pair_month_contrasts"].append(v)
        result["status"]="ORIGINAL_IL_THREE_NATIVE_MONTH_CONTRAST_COMPLETE"
    except Exception as exc:
        result["status"]="ORIGINAL_ILLINOIS_SOURCE_OR_SHARED_MONTH_HOLD"
        result["hold_stage"]=stage
        result["hold_error_type"]=type(exc).__name__
    finally:
        tmp.unlink(missing_ok=True)
    out.mkdir(parents=True,exist_ok=True)
    (out/"illinois_shared_month.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_SHARED_MONTH_ILLINOIS",json.dumps({
        "source_status":result["status"],
        "MD5":result["provider_original_file_MD5_replayed"],
        "old_q_replayed":result["prior_native_q_and_daylength_replayed_before_month_diagnosis"],
        "same_days":result["same_original_calendar_days"],
        "n_month_contrasts":len(result["source_native_pair_month_contrasts"]),
        "within_month_signs":[x["within_calendar_month_analysis"]["nominal_interval_sign"]
                              for x in result["source_native_pair_month_contrasts"]]
    },sort_keys=True),flush=True)
    return result


def wyoming(out):
    contract,reg,first=frozen("Wyoming")
    spec,prior_absolute,prior_delta=freeze_wyoming()
    old,new,_,_=verify_frame()
    source,taxa,reference,_=source_reference("Wyoming",old,new)
    if (reg["doi"]!=source["doi"] or taxa!=list(WYO_TAXA) or
            first["shared_support"]["original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"]!=reg["shared_native_calendar_days"] or
            first["calendar_date_paired_source_difference"]["shared_year_month_blocks"]!=reg["year_month_blocks"]):
        raise ValueError("Original Wyoming study native taxonomy or shared calendar changed")
    result={
        "version":"RSDS_WYOMING_WITHIN_MONTH_SHARED_CALENDAR_PHOTOPERIOD_V1",
        "region":"Wyoming","source_DOI":reg["doi"],
        "native_taxa":list(WYO_TAXA),"status":"ORIGINAL_SOURCE_HOLD",
        "provider_original_file_MD5_replayed":False,
        "prior_native_q_and_daylength_replayed_before_month_diagnosis":False,
        "same_original_calendar_days":None,"source_year_month_blocks":None,
        "source_native_pair_month_contrasts":[],
        "historic_12_ERA5_source_HOLDS_UNCHANGED":True,
        "historic_48_external_sunlight_signed_nonreplication_UNCHANGED":NEG
    }
    tmp=out/"private_month_wyoming_original.csv"
    stage="SOURCE_MD5"
    try:
        fs=[f for f in original_files(source["package_uuid"])
            if f["bitstream_uuid"]==source["original_bitstream_uuid"]]
        if (len(fs)!=1 or fs[0]["bytes"]!=source["original_bytes"] or
                fs[0]["checksum"]["value"].lower()!=source["original_md5"].lower()):
            raise ValueError("Original Wyoming source MD5/bitstream identity mismatch")
        verified=checksum_and_download(fs[0],tmp)
        if not verified["provider_checksum_verified"]:
            raise ValueError("Original Wyoming provider bytes are not MD5-verified")
        result["provider_original_file_MD5_replayed"]=True
        stage="ORIGINAL_NATIVE_Q_SOURCE_REPLAY"
        tracks,coverage,total=streaming_time_selected_csv(tmp,taxa)
        source_baseline_check(tracks,reference)
        replay_original_triple(tracks,prior_absolute,prior_delta)
        result["prior_native_q_and_daylength_replayed_before_month_diagnosis"]=True
        stage="ORIGINAL_SHARED_CALENDAR_REPLAY"
        private={sp:as_three_original_calendar(tracks[sp]) for sp in WYO_TAXA}
        aligned,dates,loops=matched_two_taxon_calendar(private)
        nblocks=len({(year,day.month) for year,day in dates})
        if len(dates)!=reg["shared_native_calendar_days"] or nblocks!=reg["year_month_blocks"]:
            raise ValueError("Original Wyoming two-native shared dates or month support changed")
        result["same_original_calendar_days"]=len(dates)
        result["source_year_month_blocks"]=nblocks
        stage="WITHIN_MONTH_SOURCE_NATIVE_CONTRAST"
        diffs=same_calendar_date_difference(aligned,dates)
        for i,sched in enumerate(WYO_SCHEDULES):
            date_q={d:row[sched] for d,row in diffs.items()}
            previous=first["calendar_date_paired_source_difference"]["slopes"][sched]
            v=one_pair({WYO_TAXA[0]:date_q,
                        WYO_TAXA[1]:{k:0 for k in date_q}},
                       WYO_TAXA[0],WYO_TAXA[1],sched,
                       previous["mule_q_minus_bighorn_q_vs_shared_reference_daylength"],
                       previous["status"],20261029+10+i)
            result["source_native_pair_month_contrasts"].append(v)
        result["status"]="ORIGINAL_WY_TWO_NATIVE_MONTH_CONTRAST_COMPLETE"
    except Exception as exc:
        result["status"]="ORIGINAL_WYOMING_SOURCE_OR_SHARED_MONTH_HOLD"
        result["hold_stage"]=stage
        result["hold_error_type"]=type(exc).__name__
    finally:
        tmp.unlink(missing_ok=True)
    out.mkdir(parents=True,exist_ok=True)
    (out/"wyoming_shared_month.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_SHARED_MONTH_WYOMING",json.dumps({
        "source_status":result["status"],
        "MD5":result["provider_original_file_MD5_replayed"],
        "old_q_replayed":result["prior_native_q_and_daylength_replayed_before_month_diagnosis"],
        "same_days":result["same_original_calendar_days"],
        "n_month_contrasts":len(result["source_native_pair_month_contrasts"]),
        "within_month_signs":[x["within_calendar_month_analysis"]["nominal_interval_sign"]
                              for x in result["source_native_pair_month_contrasts"]]
    },sort_keys=True),flush=True)
    return result


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--source",choices=("Illinois","Wyoming"),required=True)
    cli.add_argument("--out",type=Path,required=True)
    a=cli.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    {"Illinois":illinois,"Wyoming":wyoming}[a.source](a.out)


if __name__=="__main__":
    main()
