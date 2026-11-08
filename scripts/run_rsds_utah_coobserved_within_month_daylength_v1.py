"""Utah original coobserved mule-deer/bighorn MONTH-conditioned daylight slopes.

One published original Movebank field study (Utah 10.5441/001/1.711) and
one previously fixed pair of native source taxa on EXACTLY 199 coobserved
astronomy-calendar days. Original provider bitstream MD5, all FOUR native
source q/daylength first effects and original three-cadence matched-tag
effects are replayed before month adjustment.

Compare original/even/odd GPS q_mule−q_bighorn datewise contrast to
43°N reference astronomical daylength, after year×month centering using
EXACT fixed thresholds previously used on Illinois/Wyoming.

This is post-result observational seasonal confounding sensitivity,
NOT causal photoperiodic adaptation or an independent field replication.
Only per-study aggregate outputs; no original source GPS, dates or tag IDs.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download, daylength
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_utah_diel_photoperiod_20261019 import slope
from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    source_reference, source_baseline_check, verify_frame, git_blob, NEGATIVE
)
from scripts.run_rsds_utah_coobserved_photoperiod_v1 import (
    frozen as freeze_previous_utah,
    replay_original_triple,
    as_three_original_calendar,
    matched_two_taxon_calendar,
    same_calendar_date_difference,
    TAXA, SCHEDULES, REFERENCE_LATITUDE
)
from scripts.rsds_within_month_photoperiod_shared_calendar_v1 import (
    month_conditional_slope, compare_frozen_datewise_prior,
    interpret_month_adjustment
)

CONTRACT=HERE/"research/RSDS_UTAH_COOOBSERVED_WITHIN_MONTH_DAYLENGTH_V1_CONTRACT.json"
OUTPUT=HERE/"research/RSDS_UTAH_COOOBSERVED_WITHIN_MONTH_DAYLENGTH_V1_FIRST_RESULT_20261108.json"
HISTORIC_SOURCE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
NOMINAL_SEEDS={"original":20261108,"even":20261109,"odd":20261110}


def frozen():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (spec["version"]!="RSDS_UTAH_TWO_NATIVE_SHARED_CALENDAR_WITHIN_MONTH_DAYLENGTH_V1"
            or spec["predeclared_target_taxa"]!=list(TAXA)
            or set(spec["original_schedules"])!=set(SCHEDULES)
            or spec["original_schedules"]!=list(SCHEDULES)
            or spec["astronomical_exposure"].find("43.0 N")<0):
        raise ValueError("Fixed Utah within-calendar-month source ecology contract changed")
    original_ref,previous_absolute,previous_delta=freeze_previous_utah()
    if (original_ref["original_study"]["doi"]!="10.5441/001/1.711"
            or len(original_ref["original_study"]["native_taxa_in_full_original_source"])!=4
            or len(spec["preexisting_reference"])!=3):
        raise ValueError("Original Utah four source taxa and coobserved target sample not fixed")
    sources={}
    for rec in spec["preexisting_reference"]:
        raw=(HERE/rec["path"]).read_bytes()
        if git_blob(raw)!=rec["git_blob_sha"]:
            raise ValueError("Prior Utah/four-study/monthly source result was replaced: "+rec["path"])
        sources[rec["path"]]=json.loads(raw.decode("utf-8"))
    prior=sources["research/RSDS_UTAH_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json"]
    four=sources["research/RSDS_FOUR_COOOBSERVED_CALENDAR_NICHE_V1_FIRST_RESULT_20261031.json"]
    month=sources["research/RSDS_SHARED_MONTH_PHOTOPERIOD_VALIDATED_ECOLOGY_20261029.json"]
    shared=prior["shared_support"]
    if (prior["status"]!="STUDY_COMMON_DATE_EVALUATED"
            or prior["original_provider_GPS_MD5_verified"] is not True
            or not prior["both_native_taxa_original_all_date_and_original_triple_source_estimates_replayed"]
            or prior["study_Doi"]!="10.5441/001/1.711"
            or prior["native_taxa"]!=list(TAXA)
            or shared["original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"]!=199
            or shared["date_support_fixed_point_iterations"]!=2
            or shared["original_tags_with_10_days_on_both_taxa_shared_calendar"]!={TAXA[0]:10,TAXA[1]:42}
            or shared["distinct_original_calendar_years_with_coobserved_dates"]!=2
            or shared["common_daylight_standard_reference_latitude_degrees"]!=REFERENCE_LATITUDE
            or prior["historical_48_source_signed_radiation_external_NONREPLICATED_UNCHANGED"]!=HISTORIC_SOURCE
            or prior["historical_12_ERA5_raw_JSON_HOLDS_UNCHANGED"] is not True):
        raise ValueError("Frozen Utah original paired source calendar or signed history changed")
    if (four["original_source_studies_with_identical_native_source_calendar_days"]!=4
            or four["native_original_month_conditional_pair_source_comparisons"]["Utah"]!="NOT_TESTED_WITHIN_MONTH_DO_NOT_INFER_ZERO"
            or four["old_48_source_negative_radiation_direction_external_replication_UNCHANGED"]!=HISTORIC_SOURCE
            or not four["old_12_ERA5_source_version_raw_JSON_HOLDS_UNCHANGED"]
            or month["original_source_field_studies"]!=2
            or month["source_native_taxon_and_schedule_comparisons"]!=6
            or month["old_48_signed_radiation_external_nonreplication_UNCHANGED"]!=HISTORIC_SOURCE):
        raise ValueError("Published four-study month comparison or source old HOLD changed")
    return spec,original_ref,previous_absolute,previous_delta,prior


def audit_dates_and_months(aligned, dates, previous):
    support=previous["shared_support"]
    if set(aligned)!=set(TAXA):
        raise ValueError("Preselected two native original taxa replaced")
    before=support["original_tags_with_10_days_on_both_taxa_shared_calendar"]
    if (len(dates)!=support["original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"]
            or len({(year,d.month) for year,d in dates})!=10
            or len({d.month for year,d in dates})!=support["distinct_calendar_month_numbers_in_coobserved_dates"]
            or len({year for year,d in dates})!=support["distinct_original_calendar_years_with_coobserved_dates"]
            or {t:len(aligned[t]) for t in TAXA}!=before):
        raise ValueError("Original Utah two-taxonomy shared calendar source/date/animal denominator was changed")
    return {
        "source_shared_original_calendar_days":len(dates),
        "original_source_year_month_blocks":len({(year,d.month) for year,d in dates}),
        "distinct_calendar_month_numbers":len({d.month for year,d in dates}),
        "original_calendar_year_count":len({year for year,d in dates}),
        "source_qualified_tagged_animals":before,
        "original_source_same_fixpoint_support_replayed":True,
        "fixed_original_daylight_reference_latitude_degrees":REFERENCE_LATITUDE,
        "original_private_dates_or_tag_ids_exported":False,
    }


def compare_months(daywise, previous, meta):
    """Datewise source-native pair q difference, one study, three fixed clock series."""
    if len(daywise)!=199:
        raise ValueError("Previous Utah original 199-calendar-day data frame was altered")
    old=previous["calendar_date_paired_source_difference"]
    if (old["coobserved_original_source_dates"]!=199
            or old["shared_year_month_blocks"]!=10):
        raise ValueError("Original Utah datewise 10 year-month blocks no longer fixed")
    effects=[]
    for schedule in SCHEDULES:
        values={(year,d):float(pair[schedule]) for (year,d),pair in daywise.items()}
        if len(values)!=len(daywise):
            raise ValueError("Incorrect original datewise sampling or duplicates")
        stored=old["slopes"][schedule]
        y0=stored["mule_q_minus_bighorn_q_vs_shared_reference_daylength"]
        actual=slope([(daylength(REFERENCE_LATITUDE,d),val)
                      for (year,d),val in sorted(values.items())])
        check=compare_frozen_datewise_prior(actual,y0)
        new=month_conditional_slope(values,REFERENCE_LATITUDE,daylength,
                                    seed=NOMINAL_SEEDS[schedule])
        if (new["source_calendar_dates"]!=meta["source_shared_original_calendar_days"]
                or new["source_year_month_blocks"]!=meta["original_source_year_month_blocks"]
                or new["eligible_source_month_blocks"]>meta["original_source_year_month_blocks"]):
            raise ValueError("Month-centered source observation support invalid")
        effects.append({
            "original_study_doi":"10.5441/001/1.711",
            "source_native_first":TAXA[0],
            "source_native_second":TAXA[1],
            "source_GPS_UTC_clock_schedule":schedule,
            "all_months_original_shared_daylength_slope_per_hour":y0,
            "all_months_original_source_nominal_status":stored["status"],
            "all_months_prior_source_slope_replay":check,
            "month_centered_source_daylength_response":new,
            "interpretation_vs_archived_all_months":interpret_month_adjustment(
                stored["status"],new["nominal_interval_sign"]),
            "one_original_study_not_independent_taxon_replicates":True,
        })
    return effects


def run(out):
    spec,original_ref,prior_absolute,prior_delta,previous=frozen()
    old,new,_,_=verify_frame()
    source,taxa,reference,_=source_reference("Utah",old,new)
    result={
        "version":"RSDS_UTAH_TWO_NATIVE_COOOBSERVED_MONTH_CENTERED_DAYLENGTH_V1_FIRST",
        "date":"2026-11-08",
        "original_study_doi":"10.5441/001/1.711",
        "all_original_native_study_taxa":list(taxa),
        "fixed_two_target_taxa":list(TAXA),
        "status":"ORIGINAL_SOURCE_OR_MONTH_SUPPORT_HOLD",
        "provider_original_GPS_bitstream_MD5_verified":False,
        "all_four_native_original_q_and_daylength_replayed_before_month_test":False,
        "original_two_taxa_triple_source_GPS_replayed_before_month_test":False,
        "strict_shared_source_date_support":None,
        "month_conditioned_original_pair_results":[],
        "earlier_four_study_Utah_within_month_not_tested_UNCHANGED":True,
        "old_external_48_source_negative_radiation_replication_UNCHANGED":HISTORIC_SOURCE,
        "old_ERA5_raw_JSON_12_source_HOLDS_UNCHANGED":True,
        "original_shared_month_limit":"POSTHOC_SEASONAL_MONTH_CONFOUNDING_NONCAUSAL"
    }
    out.mkdir(parents=True,exist_ok=True)
    tmp=out/"private_original_Utah_month_centered_two_taxa.zip"
    stage="FROZEN_SOURCE_MD5"
    try:
        files=[f for f in original_files(source["package_uuid"])
               if f["bitstream_uuid"]==source["original_bitstream_uuid"]]
        if (len(files)!=1 or files[0]["bytes"]!=source["original_bytes"]
                or files[0]["checksum"]["value"].lower()!=source["original_md5"].lower()):
            raise ValueError("Original public Utah source file UUID/size/MD5 changed")
        verified=checksum_and_download(files[0],tmp)
        if verified["provider_checksum_verified"] is not True:
            raise ValueError("Original Utah source provider MD5 was not reproduced")
        result["provider_original_GPS_bitstream_MD5_verified"]=True
        stage="ALL_FOUR_NATIVE_SOURCE_Q_AND_DAYLENGTH_REPLAY"
        tracks,counts,unknown=extract_time_only_original_source(tmp,taxa)
        source_baseline_check(tracks,reference)
        result["all_four_native_original_q_and_daylength_replayed_before_month_test"]=True
        stage="TRIPLE_GPS_SOURCE_Q_REPLAY"
        replay_original_triple(tracks,prior_absolute,prior_delta)
        result["original_two_taxa_triple_source_GPS_replayed_before_month_test"]=True
        print("RSDS_UTAH_WITHIN_MONTH_SOURCE_PREFLIGHT",json.dumps({
            "source_DOI":"10.5441/001/1.711",
            "original_provider_MD5":True,
            "four_native_q_daylength_replay":True,
            "two_original_tag_triple_cadence_replay":True
        },sort_keys=True),flush=True)
        stage="ORIGINAL_GPS_SHARED_CALENDAR_SUPPORT"
        raw={taxon:as_three_original_calendar(tracks[taxon]) for taxon in TAXA}
        aligned,dates,iterations=matched_two_taxon_calendar(raw)
        if iterations!=previous["shared_support"]["date_support_fixed_point_iterations"]:
            raise ValueError("Published original Utah common-date fixed-point iterations changed")
        meta=audit_dates_and_months(aligned,dates,previous)
        result["strict_shared_source_date_support"]=meta
        stage="MONTHLY_SOURCE_FIXEDPAIR_RESPONSE"
        daywise=same_calendar_date_difference(aligned,dates)
        result["month_conditioned_original_pair_results"]=compare_months(
            daywise,previous,meta)
        result["status"]="UTAH_SOURCE_MONTH_CENTERED_ORIGINAL_EVALUATED"
        for taxon in TAXA: tracks[taxon].clear()
    except Exception as exc:
        result["status"]="ORIGINAL_SOURCE_OR_MONTH_SUPPORT_HOLD"
        result["hold_stage"]=stage
        result["hold_error_type"]=type(exc).__name__
    finally:
        tmp.unlink(missing_ok=True)
    dest=out/"utah_shared_calendar_month_centered_v1.json"
    dest.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_UTAH_MONTH_CENTERED_ORIGINAL_SOURCE",json.dumps({
        "status":result["status"],
        "provider_MD5":result["provider_original_GPS_bitstream_MD5_verified"],
        "first_four_native_effects_replayed":
            result["all_four_native_original_q_and_daylength_replayed_before_month_test"],
        "source_original_common_dates":(result["strict_shared_source_date_support"] or {}).get(
            "source_shared_original_calendar_days"),
        "source_year_month_blocks":(result["strict_shared_source_date_support"] or {}).get(
            "original_source_year_month_blocks"),
        "within_month_nominal_signs":[x["month_centered_source_daylength_response"][
            "nominal_interval_sign"] for x in result["month_conditioned_original_pair_results"]],
        "old_48_direction_negative_UNCHANGED":True
    },sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.out)


if __name__=="__main__":
    main()
