"""Integrate FOUR original coobserved GPS field calendars' within-month daylight contrasts.

Utah was explicitly NOT_TESTED in the immutable four-study first result. This
report adds only the separate, provider-MD5 and old 199-day replay-verified
Utah mule-deer minus bighorn month-conditioned datewise contrast, without
rewriting the prior file or treating source-native pairs as independent
global field-study replications.

This source-only integration never opens original GPS tracks, tag IDs, dates
or climate records. It records SHA fingerprints of its inputs and retains
undefined month effects as UNAVAILABLE (not zero).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FOUR=ROOT/"research/RSDS_FOUR_COOOBSERVED_CALENDAR_NICHE_V1_FIRST_RESULT_20261031.json"
UTAH=ROOT/"research/RSDS_UTAH_COOOBSERVED_WITHIN_MONTH_DAYLENGTH_V1_FIRST_RESULT_20261108.json"
OUTPUT=ROOT/"research/RSDS_FOUR_ORIGINAL_SHARED_CALENDAR_WITHIN_MONTH_DAYLENGTH_V1_FIRST_RESULT_20261108.json"
FOUR_SHA="a15bdf0a6edf9e5d632ddd2b46983af7948355b0"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
SOURCES={
    "Illinois":"10.5441/001/1.649",
    "Wyoming":"10.5441/001/1.703",
    "Utah":"10.5441/001/1.711",
    "Geese":"10.5441/001/1.fk899541",
}
STATUS=("POSITIVE","NEGATIVE","INCLUDES_ZERO","UNAVAILABLE")


def blob(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()


def validate_signal(point,ci,sign):
    if sign not in STATUS:
        raise ValueError("Unknown source native pair month daylight result sign")
    if sign=="UNAVAILABLE":
        if point is not None or ci is not None:
            raise ValueError("Unsupported within-calendar-month source result was converted to numeric")
        return
    if point is None or not isinstance(ci,(list,tuple)) or len(ci)!=2:
        raise ValueError("Estimated original source month daylight contrast has no uncertainty")
    low,high=ci
    if not all(isinstance(x,(float,int)) and math.isfinite(x) for x in (low,point,high)):
        raise ValueError("Nonfinite date-block source month daylight contrast")
    if not low<=point<=high:
        raise ValueError("Source month daylight point estimate lies outside nominal block bootstrap")
    actual="POSITIVE" if low>0 else "NEGATIVE" if high<0 else "INCLUDES_ZERO"
    if actual!=sign:
        raise ValueError("Original source-month sign disagrees with its bootstrap interval")


def combine_four(old,utah,old_sha,utah_sha):
    if old_sha!=FOUR_SHA:
        raise ValueError("Original four field calendar source first result Git blob changed")
    previous=old["native_original_month_conditional_pair_source_comparisons"]
    if (old["original_source_studies_with_identical_native_source_calendar_days"]!=4
            or set(previous)!=set(SOURCES)
            or previous["Utah"]!="NOT_TESTED_WITHIN_MONTH_DO_NOT_INFER_ZERO"
            or old["old_48_source_negative_radiation_direction_external_replication_UNCHANGED"]!=NEG
            or not old["old_12_ERA5_source_version_raw_JSON_HOLDS_UNCHANGED"]):
        raise ValueError("Original coobserved four field-study source month frame changed")
    if (utah["version"]!="RSDS_UTAH_TWO_NATIVE_COOOBSERVED_MONTH_CENTERED_DAYLENGTH_V1_FIRST"
            or utah["original_study_doi"]!=SOURCES["Utah"]
            or utah["fixed_two_target_taxa"]!=["Odocoileus hemionus","Ovis canadensis"]
            or set(utah["all_original_native_study_taxa"])!={"Antilocapra americana","Cervus elaphus","Odocoileus hemionus","Ovis canadensis"}
            or not utah["earlier_four_study_Utah_within_month_not_tested_UNCHANGED"]
            or utah["old_external_48_source_negative_radiation_replication_UNCHANGED"]!=NEG
            or not utah["old_ERA5_raw_JSON_12_source_HOLDS_UNCHANGED"]):
        raise ValueError("New Utah source native taxa or previous external/source result changed")
    Utah_ok=utah["status"]=="UTAH_SOURCE_MONTH_CENTERED_ORIGINAL_EVALUATED"
    if Utah_ok and (
        not utah["provider_original_GPS_bitstream_MD5_verified"] or
        not utah["all_four_native_original_q_and_daylength_replayed_before_month_test"] or
        not utah["original_two_taxa_triple_source_GPS_replayed_before_month_test"] or
        utah["strict_shared_source_date_support"]["source_shared_original_calendar_days"]!=199 or
        utah["strict_shared_source_date_support"]["original_source_year_month_blocks"]!=10 or
        len(utah["month_conditioned_original_pair_results"])!=3):
        raise ValueError("Utah within-month effect was published without full prior original source replay")
    rows=[]
    by_study={}
    signs=Counter()
    for study in SOURCES:
        vals=[]
        if study=="Utah":
            if Utah_ok:
                orig=next(x for x in utah["month_conditioned_original_pair_results"]
                          if x["source_GPS_UTC_clock_schedule"]=="original")
                m=orig["month_centered_source_daylength_response"]
                if orig["all_months_prior_source_slope_replay"]!="ARCHIVED_UNADJUSTED_SLOPE_REPLAYED":
                    raise ValueError("Archived Utah all-month species q-daylight slope was not exactly replayed")
                vals=[{
                    "native_taxa":[orig["source_native_first"],orig["source_native_second"]],
                    "month_slope":m["within_month_daylength_slope"],
                    "month_block_bootstrap95":m["month_cluster_bootstrap95"],
                    "nominal_source_month_sign":m["nominal_interval_sign"],
                    "source_original_year_month_blocks":m["source_year_month_blocks"],
                    "source_months_used":m["eligible_source_month_blocks"],
                    "source_calendar_days_used":m["eligible_source_calendar_days"],
                    "unadjusted_prior_pair_daylength_slope":orig["all_months_original_shared_daylength_slope_per_hour"],
                    "new_Utah_3_observation_schedules_available":True,
                }]
            else:
                vals=[{
                    "native_taxa":["Odocoileus hemionus","Ovis canadensis"],
                    "month_slope":None,"month_block_bootstrap95":None,
                    "nominal_source_month_sign":"UNAVAILABLE",
                    "source_original_year_month_blocks":10,
                    "source_months_used":None,"source_calendar_days_used":None,
                    "unadjusted_prior_pair_daylength_slope":None,
                    "source_HOLD_reason":utah.get("hold_stage") or "ORIGINAL_SOURCE_OR_MONTH_SUPPORT_HOLD",
                    "new_Utah_3_observation_schedules_available":False,
                }]
        else:
            for old_pair in previous[study]:
                vals.append({
                    "native_taxa":old_pair["native_taxa"],
                    "month_slope":old_pair["conditional_reference_daylength_slope"],
                    "month_block_bootstrap95":old_pair["year_month_block_bootstrap95"],
                    "nominal_source_month_sign":old_pair["nominal_sign"],
                    "source_original_year_month_blocks":old_pair["source_year_month_blocks"],
                    "source_months_used":old_pair["source_year_month_blocks"],
                    "source_calendar_days_used":None,
                    "unadjusted_prior_pair_daylength_slope":None,
                    "new_Utah_3_observation_schedules_available":False,
                })
        for x in vals:
            validate_signal(x["month_slope"],x["month_block_bootstrap95"],x["nominal_source_month_sign"])
            signs[x["nominal_source_month_sign"]]+=1
            if len(x["native_taxa"])!=2:
                raise ValueError("Original source study must compare exactly two registered taxa")
            rows.append({
                "original_source_field_study":study,
                "original_published_study_doi":SOURCES[study],
                "original_native_pair":x["native_taxa"],
                **{k:v for k,v in x.items() if k!="native_taxa"}
            })
        by_study[study]={
            "original_study_doi":SOURCES[study],
            "native_pair_rows":len(vals),
            "source_has_any_nominal_month_slope":
                any(x["nominal_source_month_sign"] in ("POSITIVE","NEGATIVE") for x in vals),
            "source_signed_month_direction_set":sorted({
                x["nominal_source_month_sign"] for x in vals if x["nominal_source_month_sign"] in ("POSITIVE","NEGATIVE")
            }),
            "original_study_one_unit_not_native_pair_replicates":True,
        }
    if len(rows)!=8 or len(by_study)!=4:
        raise ValueError("Original four-source eight-native-pair month sample frame incomplete")
    return {
        "version":"RSDS_FOUR_ORIGINAL_STUDY_MONTH_CENTERED_NATIVE_DAYLIGHT_CONTRAST_V1_FIRST",
        "date":"2026-11-08",
        "study_design":"POST_FIRST_MONTHLY_SOURCE_RESULT_ORIGINAL_4_FIELD_STUDIES",
        "original_field_studies":4,
        "original_study_native_pair_rows":8,
        "old_original_three_studies_native_pair_rows":7,
        "new_Utah_month_source_pair_rows":1,
        "new_Utah_source_status":utah["status"],
        "source_month_nominal_bootstrap_sign_counts":dict(signs),
        "original_study_equal_sign_support":by_study,
        "all_four_studies_native_pair_month_results":rows,
        "original_old_four_study_git_blob_sha":old_sha,
        "Utah_new_original_MD5_verified_source_result_git_blob_sha":utah_sha,
        "earlier_4_study_first_UTAH_UNTESTED_statement_UNCHANGED":True,
        "original_12_ERA5_raw_source_SHA_HOLDS_UNCHANGED":True,
        "previous_48_source_external_signed_radiation_NONREPLICATED_UNCHANGED":NEG,
        "science_ceiling":[
            "The original Utah month contrast is a separate post-result reanalysis; old four-study result remains originally NOT_TESTED. No previous first-effect values are rewritten.",
            "Four DOIs represent four original field studies, not eight independently sampled native taxon pairs or independent daily GPS events.",
            "Only original within-year×calendar-month species q differences and their month-block descriptive intervals are compared. No numerical mean across differing sources/latitude/time exposure is justified.",
            "Month adjustment cannot separate photoperiod from temperature, forage, predators, migrating behavior or other within-month seasonal co-variation.",
            "All bootstrap95 intervals are exploratory and unadjusted for paired taxa and three source UTC schedules; only original source schedule included in four-study cross-source comparison.",
            "Source q is observed GPS minimum straight-line NIGHT−DAY displacement proxy, not animal true activity-duration or chronotype physiology.",
            "Original source/tag date filters, rare original Utah mule tags, available month-year blocks and GPS recording phase bias limit generalization.",
            "Old ERA5 source raw SHA version holds, external 48 original-source signed radiation nonreplication, and any non-estimable source month effects all remain unchanged.",
            "Original animal tag IDs, exact GPS positions, UTC times, individual effects and actual source calendar dates are never shared."
        ]
    }


def run():
    raw=FOUR.read_bytes()
    s=blob(raw)
    if s!=FOUR_SHA:
        raise ValueError("Original four-source ecological first result changed")
    if not UTAH.is_file():
        raise ValueError("New Utah provider-MD5 first source month result has not been published")
    q=UTAH.read_bytes()
    return combine_four(json.loads(raw.decode("utf-8")),json.loads(q.decode("utf-8")),
                        old_sha=s,utah_sha=blob(q))


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,default=OUTPUT)
    args=cli.parse_args()
    d=run()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_FOUR_ORIGINAL_SOURCE_WITHIN_MONTH",json.dumps({
        "studies":d["original_field_studies"],"pairs":d["original_study_native_pair_rows"],
        "Utah_status":d["new_Utah_source_status"],
        "source_month_signs":d["source_month_nominal_bootstrap_sign_counts"],
        "old_negative_holds":d["previous_48_source_external_signed_radiation_NONREPLICATED_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
