"""Two original Movebank sources, SAME mule deer/bighorn shared-date q~L contrast.

Old Wyoming result and NEW source-fixed Utah outcome are evaluated at an
identical predeclared 43 N astronomical daylength coordinate, with the SAME
three source GPS observation-clock schedules in each source. Original Utah
outcome is pinned by a separate first-result Git blob receipt before this
comparator runs. Zero/unsupported cases are never counted as nonsignificant
reproductions; the decision to test Utah followed opening Wyoming.

This checks at most two purposively selected ORIGINAL field studies, not a
phylogenetic/general-world parameter and not causal irradiance selection.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
WYOMING="research/RSDS_WYOMING_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json"
UTAH="research/RSDS_UTAH_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json"
PINS="research/RSDS_UTAH_ORIGINAL_SHARED_DATE_FIRST_BLOB_V1.json"
WYOMING_SHA="9b04b18ee73227b4e7a1715ebfb97f610af66952"
SOURCE=("Odocoileus hemionus","Ovis canadensis")
SCHEDULES=("original","even","odd")
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def load_blob(path,sha):
    buf=(ROOT/path).read_bytes()
    actual=hashlib.sha1(("blob "+str(len(buf))+"\0").encode()+buf).hexdigest()
    if actual!=sha:raise ValueError("A frozen original Movebank first-result blob changed: "+path)
    return json.loads(buf.decode("utf-8"))


def validate_study(value,doi,region):
    if value["study_Doi"]!=doi or value["native_taxa"]!=list(SOURCE):
        raise ValueError("Wrong original study DOI or fixed two native taxa: "+region)
    if (not value["historical_12_ERA5_raw_JSON_HOLDS_UNCHANGED"] or
        value["historical_48_source_signed_radiation_external_NONREPLICATED_UNCHANGED"]!=NEGATIVE):
        raise ValueError("Old original source weather/shortwave failed external endpoint reclassified")
    if (value.get("original_provider_GPS_MD5_verified") is not True or
        value.get("both_native_taxa_original_all_date_and_original_triple_source_estimates_replayed") is not True):
        return "ORIGINAL_PROVIDER_OR_SOURCE_REPLAY_HOLD"
    if value["status"]=="STUDY_COMMON_DATE_EVALUATED":
        if len(value["source_taxon_results"])!=2 or not isinstance(value["shared_support"],dict):
            raise ValueError("Source result success with incomplete source native study")
        if (value["shared_support"]["common_daylight_standard_reference_latitude_degrees"]!=43.0 or
            value["shared_support"]["original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"]<10):
            raise ValueError("Source shared-date study changed reference daylength or date denominator")
        if {r["taxon"] for r in value["source_taxon_results"]}!=set(SOURCE):
            raise ValueError("Source first result missing one of two original native taxa")
        return "SHARED_DATE_SOURCE_EVALUATED"
    if value["status"]=="INSUFFICIENT_COOCCURRENT_SOURCE_DATE_SUPPORT":
        return "SHARED_DATE_SOURCE_INSUFFICIENT"
    raise ValueError("Software/provider source HOLD must not be interpreted as biological result")


def component(study,schedule):
    """Native q LEVEL versus q~daylength SLOPE are distinct inference axes."""
    by_taxon={r["taxon"]:r for r in study["source_taxon_results"]}
    native_slopes={}
    native_q={}
    for name in SOURCE:
        record=by_taxon[name]
        for endpoint,destination,label in (
            ("common_date_photoperiod_q_slope_by_schedule",native_slopes,"q~L"),
            ("common_date_temporal_q_by_schedule",native_q,"q")
        ):
            src=record[endpoint]
            if src["status"]!="ESTIMABLE" or src["original_tag_count"]<5:
                return {"status":"INCOMPLETE_ORIGINAL_TAG_SUPPORT","native":None,
                        "contrast":None,"native_q_levels":None,"q_contrast":None}
            val=src["schedules"][schedule]
            ci=val["animal_bootstrap95"]
            if (val["nominal_sign"] not in ("POSITIVE","NEGATIVE","INCLUDES_ZERO") or
                not isinstance(ci,list) or len(ci)!=2 or
                not ci[0]<=val["mean"]<=ci[1]):
                raise ValueError("Tagged-source original "+label+" effect sign/interval invalid")
            destination[name]={"mean":val["mean"],"ci":ci,
                               "sign":val["nominal_sign"],
                               "n_source_original_tags":src["original_tag_count"]}
    contrasts={}
    for sourcekey,native,name in (
        ("q_daylength_slope_mule_minus_bighorn_by_schedule",native_slopes,"slope"),
        ("q_mule_minus_bighorn_by_schedule",native_q,"q")
    ):
        obj=study["same_study_between_native_taxon"][sourcekey]["schedules"][schedule]
        value=obj["mule_minus_bighorn_slope"]
        ci=obj["descriptive_original_tag_bootstrap95"]
        if not isinstance(ci,list) or len(ci)!=2 or not ci[0]<=value<=ci[1]:
            raise ValueError("Original source native two-taxon "+name+" contrast/CI invalid")
        if abs(value-(native[SOURCE[0]]["mean"]-native[SOURCE[1]]["mean"]))>1e-12:
            raise ValueError("Within-source "+name+" contrast fails original native algebra")
        contrasts[name]={"mule_minus_bighorn":value,"ci":ci,
                         "nominal_sign":obj["status"]}
    datewise=(study.get("calendar_date_paired_source_difference") or {}).get("slopes",{}).get(schedule)
    date_slope=None
    if datewise:
        date_slope={
            "source_date_mean_q_difference_vs_reference_daylength":
                datewise["mule_q_minus_bighorn_q_vs_shared_reference_daylength"],
            "year_month_block_descriptive95":datewise["year_month_block_descriptive95"],
            "status":datewise["status"],
        }
    return {
        "status":"EVALUATED",
        "native":native_slopes,
        "contrast":contrasts["slope"],
        "native_q_levels":native_q,
        "q_contrast":contrasts["q"],
        "date_level_q_difference_photoperiod_slope":date_slope,
        "positive_mule_negative_bighorn_supported":(
            native_slopes[SOURCE[0]]["sign"]=="POSITIVE" and
            native_slopes[SOURCE[1]]["sign"]=="NEGATIVE" and
            contrasts["slope"]["nominal_sign"]=="NOMINAL_POSITIVE"),
        "mule_higher_observed_night_relative_q_supported":(
            contrasts["q"]["nominal_sign"]=="NOMINAL_POSITIVE"),
    }


def compare(utah,wyoming):
    statuses={
        "Utah":validate_study(utah,"10.5441/001/1.711","Utah"),
        "Wyoming":validate_study(wyoming,"10.5441/001/1.703","Wyoming")
    }
    records=[]
    for name,st in (("Wyoming",wyoming),("Utah",utah)):
        sc=st.get("shared_support") or {}
        outcome={
            "original_field_study":name,
            "original_study_DOI":st["study_Doi"],
            "shared_calendar_source_status":statuses[name],
            "both_native_taxa_original_tags_on_shared_calendar":
                sc.get("original_tags_with_10_days_on_both_taxa_shared_calendar"),
            "source_study_shared_calendar_dates":
                sc.get("original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"),
            "source_study_year_month_blocks":
                (st.get("calendar_date_paired_source_difference") or {}).get("shared_year_month_blocks"),
            "same_reference_astronomical_daylength_lat":sc.get(
                "common_daylight_standard_reference_latitude_degrees"),
            "three_observation_schedules":{}
        }
        if statuses[name]=="SHARED_DATE_SOURCE_EVALUATED":
            for schedule in SCHEDULES:
                outcome["three_observation_schedules"][schedule]=component(st,schedule)
        else:
            outcome["three_observation_schedules"]={
                s:{"status":"UNAVAILABLE","native":None,"contrast":None}
                for s in SCHEDULES
            }
        records.append(outcome)
    complete=[r for r in records if all(
        r["three_observation_schedules"][s]["status"]=="EVALUATED" for s in SCHEDULES)]
    shared_sign=[
        r["original_field_study"] for r in complete if
        all(r["three_observation_schedules"][s]["positive_mule_negative_bighorn_supported"]
            for s in SCHEDULES)
    ]
    q_order_studies=[
        r["original_field_study"] for r in complete if
        all(r["three_observation_schedules"][s]["mule_higher_observed_night_relative_q_supported"]
            for s in SCHEDULES)
    ]
    return {
        "version":"RSDS_TWO_ORIGINAL_STUDIES_SAME_NATIVE_TAXA_COMMON_CALENDAR_PHOTOPERIOD_V1",
        "date":"2026-10-29",
        "stage":"POST_WYOMING_POSITIVE_RESULT_ADDITIONAL_FIXED_UTAH_SOURCE_DIAGNOSTIC",
        "fixed_original_study_DOIs":2,
        "comparison_native_taxa":list(SOURCE),
        "source_43N_astronomical_reference_shared":True,
        "source_original_studies_full_calendar_comparison_identifiable":len(complete),
        "original_studies_with_same_opposite_native_daylength_sign_through_three_clocks":shared_sign,
        "both_original_studies_show_same_source_native_sign_if_supported":(
            len(complete)==2 and len(shared_sign)==2),
        "original_studies_with_mule_higher_observed_night_relative_q_through_three_clocks":q_order_studies,
        "both_studies_preserve_native_diel_q_ordering_if_supported":(
            len(complete)==2 and len(q_order_studies)==2),
        "all_study_source_original_native_photoperiod_sign_results":records,
        "source_sample_is_random_world_ecosystem":False,
        "previous_12_ERA5_source_raw_HOLDS_UNCHANGED":True,
        "previous_48_source_shortwave_direction_external_nonreplication_UNCHANGED":NEGATIVE,
        "interpretation_ceiling":[
            "Utah pair of mule deer/bighorn species and the shared 43N daylength coordinate were chosen AFTER observing the Wyoming result, not blinded preregistration.",
            "Same public taxa across two distinct Movebank source DOI represent at most two field-study contexts, not four independent ecological/phylogenetic species-level replications.",
            "Within each source the calendar year dates are shared between species, but tags/duty cycles/habitats and actual irradiance may differ; across Utah versus Wyoming exact calendar dates/years need NOT match.",
            "An insufficient Utah source shared-date sample is a missing/attrition endpoint, not a negative species behavioral biological response.",
            "Original q is night minus day observed GPS straight-line displacement proxy, not actual duration of activity or within-individual physiological photoperiod sensitivity.",
            "Source-wide mule-deer relative night-minus-day observed GPS displacement q ordering is descriptive and a separate estimand from the seasonal q~daylength slope; any cross-study recurrence is post-Wyoming-result discovery, not prospectively preregistered confirmation.",
            "A zero-spanning nominal bootstrap CI is unresolved and not proof of absence or sign equivalence. No multiplicity-corrected or external population confidence claim.",
            "Historical 12 ERA5 raw response source SHA HOLDs and separate 48-source negative signed radiation external validation unchanged; no GPS telemetry/IDs/date-level time/location redistributed."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,default=ROOT/"research/RSDS_UTAH_WYOMING_TWO_STUDY_SHARED_DATE_PHOTOPERIOD_V1_RESULT_20261029.json")
    args=cli.parse_args()
    rec=json.loads((ROOT/PINS).read_text(encoding="utf-8"))
    if rec.get("stage")!="UTAH_FIRST_RESULT_SEALED_BEFORE_SECOND_STUDY_COMPARISON":
        raise ValueError("Utah first result Git blob SHA has not been sealed")
    utah=load_blob(UTAH,rec["utah_result_git_blob_sha"])
    wyoming=load_blob(WYOMING,WYOMING_SHA)
    result=compare(utah,wyoming)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RSDS_UTAH_WYOMING_NATIVE_PHOTOPERIOD",json.dumps({
        "distinct_original_source_DOI":result["fixed_original_study_DOIs"],
        "fully_eligible_sources":result["source_original_studies_full_calendar_comparison_identifiable"],
        "same_native_opposite_photoperiod_sources":
            result["original_studies_with_same_opposite_native_daylength_sign_through_three_clocks"],
        "old_external_nonreplication":result["previous_48_source_shortwave_direction_external_nonreplication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
