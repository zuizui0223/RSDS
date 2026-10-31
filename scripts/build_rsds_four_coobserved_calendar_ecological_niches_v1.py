"""Original calendar-matched native diel displacement vs photoperiod across four studies.

Four distinct public Movebank field-study DOIs provide exact co-observed
source calendar day frames: Illinois (3 taxa, 169 days), Wyoming (2,193),
Utah (same native 2,199), and geese (3,120). Each input is an immutable
original-GPS result already original-provider MD5 and source q/daylength
verified. No original GPS events, tags, dates or weather retrieved here.

The point of the comparative source study is to separate *native relative
q position on a common calendar* from *source reference daylight-q slopes*.
Source is DOI not taxon; between-study daylength reference latitudes and
sampling opportunities differ, making a global slope pooling invalid.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"research/RSDS_FOUR_COOOBSERVED_CALENDAR_NICHE_V1_CONTRACT.json"
OUTPUT=ROOT/"research/RSDS_FOUR_COOOBSERVED_CALENDAR_NICHE_V1_FIRST_RESULT_20261031.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def pinned(path,expected):
    raw=(ROOT/path).read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
    if sha!=expected:
        raise ValueError("Previously frozen source original calendar result Git blob changed: "+path)
    return json.loads(raw.decode("utf-8"))


def inputs():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if spec["version"]!="RSDS_ORIGINAL_FOUR_COOOBSERVED_CALENDAR_NICHE_V1":
        raise ValueError("Not fixed source coobserved original study calendar comparison")
    arr=[pinned(x["path"],x["git_blob_sha"]) for x in spec["source_result_files"]]
    if len(arr)!=5 or len(spec["source_DOI_frame"])!=4:
        raise ValueError("Coobserved original study source frame changed")
    ill,twogoose,goose,month,macro=arr
    if not (ill["MD5_verified"] and
           ill["first_all_native_q_daylength_and_solar_source_replayed"] and
           ill["old_12_ERA5_raw_SHA_HOLDS_UNCHANGED"] and
           ill["prior_48_negative_shortwave_direction_external_UNCHANGED"]==NEGATIVE and
           twogoose["source_original_studies_full_calendar_comparison_identifiable"]==2 and
           twogoose["previous_12_ERA5_source_raw_HOLDS_UNCHANGED"] and
           twogoose["previous_48_source_shortwave_direction_external_nonreplication_UNCHANGED"]==NEGATIVE and
           goose["source_original_MD5_verified"]==4 and
           goose["new_geese_original_same_source_calendar_days"]==120 and
           goose["previous_48_source_signed_shortwave_direction_EXTERNAL_NONREPLICATED_UNCHANGED"]==NEGATIVE and
           month["original_studies"]==2 and
           month["old_48_signed_radiation_external_nonreplication_UNCHANGED"]==NEGATIVE and
           macro["all_original_studies"]==12 and
           macro["external_48_source_negative_shortwave_direction_nonreplication_UNCHANGED"]==NEGATIVE and
           macro["original_12_ERA5_raw_JSON_source_HOLDS_UNCHANGED"]):
        raise ValueError("One original source q/calendar/MD5 or historical source HOLD changed")
    return spec,ill,twogoose,goose,month,macro


def signed(point,ci,reported):
    if point is None or not isinstance(ci,(list,tuple)) or len(ci)!=2:
        raise ValueError("Estimable source original q or daylight effect lacks point/CI")
    lo,hi=ci
    if not all(type(x) in (int,float) and math.isfinite(x) for x in (point,lo,hi)) or not lo<=point<=hi:
        raise ValueError("Source original native effect not within original tagged-animal uncertainty")
    s="POSITIVE" if lo>0 else "NEGATIVE" if hi<0 else "INCLUDES_ZERO"
    known=reported.removeprefix("NOMINAL_")
    if known!=s:
        raise ValueError("Native/source original signed study result did not replay")
    return {"mean":point,"original_source_animal_bootstrap95":list(ci),
            "nominal_sign":s}


def one_native(taxon,n,q,day):
    if n<5:
        raise ValueError("Original source calendar native effect without five tagged animals")
    return {"taxon":taxon,"original_tagged_animals_q":n,"source_q":q,
            "source_reference_daylength_slope_per_hour":day}


def normalize_pair(a,b,qval,qci,oldsign,additional=None):
    out={
        "first":a,"second":b,
        "source_shared_calendar_relative_q_difference":signed(qval,qci,oldsign),
        "not_independent_field_study":True
    }
    if additional:
        out.update(additional)
    return out


def make_illinois(i,expect):
    observed=i["source_all_three_native_support"]
    taxa=i["native_taxa"]
    if (i["study_DOI"]!=expect["doi"] or set(taxa)!=set(expect["registered_taxa"]) or
        observed["shared_source_calendar_dates_all_three_native_taxa"]!=expect["shared_source_calendar_days"] or
        observed["shared_source_year_month_blocks"]!=expect["year_month_blocks"]):
        raise ValueError("Previous Illinois same-calendar native source frame changed")
    native=[]
    for taxon in taxa:
        obj=i["source_all_native_q_levels_on_same_original_calendar"][taxon]
        a,b=obj["q_original"],obj["daylength_original"]
        if a["original_tagged_animals"]!=obj["q_source_tagged_animals"] or b["original_tagged_animals"]!=obj["daylength_source_tagged_animals"]:
            raise ValueError("Illinois native same-calendar original tag denominators changed")
        native.append(one_native(taxon,obj["q_source_tagged_animals"],
                                 signed(a["source_animal_mean"],a["source_tag_bootstrap95"],a["nominal_sign"]),
                                 signed(b["source_animal_mean"],b["source_tag_bootstrap95"],b["nominal_sign"])))
    pairs=[]
    for key,z in i["three_source_native_pair_q_and_daylength"].items():
        if {z["native_first"],z["native_second"]}.issubset(set(taxa)) is False:
            raise ValueError("Unregistered Illinois original native pair")
        p=z["q"]
        pairs.append(normalize_pair(z["native_first"],z["native_second"],
                                    p["first_minus_second_mean"],p["original_tag_bootstrap95"],p["status"]))
    if len(pairs)!=3:
        raise ValueError("All three Illinois original native study pairs not retained")
    return native,pairs


def make_ungulate(source,expect):
    taxa=expect["registered_taxa"]
    if (source["original_study_DOI"]!=expect["doi"] or
        set(source["both_native_taxa_original_tags_on_shared_calendar"])!=set(taxa) or
        source["source_study_shared_calendar_dates"]!=expect["shared_source_calendar_days"] or
        source["source_study_year_month_blocks"]!=expect["year_month_blocks"]):
        raise ValueError("One native ungulate field study/source calendar changed")
    originals=source["three_observation_schedules"]["original"]
    native=[]
    for taxon in taxa:
        q=originals["native_q_levels"][taxon]
        d=originals["native"][taxon]
        if q["n_source_original_tags"]!=source["both_native_taxa_original_tags_on_shared_calendar"][taxon] or d["n_source_original_tags"]!=q["n_source_original_tags"]:
            raise ValueError("Ungulate native source tag support changed")
        native.append(one_native(taxon,q["n_source_original_tags"],
                                 signed(q["mean"],q["ci"],q["sign"]),
                                 signed(d["mean"],d["ci"],d["sign"])))
    p=originals["q_contrast"]
    pair=[normalize_pair(taxa[0],taxa[1],p["mule_minus_bighorn"],p["ci"],p["nominal_sign"],
         {"three_fixed_original_UTC_clock_q_contrast_signs":{
             clk:source["three_observation_schedules"][clk]["q_contrast"]["nominal_sign"].removeprefix("NOMINAL_")
             for clk in ("original","even","odd")},
          "reference_daylength_native_contrast":{
             "mean":originals["contrast"]["mule_minus_bighorn"],
             "nominal_sign":originals["contrast"]["nominal_sign"].removeprefix("NOMINAL_"),
             "original_source_tag_bootstrap95":originals["contrast"]["ci"]
          }})]
    for clk in ("original","even","odd"):
        if pair[0]["three_fixed_original_UTC_clock_q_contrast_signs"][clk]!="POSITIVE":
            raise ValueError("Original mule vs bighorn q rank across two source studies changed")
    return native,pair


def make_geese(g,expect):
    if (g["new_coobserved_geese_original_study_DOI"]!=expect["doi"] or
        g["new_geese_original_same_source_calendar_days"]!=expect["shared_source_calendar_days"] or
        g["new_geese_source_year_month_blocks"]!=expect["year_month_blocks"] or
        g["new_coobserved_geese_source_study_native_taxa"]!=3):
        raise ValueError("New goose original shared study/calendar provenance changed")
    native=[]
    for t in g["goose_native_q_and_study_reference_daylength"]:
        native.append(one_native(t["native_taxon"],t["original_same_calendar_q_tag_count"],
                                 signed(t["observed_NIGHT_minus_DAY_GPS_q"],
                                        t["q_original_tag_bootstrap95"],t["q_nominal_interval_sign"]),
                                 signed(t["q_per_hour_study_reference_daylength_slope"],
                                        t["daylength_original_tag_bootstrap95"],t["daylength_nominal_interval_sign"])))
    pairs=[]
    for p in g["all_three_frozen_original_goose_source_native_pairs"]:
        a,b=p["taxa"]
        pairs.append(normalize_pair(a,b,p["source_NIGHT_DAY_q_difference"],
                                    p["source_original_individual_pair_bootstrap95"],
                                    p["pair_q_nominal_sign"],
                                    {"within_month_study_reference_daylength_slope":{
                                        "mean":p["within_month_study_reference_daylength_slope"],
                                        "month_block_bootstrap95":p["descriptive_source_month_block_bootstrap95"],
                                        "nominal_sign":p["month_nominal_sign"],
                                        "year_month_blocks":p["source_year_month_blocks"]
                                    }}))
    if len(native)!=3 or len(pairs)!=3:
        raise ValueError("New coobserved goose three taxa/three pairs not intact")
    return native,pairs


def monthly(month,goose):
    # Month-adjusted old original source study results are exploratory
    # and NOT additional independent study/calendar units.
    results={}
    for source in month["source_study_results"]:
        region=source["original_study_region"]
        if region not in ("Illinois","Wyoming") or not source["MD5_and_native_source_effects_replayed"]:
            raise ValueError("Invalid earlier within-source month conditional original GPS frame")
        vals=[]
        for p in source["native_pair_schedule_results"]:
            if p["schedule"]!="original":
                continue
            v=p["within_year_month_demeaned_q_daylength"]
            vals.append({"native_taxa":[p["native_first"],p["native_second"]],
                         "conditional_reference_daylength_slope":v["within_month_daylength_slope"],
                         "year_month_block_bootstrap95":v["month_cluster_bootstrap95"],
                         "nominal_sign":v["nominal_interval_sign"],
                         "source_year_month_blocks":v["source_year_month_blocks"],
                         "not_a_new_field_study":True})
        if region in results or len(vals)!=(3 if region=="Illinois" else 1):
            raise ValueError("Unexpected old native pair/month adjusted source denominator")
        results[region]=vals
    results["Geese"]=[
        {"native_taxa":[p["first"],p["second"]],
         "conditional_reference_daylength_slope":p["within_month_study_reference_daylength_slope"]["mean"],
         "year_month_block_bootstrap95":p["within_month_study_reference_daylength_slope"]["month_block_bootstrap95"],
         "nominal_sign":p["within_month_study_reference_daylength_slope"]["nominal_sign"],
         "source_year_month_blocks":p["within_month_study_reference_daylength_slope"]["year_month_blocks"],
         "not_a_new_field_study":True}
        for p in goose
    ]
    results["Utah"]="NOT_TESTED_WITHIN_MONTH_DO_NOT_INFER_ZERO"
    return results


def analyze():
    spec,il,ung,goose,month,macro=inputs()
    byregion={r["region"]:r for r in spec["source_DOI_frame"]}
    s_il,p_il=make_illinois(il,byregion["Illinois"])
    un={}
    for a in ung["all_study_source_original_native_photoperiod_sign_results"]:
        region=a["original_field_study"]
        if region not in ("Wyoming","Utah") or region in un:
            raise ValueError("Previous field-study ungulate source identities absent/duplicated")
        un[region]=make_ungulate(a,byregion[region])
    if set(un)!={"Wyoming","Utah"}:
        raise ValueError("Two original ungulate field-study sources missing")
    s_goose,p_goose=make_geese(goose,byregion["Geese"])
    study=[]
    for region,native,pairs in (
        ("Illinois",s_il,p_il),
        ("Wyoming",*un["Wyoming"]),
        ("Utah",*un["Utah"]),
        ("Geese",s_goose,p_goose)):
        s=byregion[region]
        if len(native)!=len(s["registered_taxa"]) or set(z["taxon"] for z in native)!=set(s["registered_taxa"]):
            raise ValueError("Native calendar source study registered taxonomy changed")
        if any(x["source_shared_calendar_relative_q_difference"]["nominal_sign"] not in
               ("POSITIVE","NEGATIVE","INCLUDES_ZERO") for x in pairs):
            raise ValueError("Unknown source native pair signed interval")
        signs={x["source_reference_daylength_slope_per_hour"]["nominal_sign"] for x in native}
        study.append({
            "original_source_field_study":region,
            "original_study_DOI":s["doi"],
            "native_source_taxon_rows":len(native),
            "source_shared_calendar_days":s["shared_source_calendar_days"],
            "source_year_month_blocks":s["year_month_blocks"],
            "daylength_reference":s["astronomical_reference"],
            "all_source_native_taxon_q_and_season":native,
            "all_within_original_study_native_pair_q_rank_contrasts":pairs,
            "source_native_pairs_with_resolved_observed_relative_q_rank":
                sum(p["source_shared_calendar_relative_q_difference"]["nominal_sign"]!="INCLUDES_ZERO"
                    for p in pairs),
            "source_contains_nominally_positive_AND_negative_native_photoperiod_slopes":
                {"POSITIVE","NEGATIVE"}.issubset(signs),
            "no_proven_spatial_sympatry_or_biological_solar_photoperiod_causality":True
        })
    if (len(study)!=4 or len({s["original_study_DOI"] for s in study})!=4 or
        sum(s["native_source_taxon_rows"] for s in study)!=10 or
        len({x["taxon"] for s in study for x in s["all_source_native_taxon_q_and_season"]})!=8):
        raise ValueError("Original four field-study ten native source cohorts not preserved")
    positive=sum(s["source_contains_nominally_positive_AND_negative_native_photoperiod_slopes"]
                 for s in study)
    rank=sum(s["source_native_pairs_with_resolved_observed_relative_q_rank"]>=1 for s in study)
    if positive!=3 or rank!=4:
        raise ValueError("Original source-calendar q ranking / opposing daylight trend evidence changed")
    normal_month=monthly(month,p_goose)
    coobserved=set(x["original_study_DOI"] for x in study)
    across_12={r["original_study_doi"] for r in macro["all_original_study_native_taxon_rows_including_HOLD"]}
    return {
        "version":"RSDS_FOUR_COOOBSERVED_ORIGINAL_STUDY_CALENDAR_DIEL_PHOTOPERIOD_V1_FIRST",
        "date":"2026-10-31",
        "research_stage":"POST_FIRST_ORIGINAL_CALENDAR_STUDY_RESULTS_COMPARATIVE_ECOLOGY",
        "original_source_studies_with_identical_native_source_calendar_days":4,
        "native_study_taxon_rows":10,
        "distinct_native_taxa_labels":8,
        "study_calendar_days_in_separate_study_frames_NOT_INDEPENDENT_DAILY_REPLICATES":
            sum(s["source_shared_calendar_days"] for s in study),
        "original_studies_with_at_least_one_resolved_native_relative_q_rank":rank,
        "original_studies_with_opposite_supported_native_daylength_slope_signs":positive,
        "original_study_calendar_q_daylength_comparisons":study,
        "same_native_mule_deer_minus_bighorn_repeated_pair_across_two_DOIs":{
            "source_studies":2,
            "two_published_source_DOIs":[byregion["Wyoming"]["doi"],byregion["Utah"]["doi"]],
            "all_three_source_GPS_UTC_schedules_show_positive_q_pair_rank":True,
            "opposite_native_daylength_signs_supported_only_in_Wyoming":True,
            "Utah_native_daylength_signs_both_positive_and_pair_contrast_unresolved":True,
            "old_Wyoming_source_inspected_before_Utah_pair_selection":True
        },
        "native_original_month_conditional_pair_source_comparisons":normal_month,
        "four_calendar_study_DOIs_already_in_twelve_study_full_source_macro":len(coobserved&across_12),
        "Illinois_original_DOI_is_NOT_part_of_twelve_full_source_study_frame":byregion["Illinois"]["doi"] not in across_12,
        "old_12_ERA5_source_version_raw_JSON_HOLDS_UNCHANGED":True,
        "old_48_source_negative_radiation_direction_external_replication_UNCHANGED":NEGATIVE,
        "source_result_git_blob_receipts":{x["path"]:x["git_blob_sha"] for x in spec["source_result_files"]},
        "limits":[
            "Four published study DOIs are the replication unit, not ten native study-taxon rows or 681 separately sampled calendar-day opportunities. Some source field studies or animals may share ecological history.",
            "The same calendar days among source taxa DO NOT prove spatial sympatry, identical habitat, exposure, monitoring effort, light interception, shared resource availability or competition.",
            "Native q is NIGHT−DAY log1p observed straight-line GPS displacement velocity, NOT true biological activity duration or behavior state.",
            "A study-reference astronomical daylength is deterministic with season/date and differs between goose study and 43N ungulate/Illinois references. Comparing coefficients as identical exposures is unjustified.",
            "Opposing nominal source tagged-animal daylight slopes may arise from sample selection, migration route, annual season, GPS duty-cycle, ecology and observation processes; the analysis is exploratory and unadjusted for multiple original native pairs.",
            "A monthly block bootstrap can have only five to nine blocks. Results are source conditional, not worldwide/phylogenetic or causal evidence.",
            "The two-study repeated mule deer–bighorn q rank is an observational replication of a signed within-study position, NOT a blind external validation of daylength mechanisms; Utah pair selected after prior Wyoming result.",
            "Illinois coobserved study is an additional original DOI outside the separately defined 12-study full-source q macro cohort. Do not silently call this a combined 13-study homogeneous sample.",
            "Previous independently nonreplicated 48-source radiation direction and historic ERA5 raw JSON 12 context source HOLDs remain unchanged.",
            "No private animal IDs, GPS coordinates, original observation dates, tracks or individual effects are published."
        ],
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=OUTPUT)
    args=parser.parse_args()
    d=analyze()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_FOUR_COOOBSERVED_CALENDAR_ECOLOGY",json.dumps({
        "original_studies":d["original_source_studies_with_identical_native_source_calendar_days"],
        "native_taxa":d["native_study_taxon_rows"],
        "study_with_q_ranks":d["original_studies_with_at_least_one_resolved_native_relative_q_rank"],
        "study_with_opposing_daylight":d["original_studies_with_opposite_supported_native_daylength_slope_signs"],
        "mule_bighorn_studies":d["same_native_mule_deer_minus_bighorn_repeated_pair_across_two_DOIs"]["source_studies"],
        "old_neg":d["old_48_source_negative_radiation_direction_external_replication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
