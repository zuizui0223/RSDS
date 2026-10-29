"""SHA-pinned original four-study multi-native GPS ecological evidence audit.

Source-first four studies and ten registered taxa are retained. Only ONE
independent original goose field study supports >1 native q taxa on all-taxa
shared 120 calendar dates under the fixed 3 DAY/3 NIGHT GPS step/day rule.
A zero count for exact Galápagos registered canonical taxon labels is a
SOURCE TAXON-SCHEMA HOLD, not proof of no tortoises or no GPS movement.

No unarchived animal GPS/events, source coordinates, tag IDs or exact dates
accessed here. Existing Illinois/Wyoming science and 48-source negative
radiation external result are untouched.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE="research/RSDS_FOUR_NEW_MULTINATIVE_GPS_SHARED_CALENDAR_V1_FIRST_RESULT_20261029.json"
SOURCE_BLOB="4034d5318c914dc789a7fcd6535eedbf3a229d53"
METADATA="research/RSDS_FOUR_NEW_MULTINATIVE_SOURCE_METADATA_FIRST_V1_20261029.json"
METADATA_BLOB="d476c8500f99030c864c382fcf32694e0450a48e"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
OUTPUT=ROOT/"research/RSDS_FOUR_MULTINATIVE_ORIGINAL_GPS_ECOLOGICAL_VALIDATION_V1_20261029.json"


def pinned(root,path,expected):
    raw=(root/path).read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
    if sha!=expected:
        raise ValueError("Original pre-event source metadata or first GPS science blob changed")
    return json.loads(raw.decode("utf-8"))


def infer(root=ROOT):
    root=Path(root)
    result=pinned(root,SOURCE,SOURCE_BLOB)
    metadata=pinned(root,METADATA,METADATA_BLOB)
    if (metadata["original_published_study_count"]!=4
            or metadata["source_only_admissible_original_files"]!=4
            or metadata["animal_GPS_events_read_by_source_only_preflight"]!=0
            or result["fixed_original_published_studies"]!=4
            or result["registered_original_study_native_taxon_rows"]!=10
            or not result["prior_12_ERA5_version_source_HOLDS_UNCHANGED"]
            or result["prior_48_source_shortwave_direction_NEGATIVE_EXTERNAL_REPLICATION_UNCHANGED"]!=NEG):
        raise ValueError("Original four source native taxonomy or negative study outcome modified")
    by={r["published_DOI"]:r for r in result["all_frozen_source_study_original_results"]}
    fixed={
      "10.5441/001/1.2cp86266":3,
      "10.5441/001/1.3fv21n7m":2,
      "10.5441/001/1.mf903197":2,
      "10.5441/001/1.fk899541":3,
    }
    if set(by)!=set(fixed) or any(len(by[k]["source_registered_native_taxa"])!=v for k,v in fixed.items()):
        raise ValueError("Selected four original source DOI/native species cohort frame changed")
    if any(not r["source_provider_original_MD5_replayed"] for r in by.values()):
        raise ValueError("One of four original source provider MD5 files not verified")
    tortoise=by["10.5441/001/1.2cp86266"]
    if (tortoise["source_status"]!="SOURCE_NO_VALID_GPS_FIXES" or
        any(t.get("original_source_metadata_taxonomy_coverage",{}).get("original_taxon_rows")!=0
            for t in tortoise["original_native_source_q_and_daylength"])):
        raise ValueError("Original Galapagos SOURCE taxonomy mismatch/zero exact labels not preserved")
    ducks=by["10.5441/001/1.3fv21n7m"]
    if ducks["original_all_native_common_calendar_support"]["source_calendar_dates_shared_by_all_registered_taxa"]!=0:
        raise ValueError("Original duck shared native calendar noneligibility lost")
    mallard=next(x for x in ducks["original_native_source_q_and_daylength"]
                 if x["taxon"]=="Anas platyrhynchos")
    if not (mallard["original_q_full_source_tag_effect"]["n_individuals"]==21 and
            mallard["original_q_full_source_tag_effect"]["animal_bootstrap_95_percentile"][1]<0):
        raise ValueError("Observed duck baseline native source q was changed")
    vulture=by["10.5441/001/1.mf903197"]
    if (vulture["original_all_native_common_calendar_support"]["source_calendar_dates_shared_by_all_registered_taxa"]!=0
            or any(t["original_q_full_source_tag_effect"] is not None
                   for t in vulture["original_native_source_q_and_daylength"])):
        raise ValueError("Original African vulture source q support changed")
    geese=by["10.5441/001/1.fk899541"]
    shared=geese["original_all_native_common_calendar_support"]
    if (geese["source_status"]!="ORIGINAL_MULTINATIVE_SHARED_CALENDAR_COMPLETED"
            or shared["source_calendar_dates_shared_by_all_registered_taxa"]!=120
            or shared["source_yearmonth_blocks"]!=5):
        raise ValueError("New original 120 co-observed goose species dates changed")
    taxa=geese["original_native_source_q_and_daylength"]
    expected={"Branta leucopsis":62,"Anser albifrons":11,"Anser anser":10}
    if {t["taxon"]:t["same_study_all_native_common_calendar_q"]["qualified_original_tagged_animals"]
        for t in taxa}!=expected:
        raise ValueError("Goose original source tagged animal support changed")
    goose_rows=[]
    for t in taxa:
        q=t["same_study_all_native_common_calendar_q"]
        slope=t["same_study_all_native_common_calendar_daylength_slope"]
        if q["nominal_interval_sign"]!="NEGATIVE":
            raise ValueError("Co-observed goose native q sign no longer negative")
        if (len(q["source_original_tag_bootstrap95"])!=2 or
            q["source_original_tag_bootstrap95"][1]>=0):
            raise ValueError("Goose source q bootstrap uncertainty changed")
        for a in (q,slope):
            if not math.isfinite(a["mean"]):
                raise ValueError("Original goose native source effect invalid")
        goose_rows.append({
            "native_taxon":t["taxon"],
            "original_same_calendar_q_tag_count":q["qualified_original_tagged_animals"],
            "observed_NIGHT_minus_DAY_GPS_q":q["mean"],
            "q_original_tag_bootstrap95":q["source_original_tag_bootstrap95"],
            "q_nominal_interval_sign":q["nominal_interval_sign"],
            "q_per_hour_study_reference_daylength_slope":slope["mean"],
            "daylength_original_tag_bootstrap95":slope["source_original_tag_bootstrap95"],
            "daylength_nominal_interval_sign":slope["nominal_interval_sign"],
        })
    if {x["native_taxon"]:x["daylength_nominal_interval_sign"] for x in goose_rows}!={
        "Branta leucopsis":"NEGATIVE",
        "Anser albifrons":"INCLUDES_ZERO",
        "Anser anser":"POSITIVE"
    }:
        raise ValueError("Non-universal original goose q~astronomical-daylength signs were reclassified")
    pairs=geese["all_source_native_pairs_original_q_and_season"]
    if len(pairs)!=3:
        raise ValueError("One of three registered original goose pairs omitted")
    each=[]
    for z in pairs:
        q=z["source_original_animal_q_pair"]
        x=z["source_common_day_q_difference_daylength"]["within_month_reference_daylength"]
        if q["native_tag_support"]!=[expected[z["native_first"]],expected[z["native_second"]]]:
            raise ValueError("Native goose pair depends on different original tagged species support")
        if x["eligible_source_month_blocks"]!=5 or x["eligible_source_calendar_days"]!=120:
            raise ValueError("Goose same-study within-month calendar source support changed")
        each.append({
            "taxa":[z["native_first"],z["native_second"]],
            "source_NIGHT_DAY_q_difference":q["first_minus_second_q"],
            "source_original_individual_pair_bootstrap95":q["original_independent_tag_bootstrap95"],
            "pair_q_nominal_sign":q["nominal_interval_sign"],
            "within_month_study_reference_daylength_slope":x["within_month_daylength_slope"],
            "source_year_month_blocks":x["eligible_source_month_blocks"],
            "descriptive_source_month_block_bootstrap95":x["month_cluster_bootstrap95"],
            "month_nominal_sign":x["nominal_interval_sign"],
            "source_month_leave_one_out_slope_range":x["leave_one_month_block_out_range"]
        })
    qcounts={sign:sum(x["pair_q_nominal_sign"]==sign for x in each) for sign in ("POSITIVE","NEGATIVE","INCLUDES_ZERO")}
    slopes={sign:sum(x["month_nominal_sign"]==sign for x in each) for sign in ("POSITIVE","NEGATIVE","INCLUDES_ZERO")}
    if qcounts!={"POSITIVE":0,"NEGATIVE":2,"INCLUDES_ZERO":1} or slopes!={"POSITIVE":0,"NEGATIVE":1,"INCLUDES_ZERO":2}:
        raise ValueError("Original goose source native pair ecological sign pattern changed")
    return {
        "version":"RSDS_FOUR_NEW_MULTI_NATIVE_ORIGINAL_GPS_SOURCE_SCIENTIFIC_VALIDATION_V1",
        "date":"2026-10-29",
        "original_source_studies_selected_before_new_q":4,
        "registered_native_study_taxon_rows":10,
        "source_original_MD5_verified":4,
        "original_multinative_studies_with_any_all_taxa_shared_q_calendar":1,
        "original_studies_with_no_all_native_shared_calendar":2,
        "Galapagos_exact_registered_taxon_label_mapping":"UNRESOLVED_ZERO_EXACT_REGISTERED_TAXON_ROWS_NOT_BIOLOGICAL_ZERO",
        "Galapagos_metadata_source_only_primary_file_selected_and_MD5_verified":True,
        "duck_mallard_full_study_q_eligible_tagged_animals":21,
        "duck_mallard_full_study_negative_q":"NEGATIVE_UNADJUSTED_ORIGINAL_TAG_BOOTSTRAP",
        "African_vulture_source_night_day_q":"UNSUPPORTED_NO_APPROPRIATE_SOURCE_TAGS",
        "new_coobserved_geese_original_study_DOI":"10.5441/001/1.fk899541",
        "new_coobserved_geese_source_study_native_taxa":3,
        "new_geese_original_same_source_calendar_days":120,
        "new_geese_source_year_month_blocks":5,
        "goose_native_q_and_study_reference_daylength":goose_rows,
        "all_three_frozen_original_goose_source_native_pairs":each,
        "goose_native_pair_q_nominal_interval_counts":qcounts,
        "goose_within_month_pair_daylength_nominal_interval_counts":slopes,
        "source_file_snapshot_and_original_behavior_first_result_git_blobs":{
            METADATA:METADATA_BLOB,SOURCE:SOURCE_BLOB
        },
        "previous_12_ERA5_raw_JSON_HOLDS_UNCHANGED":True,
        "previous_48_source_signed_shortwave_direction_EXTERNAL_NONREPLICATED_UNCHANGED":NEG,
        "scientific_limitations":[
            "Only one of the four additional original field studies supports all preselected native taxa on exactly the same source calendar under the predeclared q definition. This is meaningful sampling qualification, not absence of ecological differences in the other 3.",
            "Galapagos primary original file MD5 matches, but none of its three registered taxon labels appear under exact canonical taxon matching; diagnose metadata taxonomy/schema before making any animal conclusion. Do not replace the source/file or infer no turtles.",
            "The 120 goose calendar dates do not imply geese occupied the same geography, shared predators or interacted; real individual encountered sunshine, temperature and weather are not measured.",
            "All goose q values are recorded straight-line GPS displacement speed contrasts rather than biological activity hours. Instrument schedules may differ among species within the same study.",
            "A pair q difference has independent tagged animal bootstrap uncertainty, not study-global replication. Source-month conditional q-daylength association can reflect calendar/season rather than causal photoperiod.",
            "Five observed year-month calendar blocks and multiple native pair comparisons are not sufficient for a globally general ecological law or adjusted inference.",
            "Older 8-study macro source ecological results, 12 old ERA5 raw JSON source holds, and 48-study heldout signed radiation negative replication are unchanged."
        ]
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=OUTPUT)
    args=parser.parse_args()
    z=infer()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(z,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_NEW_GOOSE_TEMPORAL_NICHE",json.dumps({
        "fixed_field_studies":z["original_source_studies_selected_before_new_q"],
        "supported_coobserved_studies":z["original_multinative_studies_with_any_all_taxa_shared_q_calendar"],
        "goose_native_N":z["goose_native_q_and_study_reference_daylength"],
        "q_pair_signs":z["goose_native_pair_q_nominal_interval_counts"],
        "month_pair_signs":z["goose_within_month_pair_daylength_nominal_interval_counts"],
        "Galapagos":"HOLD_EXACT_TAXON_MAPPING"
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
