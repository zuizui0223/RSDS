"""Full original-source versus matched-cadence temporal ecological sign diagnostic.

No new individual GPS data. Compares two Git-blob-pinned ORIGINAL source
aggregate data artifacts:
(1) first 8-study all-qualified-animal/date q and q~photoperiod,
(2) 8-study q and photoperiod under exactly same tagged animal/year/dates in
    three different GPS fix-clock acquisition schedules.

The first frame and the matched frame have DIFFERENT original-tag support.
A stable nominal sign in the 3 matched schedules but not the full frame does
NOT demonstrate species-wide ecological generality; a lost sign does NOT
establish biological nonreplication because of major sample attrition.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OLD="research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json"
MATCHED="research/RSDS_TRIPLE_DATE_ECOLOGICAL_Q_DAYLENGTH_INVARIANCE_V1_FIRST_RESULT_20261028.json"
SHA={
    OLD:"2c4670b54fa5b52737993ebf6babc3a7a58b166a",
    MATCHED:"157af4cacbe4359b07f4cedbed2793f8c49972b3",
}
OLD_SIGN={
    "POSITIVE_NOMINAL_SOURCE_TAG_INTERVAL":"POSITIVE",
    "NEGATIVE_NOMINAL_SOURCE_TAG_INTERVAL":"NEGATIVE",
    "INCLUDES_ZERO":"UNRESOLVED",
    "UNAVAILABLE":"UNAVAILABLE"
}
P="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
OUTPUT=ROOT/"research/RSDS_TRIPLE_DATE_VS_ORIGINAL_ECOLOGICAL_TIME_SIGN_AUDIT_20261028.json"


def source(root,p):
    raw=(root/p).read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
    if sha!=SHA[p]:
        raise ValueError("Original source first scientific result changed: "+p)
    return json.loads(raw.decode("utf-8"))


def compare_sign(original, matched_stability):
    label=OLD_SIGN[original]
    if matched_stability is None or matched_stability=="UNAVAILABLE":
        return "UNAVAILABLE_TRIPLE_SAME_TAG_DATE"
    if matched_stability=="THREE_SCHEDULES_NOMINAL_POSITIVE":
        new="POSITIVE"
    elif matched_stability=="THREE_SCHEDULES_NOMINAL_NEGATIVE":
        new="NEGATIVE"
    elif matched_stability=="UNRESOLVED_AT_LEAST_ONE_ZERO_SPANNING_INTERVAL":
        return "UNRESOLVED_IN_ONE_OR_MORE_MATCHED_GPS_SCHEDULES"
    elif matched_stability=="OPPOSITE_NOMINAL_SOURCE_TAG_SIGN":
        return "OPPOSITE_SIGN_AMONG_MATCHED_GPS_SCHEDULES"
    else:
        raise ValueError("Unknown original tagged animal matched schedule source stability")
    if label==new:
        return "SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_"+new
    if label=="UNRESOLVED":
        return "MATCHED_ONLY_SIGN_SUPPORTED_NOT_FULL_NATIVE_SOURCE_"+new
    if label=="UNAVAILABLE":
        raise ValueError("Full source q missing despite identifiable matched original tags")
    return "OPPOSITE_FULL_VS_MATCHED_NATIVE_SOURCE_SIGN"


def run(root=ROOT):
    first=source(root,OLD)
    matched=source(root,MATCHED)
    if (first["original_study_denominator"]!=8
            or first["source_taxon_rows"]!=12
            or matched["original_source_studies"]!=8
            or matched["original_study_taxon_rows"]!=12
            or not matched["historical_original_12_ERA5_source_version_HOLDS_UNCHANGED"]
            or matched["historical_48_source_signed_shortwave_external_NONREPLICATED_UNCHANGED"]!=P
            or first["old_48_radiation_study_external_negative_direction_test_UNCHANGED"]!=P):
        raise ValueError("Earlier source selected biological study population was altered")
    old={(x["source_study_doi"],x["taxon"]):x for x in first["all_sources_including_egret_original_metadata_HOLD"]}
    if len(old)!=12:
        raise ValueError("The original source taxonomy has duplicated identity")
    rows=[]
    totals={"q":Counter(),"daylength":Counter()}
    supported_studies={"q":set(),"daylength":set()}
    retained_three={}
    for study in matched["all_fixed_original_studies_and_native_taxa"]:
        doi=study["original_study_DOI"]
        for r in study["native_source_results"]:
            key=(doi,r["taxon"])
            origin=old.pop(key,None)
            if origin is None:
                raise ValueError("Missing first original public study/taxon source")
            prior_q=origin["q"]
            prior_l=origin["daylength_slope"]
            strict_q=r.get("original_same_date_q_source_estimates")
            strict_l=r.get("original_same_date_daylength_slope_estimates")
            if strict_q is not None and set(strict_q)!={"q_original","q_even","q_odd"}:
                raise ValueError("Missing one of original GPS q schedules")
            if strict_l is not None and set(strict_l)!={"daylength_original","daylength_even","daylength_odd"}:
                raise ValueError("Missing original q daylight schedule")
            n_prior=prior_q["source_animals"]
            n_matched=(strict_q or {}).get("q_original",{}).get("n_source_tags",0)
            n_grad=(strict_l or {}).get("daylength_original",{}).get("n_source_tags",0)
            if n_matched>n_prior or n_grad>prior_l["source_animals"]:
                raise ValueError("Identical supported original tags more numerous than full first source")
            qa=compare_sign(prior_q["status"],r.get("source_q_three_schedule_nominal_sign_stability"))
            la=compare_sign(prior_l["status"],r.get("source_daylength_three_schedule_nominal_sign_stability"))
            totals["q"][qa]+=1
            totals["daylength"][la]+=1
            if qa.startswith("SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_"):
                supported_studies["q"].add(doi)
            if la.startswith("SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_"):
                supported_studies["daylength"].add(doi)
            rows.append({
                "original_study_DOI":doi,"native_taxon":r["taxon"],
                "original_source_status":study["status"],"original_triple_support_status":r["status"],
                "full_original_source_q_n_tagged_animals":n_prior,
                "original_identical_triple_matched_q_n_tagged_animals":n_matched,
                "original_same_tag_day_retention_fraction_of_full_q_tags":
                    n_matched/n_prior if n_prior else None,
                "full_original_source_q":prior_q,
                "original_same_date_q_level":strict_q["q_original"] if strict_q else None,
                "same_tag_three_schedule_q_sign":r.get("source_q_three_schedule_nominal_sign_stability"),
                "full_vs_matched_source_q_inference_class":qa,
                "full_original_source_q_daylength_slope":prior_l,
                "original_same_date_q_daylength_slope":strict_l["daylength_original"] if strict_l else None,
                "same_tag_three_schedule_daylength_sign":r.get("source_daylength_three_schedule_nominal_sign_stability"),
                "full_vs_matched_source_daylength_inference_class":la,
            })
    if old or len(rows)!=12:
        raise ValueError("One original native study/taxon missing from full versus matched source report")
    if (len(supported_studies["q"])!=5
            or len(supported_studies["daylength"])!=2
            or totals["q"]["SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_NEGATIVE"]!=2
            or totals["q"]["SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_POSITIVE"]!=3
            or totals["daylength"]["SAME_SIGN_SUPPORTED_IN_FULL_AND_ALL_THREE_SCHEDULES_POSITIVE"]!=2):
        raise ValueError("Sign stability of first original ecological source study changed")
    return {
        "version":"RSDS_FULL_SOURCE_VS_MATCHED_GPS_TEMPORAL_ECOLOGY_SIGN_AUDIT_V1",
        "date":"2026-10-28",
        "research_stage":"POST_BOTH_FIRST_RESULTS_SOURCE_SUPPORT_REWEIGHTING_VALIDATION",
        "original_field_study_doi_count":8,
        "source_taxon_rows":12,
        "first_full_source_git_blob_receipts":SHA,
        "full_vs_triple_q_support_sign_count":dict(totals["q"]),
        "full_vs_triple_daylength_slope_sign_count":dict(totals["daylength"]),
        "n_distinct_original_field_studies_with_q_sign_supported_in_full_and_all_three_schedules":len(supported_studies["q"]),
        "n_distinct_original_field_studies_with_daylength_sign_supported_in_full_and_all_three_schedules":len(supported_studies["daylength"]),
        "full_source_and_strict_same_original_tag_date_results":rows,
        "old_12_ERA5_serialized_raw_SHA_HOLDS_UNCHANGED":True,
        "old_48_source_external_signed_radiation_nonreplication_UNCHANGED":P,
        "science_ceiling":[
            "Full source q and strict triple-identical original-tag-year-date q are different supported source populations/periods; a change of CI or sign can arise from sampling availability, not a biological response to clock time.",
            "The same source nominal sign in full and three schedule subsamples is descriptive conditional evidence of q displacement proxy sign robustness, not ecological identity or globally transferable species/phylogenetic parameter.",
            "Only source original tagged animal bootstrap unadjusted 95% intervals; prior different source studies and local habitats are not globally randomized and taxa within same source study not independent ecological replications.",
            "In the triple schedule, original GPS collection cadence and observed route chord remain noncomparable 0.5-2hr vs 1.5-2.5hr; q is observed displacement not true time in activity.",
            "The original turkey and fishing cat have full source q estimates but NO strict original 3-schedule supported tags; cannot call their original effects refuted or absent.",
            "The original Great Egret metadata source HOLD and old ERA5 provider raw-SHA source-version HOLDS stay untouched. Older negative-direction shortwave external holdout source test remains NOT REPLICATED.",
            "No original GPS animal IDs, exact original calendar dates/UTC times, coordinates, raw telemetry or per-tag slopes made public."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,default=OUTPUT)
    args=cli.parse_args()
    result=run()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_FULL_VS_MATCHED_ORIGINAL_TIME_ECOLOGY",json.dumps({
        "studies":result["original_field_study_doi_count"],
        "source_taxa":result["source_taxon_rows"],
        "q_sign_results":result["full_vs_triple_q_support_sign_count"],
        "daylength_sign_results":result["full_vs_triple_daylength_slope_sign_count"],
        "q_source_studies":result["n_distinct_original_field_studies_with_q_sign_supported_in_full_and_all_three_schedules"],
        "daylength_source_studies":result["n_distinct_original_field_studies_with_daylength_sign_supported_in_full_and_all_three_schedules"],
        "historical_12_ERA5_source_HOLDS":result["old_12_ERA5_serialized_raw_SHA_HOLDS_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
