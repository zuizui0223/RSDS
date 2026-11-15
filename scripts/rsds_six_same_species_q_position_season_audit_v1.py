"""Derivative-only comparison of same-species q niche position and photoperiod slope.

Inputs are SHA-frozen ORIGINAL source-native study aggregates. Never re-read
original Movebank GPS or provider files. Sign interpretation is post-outcome;
within-source animal-bootstrap CIs do not infer source-to-source differences.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned

SPEC=HERE/"research/RSDS_SAME_SPECIES_Q_POSITION_SEASONAL_SLOPE_AUDIT_CONTRACT_20261115.json"
VERSION="RSDS_SAME_SPECIES_Q_POSITION_VS_DAYLENGTH_SLOPE_POSTOUTCOME_V1"


def check_effect(e):
    if not isinstance(e,dict):
        raise ValueError("Original source-native result missing")
    if e.get("available") is False:
        if e.get("q") is not None or e.get("nominal_source_CI_sign")!="UNAVAILABLE":
            raise ValueError("Unavailable original source q/daylength improperly recoded")
        return {"available":False,"mean":None,"ci":None,"n":None,"point":"UNAVAILABLE","supported":"UNAVAILABLE"}
    if e.get("available") is not True:
        raise ValueError("Unknown source effect availability status")
    q=e["q"]
    ci=e.get("source_tag_bootstrap95")
    n=e.get("original_source_tagged_animals")
    if (not isinstance(q,(int,float)) or not isinstance(n,int) or n<5
        or not isinstance(ci,list) or len(ci)!=2 or not all(isinstance(x,(float,int)) for x in ci)
        or not ci[0]<=ci[1]):
        raise ValueError("Invalid original source-tag q/photoperiod CI")
    point="POSITIVE" if q>0 else "NEGATIVE" if q<0 else "ZERO"
    supported="POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
    if e["point_sign"]!=point or e["nominal_source_CI_sign"]!=supported:
        raise ValueError("Original q/photoperiod bootstrap sign changed")
    return {"available":True,"mean":float(q),"ci":list(ci),"n":n,"point":point,"supported":supported}


def read_inputs():
    c=json.loads(SPEC.read_text(encoding="utf-8"))
    if (c.get("version")!=VERSION
        or c.get("prior_original_63_DOIs")!=63
        or c.get("prior_source_native_records")!=88
        or c.get("new_source_IDs")!=6
        or c.get("new_file_MD5_admissible")!=4
        or c.get("new_original_q_estimated")!=3
        or c.get("new_original_daylength_slope_estimated")!=3
        or c.get("no_new_original_animal_GPS_reads") is not True
        or c.get("old_48_source_sunlight_heading_direction_external_not_replicated") is not True
        or c.get("old_12_ERA5_raw_JSON_SHA_hold") is not True):
        raise ValueError("Frozen derivative-only science contract changed")
    old=pinned(c["prior_63_source_native_q_path"],c["prior_63_source_native_q_git_blob"])
    new=pinned(c["six_original_first_q_result_path"],c["six_original_first_q_result_git_blob"])
    prior_contract=pinned(c["six_original_q_pre_event_contract_path"],c["six_original_q_pre_event_contract_git_blob"])
    exposure=pinned(c["six_original_prior_source_exposure_path"],c["six_original_prior_source_exposure_git_blob"])
    source=pinned(c["six_original_source_only_manifest_path"],c["six_original_source_only_manifest_git_blob"])
    metadata=pinned(c["six_original_provider_file_meta_path"],c["six_original_provider_file_meta_git_blob"])
    if (old["distinct_original_archived_source_DOIs_and_study_IDs"]!=63
        or old["all_registered_original_source_native_rows"]!=88
        or len(old["same_species_source_study_q_comparison_candidates"])!=6
        or len(new["all_six_original_source_taxon_results_including_HOLD"])!=6
        or new["original_previous_archive_source_DOIs"]!=63
        or new["new_original_source_q_supported"]!=3
        or new["new_original_source_daylength_supported"]!=3
        or len(exposure["original_followup_source_exposure_rows"])!=6
        or exposure["global_prior_GPS_exposure_audit_completed"] is not False
        or len(metadata["original_source_metadata_rows"])!=6
        or metadata["source_original_metadata_status_counts"]!={"ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":4,"PREVIOUS_SOURCE_METADATA_HOLD_PRESERVED":1,"SOURCE_METADATA_HOLD":1}
        or prior_contract["source_candidates_original"]!=6
        or prior_contract["source_original_file_MD5_admissible"]!=4
        or prior_contract["source_original_file_HOLD"]!=1
        or prior_contract["prior_original_source_metadata_HOLD_preserved"]!=1
        or prior_contract["original_48_signed_irradiance_heading_external_nonreplicated"] is not True
        or new["previous_original_48_solar_direction_signed_heldout_nonreplication"] is not True
        or new["previous_original_12_ERA5_raw_source_HOLD_unchanged"] is not True):
        raise ValueError("Original past/native and first follow-up sources changed")
    return c,old,new,prior_contract,exposure,source,metadata


def audit():
    c,old,new,first_contract,exposure,selection,files=read_inputs()
    target={x["scientific_name"]:x for x in old["same_species_source_study_q_comparison_candidates"]}
    if sorted(target)!=sorted(c["target_species_names"]):
        raise ValueError("Prior six target source-native species changed")
    registered=selection["original_source_metadata_eligible_six"]
    if registered!=6:
        raise ValueError("Original six provider source identities changed")
    s_metadata={x["original_source_DOI"]:x for x in files["original_source_metadata_rows"]}
    s_exposure={x["original_DOI"]:x for x in exposure["original_followup_source_exposure_rows"]}
    new_sources=new["all_six_original_source_taxon_results_including_HOLD"]
    if len(set(z["original_source_DOI"] for z in new_sources))!=6:
        raise ValueError("Duplicate new study DOI")
    all_expected=first_contract["source_q_primary_prediction_sign_by_DOI"]
    if len(all_expected)!=4:
        raise ValueError("First q sign expectations changed")
    out=[]
    for r in new_sources:
        taxon=r["scientific_name"]
        doi=r["original_source_DOI"]
        if taxon not in target or doi not in s_metadata or doi not in s_exposure:
            raise ValueError("Unregistered target study or source exposure")
        if r["source_index"] not in range(6):
            raise ValueError("Out-of-selection source index")
        sm=s_metadata[doi]
        se=s_exposure[doi]
        if (sm["original_repository_study_id"]!=r["original_source_repository_study_ID"]
            or sm["original_source_DOI"]!=doi
            or se["original_DOI"]!=doi
            or r["known_previous_original_GPS_exposure_status"]!=se["prior_data_exposure_status"]
            or r["source_status"]=="SOURCE_Q_ESTIMABLE" and sm["source_metadata_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"):
            raise ValueError("Source original ID, provider file or prior exposure changed")
        prior=target[taxon]
        priors=[{
            "DOI":z["original_DOI"],
            "native_q":check_effect(z["original_native_q"]),
            "native_daylength":check_effect(z["original_native_daylength"]),
            "original_archive_title":z["original_catalog_title"]
        } for z in prior["source_native_rows"] if z["original_native_q"]["available"]]
        if len(priors)!=prior["source_native_q_estimable_DOI"]:
            raise ValueError("Old q-supported source count changed")
        nq=check_effect(r["original_new_source_q"])
        nd=check_effect(r["original_new_source_q_daylength_slope"])
        required=all_expected.get(doi)
        if required is not None:
            if required not in {"POSITIVE","NEGATIVE"} or r["prediction_frozen_from_past_results"]!=required:
                raise ValueError("Frozen old-q sign prediction changed")
        elif r["prediction_frozen_from_past_results"] is not None:
            raise ValueError("Previously no original q sign existed but a prediction was added")
        primary_verdict=(
            "NOMINAL_SUPPORTED_SAME_SIGN" if nq["available"] and
            nq["supported"]==required and required in ("POSITIVE","NEGATIVE")
            else "NOMINAL_OPPOSITE" if nq["available"] and
            nq["supported"] not in ("INCLUDES_ZERO","UNAVAILABLE") and
            nq["supported"]!=required and required is not None
            else "NOMINAL_CI_INCLUDES_ZERO" if nq["available"] and nq["supported"]=="INCLUDES_ZERO"
            else "UNAVAILABLE_OR_ORIGINAL_SOURCE_HOLD"
        )
        if r["nominal_source_CI_repeatability_verdict"]!=primary_verdict:
            raise ValueError("Previously fixed first q verdict reclassified")
        old_supported=[z for z in priors if z["native_daylength"]["supported"] in ("POSITIVE","NEGATIVE")]
        opposite_slope=[z for z in old_supported
                        if nd["supported"] in ("POSITIVE","NEGATIVE") and
                        z["native_daylength"]["supported"]!=nd["supported"]]
        same_slope=[z for z in old_supported
                    if nd["supported"] in ("POSITIVE","NEGATIVE") and
                    z["native_daylength"]["supported"]==nd["supported"]]
        # "Opposite" here only means that two independent ORIGINAL within-source
        # tag-bootstrap intervals exclude zero on opposite sides. It does NOT
        # certify a source-to-source difference, or tag/population independence.
        out.append({
            "original_source_index":r["source_index"],
            "species":taxon,
            "original_new_DOI":doi,
            "original_new_source_status":r["source_status"],
            "new_q":nq,
            "new_q_daylength":nd,
            "frozen_expected_q_point_sign":required,
            "first_new_source_nominal_q_verdict":primary_verdict,
            "prior_q_supported_archives":priors,
            "prior_q_supported_archive_count":len(priors),
            "prior_daylength_nominal_supported_count":len(old_supported),
            "prior_season_opposite_new_sign_supported_source_DOIs":[z["DOI"] for z in opposite_slope],
            "prior_season_same_new_sign_supported_source_DOIs":[z["DOI"] for z in same_slope],
            "new_daylength_opposite_prior_nominal_sign":bool(opposite_slope),
            "prior_original_solar_GPS_exposure_status":r["known_previous_original_GPS_exposure_status"],
            "provider_MD5_original_source_verified":r["original_source_MD5_verified"],
            "original_animal_tags_or_population_independent_verified":False,
        })
    out.sort(key=lambda z:z["original_source_index"])
    if [z["original_source_index"] for z in out]!=list(range(6)):
        raise ValueError("Missing one of six original source IDs")
    status=Counter(z["first_new_source_nominal_q_verdict"] for z in out)
    q_supported=[z for z in out if z["new_q"]["available"]]
    seasonal=[z for z in q_supported if z["new_daylength"]["available"]]
    opposed=[z for z in seasonal if z["new_daylength_opposite_prior_nominal_sign"]]
    no_prior_fox=[z for z in q_supported if "KNOWN_PREVIOUS_ORIGINAL_GPS" not in z["prior_original_solar_GPS_exposure_status"]]
    if (len(q_supported)!=3 or len(seasonal)!=3 or len(opposed)!=2 or
        status!={"NOMINAL_SUPPORTED_SAME_SIGN":3,"UNAVAILABLE_OR_ORIGINAL_SOURCE_HOLD":3} or
        {z["species"] for z in opposed}!={"Anser albifrons","Odocoileus hemionus"}):
        raise ValueError("Frozen first native q/photoperiod sign support changed")
    return {
        "version":VERSION,
        "stage":"POST_FIRST_Q_RESULTS_SOURCE_NATIVE_POSITION_VS_DAYLENGTH_SLOPE_CONTRAST",
        "original_previous_source_studies":63,
        "previous_q_supported_native_rows":38,
        "original_new_source_study_slots":6,
        "original_source_file_admissible":4,
        "original_source_file_HOLD":2,
        "original_q_supported_new":3,
        "original_daylength_supported_new":3,
        "new_source_bootstrap_primary_q_sign_verdicts":dict(status),
        "q_same_sign_supported_new_species":[z["species"] for z in q_supported
                                            if z["first_new_source_nominal_q_verdict"]=="NOMINAL_SUPPORTED_SAME_SIGN"],
        "new_q_and_daylength_supported_source_studies":len(seasonal),
        "species_with_opposite_nominal_daylength_slope_support_across_archived_sources":[z["species"] for z in opposed],
        "opposite_nominal_daylength_slope_support_species_count":len(opposed),
        "new_same_q_sign_with_prior_GPS_exposure_unresolved_count":len(no_prior_fox),
        "new_same_q_sign_with_known_prior_GPS_exposure_count":len(q_supported)-len(no_prior_fox),
        "all_original_new_source_and_prior_q_source_daylength_rows":out,
        "old_48_source_signed_solar_direction_heldout_NOT_REPLICATED":True,
        "old_12_ERA5_provider_raw_SHA_HOLD":True,
        "no_new_GPS_or_original_provider_metadata_access":True,
        "post_outcome_comparison_not_a_registered_daylength_reversal_test":True,
        "source_to_source_tag_and_population_independence_unverified":True,
        "source_to_source_slope_difference_interval_not_estimated":True,
        "no_common_seasonal_calendar_control":True,
        "no_global_random_wildlife_species_estimate":True,
        "interpretation_limit":c["source_independence_limitation"],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    d=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_Q_POSITION_SEASON_CONTRAST",json.dumps({
        "new_q_support":d["original_q_supported_new"],
        "sign_support":d["new_source_bootstrap_primary_q_sign_verdicts"],
        "opposite_nominal_daylength_species":d["species_with_opposite_nominal_daylength_slope_support_across_archived_sources"],
        "prior_GPS_exposed_q_new":d["new_same_q_sign_with_known_prior_GPS_exposure_count"],
        "original_animal_GPS_rows_reopened":0,
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
