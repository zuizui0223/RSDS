"""Study-ID source-only first q comparison for six fixed post-discovery targets.

The four MD5-admissible original GPS sources are processed in separate
temporary runners before this aggregator. Never reopen GPS here.
Every original source effect may be supported, inconclusive, contradictory,
or unavailable. Never treat unavailable as q=0 or replace source slots.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned
from scripts.rsds_original_source_native_same_species_repeatability_v1 import q_effect
from scripts.run_rsds_six_same_species_followup_original_gps_q_v1 import (
    freeze_original,RESULT_VERSION
)

def compare(e,predicted):
    if not e["available"]:
        return "UNAVAILABLE_OR_ORIGINAL_SOURCE_HOLD"
    if predicted is None:
        return "NOT_PREDECLARED_SIGN_TEST"
    sign=e["point_sign"]
    ci=e["nominal_source_CI_sign"]
    if sign=="ZERO":
        return "POINT_EXACT_ZERO"
    if sign==predicted:
        return "NOMINAL_SUPPORTED_SAME_SIGN" if ci==predicted else "SAME_POINT_SIGN_CI_INCLUDES_ZERO"
    return "NOMINAL_SUPPORTED_OPPOSITE_SIGN" if ci==sign else "OPPOSITE_POINT_SIGN_CI_INCLUDES_ZERO"

def aggregate(original):
    if len(original)!=6:
        raise ValueError("Exactly six original source target slots required")
    contract,sources=freeze_original()
    prior=pinned(contract["original_previous_63_source_native_atlas_path"],
                 contract["original_previous_63_source_native_atlas_git_blob"])
    candidates={x["scientific_name"]:x for x in
                prior["same_species_source_study_q_comparison_candidates"]}
    if (len(candidates)!=6 or
        prior["distinct_original_archived_source_DOIs_and_study_IDs"]!=63 or
        prior["repeated_species_labels_with_two_or_more_estimable_source_native_q"]!=6):
        raise ValueError("Original cross-source repeatability discovery result changed")
    rows=[]
    for i,(s,rec) in enumerate(zip(sources,original)):
        doi=s["original_source_DOI"]
        sp=s["registered_original_native_labels"][0]
        if (rec.get("version")!=RESULT_VERSION or
            rec.get("source_index")!=i or
            rec.get("source_original_DOI")!=doi or
            rec.get("source_original_repository_study_id")!=s["original_repository_study_id"] or
            rec.get("source_original_public_package_uuid")!=s["original_public_package_uuid"] or
            rec.get("source_original_registered_native_taxa")!=s["registered_original_native_labels"] or
            rec.get("source_original_file_metadata_status")!=s["source_metadata_status"] or
            rec.get("prior_original_source_exposure_status")!=s["prior_original_source_exposure_status"] or
            rec.get("prior_original_archive_selected_prediction_sign")!=contract["source_q_primary_prediction_sign_by_DOI"].get(doi) or
            rec.get("no_raw_animal_GPS_IDs_coordinates_exact_dates_or_tag_effects_exported") is not True or
            rec.get("earlier_original_48_solar_direction_signed_external_not_replicated") is not True or
            rec.get("earlier_original_12_ERA5_source_SHA_HOLD") is not True):
            raise ValueError("Original source/taxon/endpoint/HOLD/prediction changed")
        native=rec["native_group_results"]
        if len(native)!=1 or native[0]["taxon"]!=sp:
            raise ValueError("Original registered source native taxon changed")
        status=rec["source_q_science_status"]
        if status not in ("SOURCE_Q_ESTIMABLE","SOURCE_Q_UNAVAILABLE_NOT_ZERO",
                          "SOURCE_GPS_HOLD","SOURCE_FILE_METADATA_HOLD",
                          "HUMAN_LOCATION_PRIVACY_HOLD"):
            raise ValueError("Unknown original source outcome status")
        if s["source_metadata_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":
            if (status!="SOURCE_FILE_METADATA_HOLD"
                or rec["original_file_MD5_bytes_verified"]
                or rec["original_animal_GPS_event_rows_opened"]!=0):
                raise ValueError("Previously frozen original metadata HOLD was reopened")
        else:
            if status=="SOURCE_FILE_METADATA_HOLD":
                raise ValueError("Admissible MD5 source was reclassified to metadata HOLD")
            if status in ("SOURCE_Q_ESTIMABLE","SOURCE_Q_UNAVAILABLE_NOT_ZERO") and not rec["original_file_MD5_bytes_verified"]:
                raise ValueError("Original source q without original provider MD5 confirmation")
        q=q_effect(native[0]["q"])
        day=q_effect(native[0]["daylength"])
        if (q["available"]!=(status=="SOURCE_Q_ESTIMABLE") or
            status=="SOURCE_GPS_HOLD" and (q["available"] or day["available"])):
            raise ValueError("Original source native q mean/support inconsistent with source status")
        if sp not in candidates:
            raise ValueError("Post-discovery source target species was substituted")
        prev=candidates[sp]
        prev_supported=[z for z in prev["source_native_rows"] if
                        z["original_native_q"]["available"]]
        if len(prev_supported)<2 or any(z["original_DOI"]==doi for z in prev_supported):
            raise ValueError("Post-outcome source is not new distinct archived study ID")
        predicted=contract["source_q_primary_prediction_sign_by_DOI"].get(doi)
        if predicted and (not all(z["original_native_q"]["point_sign"]==predicted
                                  for z in prev_supported)):
            raise ValueError("Frozen prior source-native predicted same sign not justified")
        rows.append({
            "source_index":i,
            "original_source_DOI":doi,
            "original_source_repository_study_ID":s["original_repository_study_id"],
            "scientific_name":sp,
            "source_status":status,
            "original_source_MD5_verified":rec["original_file_MD5_bytes_verified"],
            "original_provider_GPS_rows_read_only_privately":rec["original_animal_GPS_event_rows_opened"],
            "original_prior_archive_q_source_count":len(prev_supported),
            "original_prior_archive_q_point_signs":prev["q_point_signs"],
            "prediction_frozen_from_past_results":predicted,
            "original_new_source_q":q,
            "original_new_source_q_daylength_slope":day,
            "nominal_source_CI_repeatability_verdict":compare(q,predicted),
            "known_previous_original_GPS_exposure_status":s["prior_original_source_exposure_status"],
            "source_external_GPS_prior_access_fully_excluded":False,
            "original_animal_tag_population_independence_confirmed":False
        })
    if len({r["original_source_DOI"] for r in rows})!=6:
        raise ValueError("Original source study-DOI repeated")
    counts=Counter(r["nominal_source_CI_repeatability_verdict"] for r in rows)
    return {
        "version":"RSDS_SIX_SAME_SPECIES_FOLLOWUP_ORIGINAL_GPS_Q_FIRST_META_V1",
        "stage":"FIRST_NEW_ORIGINAL_SOURCE_Q_OUTCOMES_AFTER_SPECIES_DISCOVERY_CONDITIONAL_VALIDATION",
        "original_previous_archive_source_DOIs":63,
        "original_followup_source_DOIs_all":6,
        "original_source_file_admissible_pre_event":4,
        "original_pre_event_metadata_HOLD":2,
        "new_original_source_q_supported":sum(r["original_new_source_q"]["available"] for r in rows),
        "new_original_source_daylength_supported":sum(r["original_new_source_q_daylength_slope"]["available"] for r in rows),
        "new_source_effect_repeatability_nominal_verdict_counts":dict(sorted(counts.items())),
        "all_six_original_source_taxon_results_including_HOLD":rows,
        "known_prior_solar_heading_GPS_source_used":1,
        "previous_mallard_metadata_source_HOLD_unchanged":True,
        "new_source_study_IDs_selected_after_previous_same_species_results":True,
        "source_taxon_same_original_tag_population_independence_unverified":True,
        "no_original_animal_ID_GPS_coordinates_exact_dates_or_per_tag_coefficients_exported":True,
        "previous_original_48_solar_direction_signed_heldout_nonreplication":True,
        "previous_original_12_ERA5_raw_source_HOLD_unchanged":True,
        "scientific_ceilings":[
            "All four prediction directions were fixed before these original GPS q outcomes but their species targets were selected after old q results; not a completely prospective global hypothesis test.",
            "A source marked no-match in eight archived solar result files may have prior GPS uses elsewhere. The fox source definitely had earlier GPS opened for solar direction.",
            "Separate source DOI and study ID are not proof of new ecological populations or disjoint original animal tags, and telemetry observation process differs.",
            "Source-native animal-bootstrap nominal 95 CIs cannot be combined into confirmatory cross-study meta-level p values or independent source-effect CIs.",
            "q represents NIGHT−DAY original observed straight-line GPS displacement speed, not real activity fraction or cause of photoperiodic changes.",
        ]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    files=sorted(a.directory.glob("source_*.json"))
    if len(files)!=6:
        raise ValueError("All six first source attempt artifacts including Holds required")
    raw=[json.loads(z.read_text()) for z in files]
    obj=aggregate(raw)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("SIX_TARGET_SOURCE_ORIGINAL_Q_FIRST",json.dumps({
        "source_q_supported":obj["new_original_source_q_supported"],
        "source_daylength_supported":obj["new_original_source_daylength_supported"],
        "source_nominal_verdict_counts":obj["new_source_effect_repeatability_nominal_verdict_counts"],
        "source": [{"doi":z["original_source_DOI"],"taxon":z["scientific_name"],
                    "status":z["source_status"],"q":z["original_new_source_q"]["q"],
                    "q_ci":z["original_new_source_q"]["source_tag_bootstrap95"],
                    "repeatability":z["nominal_source_CI_repeatability_verdict"],
                    "prior_GPS_exposure":z["known_previous_original_GPS_exposure_status"]}
                   for z in obj["all_six_original_source_taxon_results_including_HOLD"]],
    },sort_keys=True))

if __name__=="__main__":
    main()
