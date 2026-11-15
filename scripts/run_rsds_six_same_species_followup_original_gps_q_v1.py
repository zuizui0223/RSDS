"""First source-native GPS NIGHT-DAY q on six after-discovery fixed study slots.

Existing source-file/HOLD receipts are SHA pinned before event reading.
Only original provider MD5-verified original animal GPS bytes are processed
ephemerally on GitHub Action runners. Output native source q, daylight slope,
tag bootstrap CI and source-specific HOLD, never original tags/coordinates.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import tempfile

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download
)
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.rsds_437_systematic_source_frame_v1 import pinned

SPEC=HERE/"research/RSDS_SIX_FOLLOWUP_ORIGINAL_GPS_Q_FIRST_V1_CONTRACT.json"
VERSION="RSDS_SIX_SAME_SPECIES_TARGETED_FOLLOWUP_ORIGINAL_GPS_Q_V1_CONTRACT"
RESULT_VERSION="RSDS_SIX_FOLLOWUP_ORIGINAL_SOURCE_NATIVE_Q_FIRST_V1"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"

def freeze_original():
    c=json.loads(SPEC.read_text())
    if (c.get("version")!=VERSION or
        c.get("source_candidates_original")!=6 or
        c.get("source_original_file_MD5_admissible")!=4 or
        c.get("source_original_file_HOLD")!=1 or
        c.get("prior_original_source_metadata_HOLD_preserved")!=1 or
        c.get("original_source_file_max_bytes")!=200_000_000 or
        c.get("original_read_max_CSV_rows")!=2_000_000 or
        c.get("source_gps_location_data_currently_opened_for_four_q_targets")!=0 or
        c.get("source_holds_cannot_be_replaced_or_modified") is not True or
        c.get("first_target_species_selected_after_original_63_source_outcomes") is not True or
        c.get("source_study_independence_not_proven") is not True or
        c.get("no_raw_animal_ID_coordinates_exact_dates_tracks_or_per_tag_coefficients_exported") is not True or
        c.get("original_48_signed_irradiance_heading_external_nonreplicated") is not True or
        c.get("original_12_ERA5_source_SHA_HOLD") is not True):
        raise ValueError("Six fixed source follow-up q/GPS contract modified")
    fixed_selection=pinned(c["selection_receipt_path"],c["selection_receipt_git_blob"])
    exposure=pinned(c["prior_source_exposure_receipt_path"],
                    c["prior_source_exposure_receipt_git_blob"])
    metadata=pinned(c["original_provider_file_receipt_path"],
                    c["original_provider_file_receipt_git_blob"])
    prior=pinned(c["original_previous_63_source_native_atlas_path"],
                 c["original_previous_63_source_native_atlas_git_blob"])
    sources=metadata["original_source_metadata_rows"]
    if (len(sources)!=6 or
        len(fixed_selection["original_target_metadata_candidates_eligible"])!=6 or
        len(exposure["original_followup_source_exposure_rows"])!=6 or
        prior["distinct_original_archived_source_DOIs_and_study_IDs"]!=63 or
        metadata["source_original_metadata_status_counts"]!={"ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":4,
                "PREVIOUS_SOURCE_METADATA_HOLD_PRESERVED":1,
                "SOURCE_METADATA_HOLD":1} or
        metadata["all_new_original_animal_GPS_events_opened"]!=0 or
        not metadata["prior_fox_original_GPS_read_known"] or
        not metadata["prior_mallard_source_metadata_HOLD_preserved"]):
        raise ValueError("Original six-source file metadata/HOLD receipt changed")
    expected={
        "10.5441/001/1.46b0mq21":"NEGATIVE",
        "10.5441/001/1.708":"NEGATIVE",
        "10.5441/001/1.72hh609t":"POSITIVE",
        "10.5441/001/1.ps244r11":"NEGATIVE",
    }
    if c["source_q_primary_prediction_sign_by_DOI"]!=expected:
        raise ValueError("Original pre-event scientific prediction map altered")
    for s,p,e in zip(sources,fixed_selection["original_target_metadata_candidates_eligible"],
                     exposure["original_followup_source_exposure_rows"]):
        if (s["original_source_DOI"]!=p["original_DOI"] or
            s["original_source_DOI"]!=e["original_DOI"] or
            s["original_repository_study_id"]!=p["original_repository_study_id"] or
            s["original_public_package_uuid"]!=p["original_public_package_uuid"] or
            s["registered_original_native_labels"]!=p["original_scientific_names"] or
            s["prior_original_source_exposure_status"]!=e["prior_data_exposure_status"] or
            len(s["registered_original_native_labels"])!=1):
            raise ValueError("Original source DOI/study/taxon/previous exposure changed")
    return c,sources

def process(index):
    c,sources=freeze_original()
    if type(index) is not int or not 0<=index<6:
        raise ValueError("Only six original fixed target source slots")
    s=sources[index]
    original_taxa=s["registered_original_native_labels"]
    doi=s["original_source_DOI"]
    predicted=c["source_q_primary_prediction_sign_by_DOI"].get(doi)
    result={
        "version":RESULT_VERSION,
        "source_index":index,
        "source_original_DOI":doi,
        "source_original_repository_study_id":s["original_repository_study_id"],
        "source_original_public_package_uuid":s["original_public_package_uuid"],
        "source_original_registered_native_taxa":original_taxa,
        "source_original_file_metadata_status":s["source_metadata_status"],
        "prior_original_source_exposure_status":s["prior_original_source_exposure_status"],
        "prior_original_archive_selected_prediction_sign":predicted,
        "original_file_MD5_bytes_verified":False,
        "source_q_science_status":"SOURCE_GPS_HOLD",
        "source_hold_stage":None,
        "source_hold_error_type":None,
        "native_group_results":[{"taxon":x,"q":None,"daylength":None,"status":"SOURCE_UNEVALUATED"} for x in original_taxa],
        "original_animal_GPS_event_rows_opened":0,
        "no_raw_animal_GPS_IDs_coordinates_exact_dates_or_tag_effects_exported":True,
        "earlier_original_48_solar_direction_signed_external_not_replicated":True,
        "earlier_original_12_ERA5_source_SHA_HOLD":True
    }
    if "Homo sapiens" in original_taxa:
        result["source_q_science_status"]="HUMAN_LOCATION_PRIVACY_HOLD"
        result["source_hold_stage"]="HUMAN_TRAJECTORY_NO_ACCESS"
        return result
    if s["source_metadata_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":
        result["source_q_science_status"]="SOURCE_FILE_METADATA_HOLD"
        result["source_hold_stage"]="ORIGINAL_PROVIDER_PRE_EVENT_METADATA"
        result["source_hold_error_type"]=s["source_metadata_hold_reason"]
        return result
    meta=s["original_file"]
    extension=".zip" if meta["filename"].lower().endswith(".zip") else (
        ".csv.gz" if meta["filename"].lower().endswith(".csv.gz") else ".csv")
    stage="ORIGINAL_BITSTREAM_UUID_BYTES_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_six_original_private_") as directory:
            filepath=Path(directory)/("source_original"+extension)
            candidate=[f for f in original_files(s["original_public_package_uuid"])
                       if f["bitstream_uuid"]==meta["bitstream_uuid"]]
            if (len(candidate)!=1 or
                candidate[0]["bytes"]!=meta["bytes"] or
                candidate[0].get("checksum",{}).get("checkSumAlgorithm","").upper()!="MD5" or
                candidate[0]["checksum"]["value"].lower()!=meta["checksum"]["value"]):
                raise ValueError("Pre-outcome source original bitstream UUID MD5 no longer matches")
            original=checksum_and_download(candidate[0],filepath)
            if original.get("provider_checksum_verified") is not True:
                raise ValueError("Original provider MD5 not verified")
            result["original_file_MD5_bytes_verified"]=True
            stage="ORIGINAL_GPS_SOURCE_NATIVE_Q_AND_ASTRONOMICAL_DAYLENGTH"
            tracks,coverage,report=load_source_taxa(filepath,meta["filename"],original_taxa)
            result["original_animal_GPS_event_rows_opened"]=report["original_total_rows"]
            native=[]
            for taxon in original_taxa:
                d,_=one_taxon(taxon,tracks[taxon])
                native.append({
                    "taxon":taxon,
                    "original_animal_count":d["total_original_animals"],
                    "original_hour_selected_fixes":coverage[taxon]["hour_selected_fixes"],
                    "original_adjacent_intervals":d["original_adjacent_hour_steps"],
                    "original_qualified_daypart_adjacent_intervals":d["qualified_daypart_step_observations"],
                    "original_q_tagged_animals":d["admitted_animals_diel"],
                    "original_daylength_tagged_animals":d["admitted_animals_photoperiod_slope"],
                    "q":d["day_night_speed_contrast"],
                    "daylength":d["daylength_association"],
                    "status":"NATIVE_Q_ESTIMABLE" if d["day_night_speed_contrast"] is not None
                        else "NATIVE_Q_UNAVAILABLE",
                })
            result["native_group_results"]=native
            result["source_q_science_status"]="SOURCE_Q_ESTIMABLE" if any(x["q"] is not None for x in native) else "SOURCE_Q_UNAVAILABLE_NOT_ZERO"
            tracks.clear()
    except Exception as exc:
        result["source_q_science_status"]="SOURCE_GPS_HOLD"
        result["source_hold_stage"]=stage
        result["source_hold_error_type"]=type(exc).__name__
        result["native_group_results"]=[{
            "taxon":sp,"q":None,"daylength":None,"status":"SOURCE_WIDE_FAIL_CLOSED"
        } for sp in original_taxa]
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--index",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z=process(a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("SIX_FOLLOWUP_FIRST_Q_ORIGINAL_SOURCE",json.dumps({
        "i":z["source_index"],"doi":z["source_original_DOI"],
        "status":z["source_q_science_status"],
        "native_q_supported":sum(x["q"] is not None for x in z["native_group_results"]),
        "source_MD5_verified":z["original_file_MD5_bytes_verified"],
        "original_event_rows_opened":z["original_animal_GPS_event_rows_opened"],
        "prior_GPS_exposure":z["prior_original_source_exposure_status"],
        "hold_type":z["source_hold_error_type"],
    },sort_keys=True))

if __name__=="__main__":
    main()
