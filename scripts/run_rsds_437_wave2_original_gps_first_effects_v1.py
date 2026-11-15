"""First GPS outcomes on nineteen previous-SHA-ranked new source slots.

All original repository study IDs and 15 source file MD5 identities were
committed before any Wave2 animal GPS event was accessed. The 3 file-size
HOLDS and 1 Homo sapiens privacy HOLD remain in the denominator. Entire
registered native species lists are never pruned based on response signs.

Outputs are source-level groups only: no tracks, dates, coordinates, private
animal names/IDs or per-animal vectors. No cross-source meta claims here.
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
from scripts.run_rsds_eight_study_gps_cadence_v1 import per_tag_year_days
from scripts.run_rsds_four_new_multinative_original_gps_temporal_v1 import (
    source_reference_latitude, common_original_calendar
)
from scripts.rsds_437_wave1_joint_source_covariance_v1 import original_joint
from scripts.rsds_437_systematic_source_frame_v1 import pinned
from scripts.rsds_437_wave2_original_file_metadata_v1 import freeze as source_selection

SPEC=HERE/"research/RSDS_437_WAVE2_ORIGINAL_GPS_FIRST_EFFECTS_V1_CONTRACT.json"
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"

def freeze_original():
    c=json.loads(SPEC.read_text(encoding="utf-8"))
    if (c.get("version")!="RSDS_437_WAVE2_19_ORIGINAL_GPS_FIRST_EFFECTS_CONTRACT_V1"
        or c.get("wave2_source_slots")!=19
        or c.get("wave2_file_metadata_admissible")!=15
        or c.get("wave2_file_size_HOLD")!=3
        or c.get("wave2_human_privacy_HOLD")!=1
        or c.get("original_file_max_bytes")!=200_000_000
        or c.get("original_read_max_CSV_rows")!=2_000_000
        or c.get("new_wave2_GPS_event_rows_opened_prior_to_contract")!=0
        or c.get("historical_48_signed_irradiance_external")!=NEG
        or c.get("historical_12_ERA5_original_raw_SHA_HOLD") is not True):
        raise ValueError("Wave2 source/GPS contract altered")
    prior=pinned(c["wave1_original_first_result"],c["wave1_first_git_blob"])
    if (prior["original_selection_wave_slots"]!=16 or
        prior["stratified_source_complete_effect_sensitivity"]["q"]["first_wave_source_complete_response"]!=4):
        raise ValueError("Wave1 first source effect immutable frame altered")
    result=pinned(c["source_provider_file_metadata_path"],
                  c["source_provider_file_metadata_git_blob"])
    chosen=source_selection()
    sources=result["source_rows"]
    if (len(chosen)!=19 or len(sources)!=19
        or result["source_status_counts"]!={
            "HUMAN_LOCATION_PRIVACY_HOLD":1,
            "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":15,
            "SOURCE_METADATA_HOLD":3}
        or result["human_GPS_or_position_data_retrieved"] is not False
        or result["all_new_wave2_animal_GPS_events_opened"]!=0
        or result["wave2_source_selection_git_blob"]!=c["source_selection_git_blob"]):
        raise ValueError("Original provider/source event-free first receipt altered")
    for a,b in zip(chosen,sources):
        if (a["wave2_index"]!=b["wave2_index"] or
            a["package_uuid"]!=b["package_uuid"] or
            a["original_repository_study_id"]!=b["original_repository_study_id"] or
            a["original_DOI_labels"]!=b["original_DOI_labels"] or
            a["scientific_name_labels"]!=b["registered_native_taxa"] or
            a["stratum"]!=b["selection_stratum"]):
            raise ValueError("Wave2 native source metadata changed after SHA source commitment")
    return c,sources

def process(index):
    _, sources=freeze_original()
    if type(index)!=int or not 0<=index<19:
        raise ValueError("Only nineteen frozen wave2 source slots")
    source=sources[index]
    taxa=source["registered_native_taxa"]
    result={
      "version":"RSDS_WAVE2_ORIGINAL_GPS_FIRST_SOURCE_Q_DAYLENGTH_V1",
      "wave2_index":index,
      "source_package_uuid":source["package_uuid"],
      "source_original_repository_study_id":source["original_repository_study_id"],
      "source_original_DOI_labels":source["original_DOI_labels"],
      "source_metadata_stratum":source["selection_stratum"],
      "original_registered_native_taxa":taxa,
      "original_file_source_status":source["source_status"],
      "original_file_MD5_bytes_verified":False,
      "source_science_status":"SOURCE_GPS_HOLD",
      "source_hold_stage":None,
      "source_hold_reason":None,
      "native_group_results":[{"taxon":x,"q":None,"daylength":None,
                               "status":"SOURCE_UNEVALUATED"} for x in taxa],
      "native_joint_calendar_covariance":{"status":"NOT_EVALUATED"},
      "original_source_GPS_event_rows_opened":0,
      "human_mobility_GPS_data_opened":False,
      "old_48_signed_solar_external_nonreplication_UNCHANGED":NEG,
      "old_12_ERA5_source_raw_SHA_HOLD_UNCHANGED":True,
      "original_raw_GPS_dates_tag_IDs_coordinates_exported":False,
    }
    # The human source stops before any provider access. Never substitute it.
    if "Homo sapiens" in taxa:
        if source["source_status"]!="HUMAN_LOCATION_PRIVACY_HOLD":
            raise ValueError("Human tracking privacy HOLD was reclassified")
        result["source_science_status"]="HUMAN_LOCATION_PRIVACY_HOLD"
        result["source_hold_stage"]="HUMAN_SOURCE_PRIVACY_NO_ACCESS"
        result["source_hold_reason"]="NO_HUMAN_TRAJECTORY_SOURCE_ACCESS"
        return result
    if source["source_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":
        result["source_science_status"]="SOURCE_FILE_METADATA_HOLD"
        result["source_hold_stage"]="SOURCE_ONLY_ORIGINAL_FILE_METADATA"
        result["source_hold_reason"]=source.get("source_hold_reason")
        return result

    meta=source["original_file"]
    ext=(".zip" if meta["filename"].lower().endswith(".zip") else
         ".csv.gz" if meta["filename"].lower().endswith(".csv.gz") else ".csv")
    stage="ORIGINAL_UUID_BYTES_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_wave2_private_gps_") as directory:
            filepath=Path(directory)/("original"+ext)
            candidates=[f for f in original_files(source["package_uuid"])
                        if f["bitstream_uuid"]==meta["bitstream_uuid"]]
            if (len(candidates)!=1 or candidates[0]["bytes"]!=meta["bytes"]
                or candidates[0].get("checksum",{}).get("checkSumAlgorithm","").upper()!="MD5"
                or candidates[0]["checksum"]["value"].lower()!=meta["checksum"]["value"]):
                raise ValueError("Original provider bitstream identity changed")
            verified=checksum_and_download(candidates[0],filepath)
            if verified.get("provider_checksum_verified") is not True:
                raise ValueError("Provider MD5 mismatch on downloaded original")
            result["original_file_MD5_bytes_verified"]=True

            stage="ORIGINAL_NATIVE_GPS_Q_AND_DAYLENGTH"
            tracks,coverage,report=load_source_taxa(filepath,meta["filename"],taxa)
            result["original_source_GPS_event_rows_opened"]=report["original_total_rows"]
            grouped=[]
            for sp in taxa:
                data,_=one_taxon(sp,tracks[sp])
                grouped.append({
                  "taxon":sp,
                  "original_animal_count":data["total_original_animals"],
                  "original_hour_selected_fixes":coverage[sp]["hour_selected_fixes"],
                  "original_adjacent_intervals":data["original_adjacent_hour_steps"],
                  "qualified_same_daypart_adjacent_intervals":data["qualified_daypart_step_observations"],
                  "n_q_tagged_animals":data["admitted_animals_diel"],
                  "n_daylength_tagged_animals":data["admitted_animals_photoperiod_slope"],
                  "q":data["day_night_speed_contrast"],
                  "daylength":data["daylength_association"],
                  "status":"NATIVE_Q_ESTIMABLE" if data["day_night_speed_contrast"] is not None
                          else "NATIVE_Q_UNAVAILABLE"
                })
            result["native_group_results"]=grouped
            result["source_science_status"]=(
                "SOURCE_AT_LEAST_ONE_NATIVE_Q_ESTIMABLE" if any(x["q"] is not None for x in grouped)
                else "SOURCE_NO_NATIVE_Q_ESTIMABLE"
            )
            if len(taxa)>1:
                stage="REGISTERED_ALL_NATIVE_SAME_CALENDAR_COVARIANCE"
                try:
                    original_days={
                        sp:per_tag_year_days(tracks[sp],"ORIGINAL_1H")[0]
                        for sp in taxa
                    }
                    aligned,dates,rounds=common_original_calendar(original_days)
                    reference=source_reference_latitude(tracks)
                    if reference is None:
                        raise ValueError("Missing original source astronomical reference")
                    result_cov=original_joint(aligned,dates,taxa,reference,
                                              seed=20261115+index*7001,replicates=1000)
                    result_cov["source_calendar_support_fixed_point_rounds"]=rounds
                    result["native_joint_calendar_covariance"]=result_cov
                except Exception as exc:
                    result["native_joint_calendar_covariance"]={
                        "status":"ALL_NATIVE_COMMON_CALENDAR_COVARIANCE_HOLD",
                        "hold_error_type":type(exc).__name__,
                        "all_original_native_taxa_retained":True
                    }
            else:
                result["native_joint_calendar_covariance"]={
                    "status":"NOT_APPLICABLE_SINGLE_NATIVE_SOURCE"}
            result["source_hold_stage"]=None
            result["source_hold_reason"]=None
            tracks.clear()
    except Exception as exc:
        result["source_science_status"]="SOURCE_GPS_HOLD"
        result["source_hold_stage"]=stage
        result["source_hold_reason"]=type(exc).__name__
        # Original source error can invalidate all previously partly parsed
        # native effects: fail closed without selecting a favorable subset.
        result["native_group_results"]=[{
            "taxon":sp,"q":None,"daylength":None,
            "status":"SOURCE_WIDE_FAIL_CLOSED"} for sp in taxa]
        result["native_joint_calendar_covariance"]={"status":"SOURCE_WIDE_HOLD"}
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--index",type=int,required=True)
    parser.add_argument("--out",type=Path,required=True)
    arg=parser.parse_args()
    outcome=process(arg.index)
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(outcome,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                       encoding="utf-8")
    print("RSDS_WAVE2_GPS_ORIGINAL_SOURCE",json.dumps({
      "index":arg.index,
      "source_status":outcome["source_science_status"],
      "native_count":len(outcome["native_group_results"]),
      "native_q_supported":sum(x["q"] is not None for x in outcome["native_group_results"]),
      "source_MD5_verified":outcome["original_file_MD5_bytes_verified"],
      "joint_covariance_status":outcome["native_joint_calendar_covariance"]["status"],
      "original_event_rows_opened":outcome["original_source_GPS_event_rows_opened"],
      "hold_stage":outcome["source_hold_stage"],
      "hold_type":outcome["source_hold_reason"]
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
