"""Outcome-blind selected 16-source Wave1 ORIGINAL GPS first effects.

One source package per original, previously frozen 437-package SHA ranking.
Three provider metadata/size HOLD slots stay in the denominator, never
replaced. For thirteen source-file candidates, original UUID, byte count and
MD5 must be replayed before original GPS opening. All registered native
taxa are processed using the same archived 1h, 0.5-2h observed q method.

No private original animal IDs, positions, individual q, event times/dates,
or GPS records are returned or written; only DOI/native group effects and
optional joint source-native month-block covariance matrices.
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
from scripts.rsds_437_systematic_source_frame_v1 import pinned
from scripts.rsds_437_wave1_original_file_metadata_v1 import (
    freeze as prior_source_wave, SELECTION_BLOB_SHA, git_blob
)
from scripts.rsds_437_wave1_joint_source_covariance_v1 import original_joint

SPEC=HERE/"research/RSDS_437_WAVE1_ORIGINAL_GPS_FIRST_EFFECTS_V1_CONTRACT.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"

def freeze_original():
    contract=json.loads(SPEC.read_text(encoding="utf-8"))
    if (contract["version"]!="RSDS_437_WAVE1_16_SOURCE_ORIGINAL_GPS_FIRST_EFFECTS_V1"
        or contract["wave1_source_slot_count"]!=16
        or contract["metadata_file_admissible"]!=13
        or contract["metadata_file_holds"]!=3
        or contract["original_file_max_bytes"]!=200_000_000
        or contract["original_event_row_max"]!=2_000_000
        or contract["old_48_radiation_external_nonreplication"]!=NEGATIVE
        or contract["old_12_ERA5_raw_SHA_HOLD"] is not True
        or contract["selection_manifest_git_blob"]!=SELECTION_BLOB_SHA):
        raise ValueError("Original wave1 pre-event source/effect contract modified")
    first=pinned(contract["first_file_preflight"],
                 contract["first_file_preflight_git_blob"])
    _,selection=prior_source_wave()
    source=first["provider_original_source_rows"]
    if (len(source)!=16 or first["original_package_slots"]!=16
        or first["original_GPS_event_records_opened"]!=0
        or first["source_status_counts"]!={
            "ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":13,
            "SOURCE_METADATA_HOLD":3}
        or not first["prior_48_shortwave_direction_external_not_replicated"]
        or not first["prior_original_ERA5_12_source_HOLD"]):
        raise ValueError("First provider original metadata hold and count altered")
    if sum(r["source_status"]=="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE"
           for r in source)!=13:
        raise ValueError("Frozen source provider admission slots changed")
    for a,b in zip(selection,source):
        if (a["ordinal"]!=b["ordinal"] or a["package_uuid"]!=b["package_uuid"]
            or a["original_repository_study_id"]!=b["repository_study_id"]
            or a["scientific_names"]!=b["registered_native_labels"]
            or a["source_DOI_label"]!=b["source_DOI_label"]):
            raise ValueError("Source IDs or native taxonomy changed since freeze")
    return contract,source

def source_effects(index):
    contract,sources=freeze_original()
    if not isinstance(index,int) or index not in range(16):
        raise ValueError("Only the sixteen preselected original source slots")
    src=sources[index]
    report={
        "version":"RSDS_437_WAVE1_ORIGINAL_GPS_FIRST_SOURCE_Q_DAYLENGTH_V1",
        "index":index,"original_ordinal":src["ordinal"],
        "source_package_uuid":src["package_uuid"],
        "original_repository_study_id":src["repository_study_id"],
        "source_DOI_label":src["source_DOI_label"],
        "source_selection_stratum":src["selection_stratum"],
        "registered_native_taxa":src["registered_native_labels"],
        "original_file_preflight_status":src["source_status"],
        "original_provider_MD5_validated":False,
        "source_status":"ORIGINAL_SOURCE_GPS_HOLD",
        "source_hold_stage":None,
        "source_hold_reason":None,
        "native_effects":[{"native_taxon":n,"q":None,"daylength":None,
                           "reason":"SOURCE_NOT_EVALUATED"}
                          for n in src["registered_native_labels"]],
        "native_common_calendar_covariance":{"status":"NOT_EVALUATED"},
        "observed_original_source_event_rows_if_opened":0,
        "old_48_signed_solar_shortwave_external_NONREPLICATED_UNCHANGED":NEGATIVE,
        "old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED":True,
        "no_private_GPS_or_individual_data_in_output":True,
    }
    if src["source_status"]!="ORIGINAL_SOURCE_FILE_METADATA_ADMISSIBLE":
        report["source_status"]="PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD"
        report["source_hold_stage"]="SOURCE_ONLY_ORIGINAL_FILE_METADATA"
        report["source_hold_reason"]=src.get("hold_reason")
        return report

    original=src["original_file"]
    suffix=(".zip" if original["filename"].lower().endswith(".zip")
            else ".csv.gz" if original["filename"].lower().endswith(".csv.gz")
            else ".csv")
    stage="PROVIDER_ORIGINAL_FILE_UUID_BYTES_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_original_wave1_private_") as directory:
            local=Path(directory)/("original"+suffix)
            f=[x for x in original_files(src["package_uuid"])
               if x.get("bitstream_uuid")==original["bitstream_uuid"]]
            if (len(f)!=1 or f[0]["bytes"]!=original["bytes"]
                or f[0].get("checksum",{}).get("checkSumAlgorithm","").upper()!="MD5"
                or f[0]["checksum"]["value"].lower()!=original["checksum"]["value"]):
                raise ValueError("Provider original file bytes/UUID/MD5 changed")
            if not checksum_and_download(f[0],local)["provider_checksum_verified"]:
                raise ValueError("Original provider MD5 bytes not verified")
            report["original_provider_MD5_validated"]=True

            stage="ORIGINAL_GPS_SOURCE_SCHEMA_NATIVE_AND_Q"
            native=src["registered_native_labels"]
            tracks,coverage,source_count=load_source_taxa(
                local,original["filename"],native
            )
            report["observed_original_source_event_rows_if_opened"]=source_count[
                "original_total_rows"]
            original_effects=[]
            for taxon in native:
                effect,_=one_taxon(taxon,tracks[taxon])
                original_effects.append({
                    "native_taxon":taxon,
                    "original_source_animal_count":effect["total_original_animals"],
                    "original_adjacent_hour_steps":effect["original_adjacent_hour_steps"],
                    "original_qualified_daypart_steps":effect["qualified_daypart_step_observations"],
                    "original_tagged_animals_with_qualified_q":effect["admitted_animals_diel"],
                    "original_tagged_animals_with_daylength_slope":effect[
                        "admitted_animals_photoperiod_slope"],
                    "source_native_hour_selected_GPS_fixes":coverage[taxon]["hour_selected_fixes"],
                    "q":effect["day_night_speed_contrast"],
                    "daylength":effect["daylength_association"],
                })
            report["native_effects"]=original_effects
            supported=sum(x["q"] is not None for x in original_effects)
            report["source_status"]=(
                "ORIGINAL_SOURCE_GPS_Q_ESTIMABLE" if supported
                else "ORIGINAL_SOURCE_GPS_Q_UNAVAILABLE_NOT_ZERO"
            )

            if len(native)>1:
                stage="ORIGINAL_NATIVE_COMMON_CALENDAR_JOINT_COVARIANCE"
                try:
                    daymaps={
                      taxon:per_tag_year_days(tracks[taxon],"ORIGINAL_1H")[0]
                      for taxon in native
                    }
                    aligned,dates,rounds=common_original_calendar(daymaps)
                    latitude=source_reference_latitude(tracks)
                    if latitude is None:
                        raise ValueError("Missing source astronomy reference")
                    result=original_joint(aligned,dates,native,latitude,
                                          seed=20261115+index*5001,replicates=1000)
                    result["source_common_calendar_fixed_point_rounds"]=rounds
                    report["native_common_calendar_covariance"]=result
                except Exception as error:
                    report["native_common_calendar_covariance"]={
                        "status":"ORIGINAL_NATIVE_COMMON_CALENDAR_COVARIANCE_HOLD",
                        "hold_error_type":type(error).__name__,
                        "all_registered_native_taxa_retained":True,
                    }
            else:
                report["native_common_calendar_covariance"]={
                    "status":"NOT_APPLICABLE_SINGLE_REGISTERED_NATIVE"}
            report["source_hold_stage"]=None
            report["source_hold_reason"]=None
            tracks.clear()
    except Exception as error:
        # A source-wide failure voids all partially observed source-native
        # coefficients. Do not promote a preferred species/subset or replace
        # the source after seeing data.
        report["source_status"]="ORIGINAL_SOURCE_GPS_HOLD"
        report["source_hold_stage"]=stage
        report["source_hold_reason"]=type(error).__name__
        report["native_effects"]=[
            {"native_taxon":n,"q":None,"daylength":None,
             "reason":"SOURCE_WIDE_FAIL_CLOSED"}
            for n in src["registered_native_labels"]]
        report["native_common_calendar_covariance"]={"status":"SOURCE_WIDE_HOLD"}
    return report

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--index",type=int,required=True)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    result=source_effects(args.index)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                        encoding="utf-8")
    print("RSDS_ORIGINAL_437_WAVE1_SOURCE",json.dumps({
        "source_index":result["index"],
        "DOI":result["source_DOI_label"],
        "status":result["source_status"],
        "MD5":result["original_provider_MD5_validated"],
        "native_rows":len(result["native_effects"]),
        "estimable_q_native":sum(x["q"] is not None for x in result["native_effects"]),
        "native_joint_covariance":result["native_common_calendar_covariance"]["status"],
        "HOLD":result["source_hold_stage"],
        "reason":result["source_hold_reason"],
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
