"""Read-only replay of ORIGINAL 35-source GPS observation-eligibility funnels.

All first source responses are previously frozen; this is a post-outcome
methodological diagnosis. This runner never creates or replaces source q or
daylength endpoints. It only returns aggregate gate counts and ratios.
A former source-wide HOLD remains HOLD and is NOT reclassified if a provider
has since made another file available. Human-source GPS never queried.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import tempfile

from scripts.movebank_native_ecology_20261019 import HERE, original_files, checksum_and_download
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.rsds_437_systematic_source_frame_v1 import pinned
from scripts.run_rsds_437_wave1_original_gps_first_effects_v1 import freeze_original as freeze_wave1
from scripts.run_rsds_437_wave2_original_gps_first_effects_v1 import freeze_original as freeze_wave2
from scripts.rsds_35_original_gps_observation_funnel_v1 import profile_native, redact_native_replay

SPEC=HERE/"research/RSDS_35_SOURCE_GPS_OBSERVATION_FUNNEL_V1_CONTRACT.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def frozen():
    c=json.loads(SPEC.read_text(encoding="utf-8"))
    if (c["version"]!="RSDS_35_ORIGINAL_SOURCE_GPS_OBSERVATION_FUNNEL_V1"
        or c["original_35_source_study_ids"]!=35
        or c["original_registered_native_source_rows"]!=50
        or c["source_MD5_original_qualified"]!=28
        or c["original_file_or_human_privacy_HOLD"]!=7
        or c["original_GPS_source_wide_HOLD"]!=3
        or c["first_result_GPS_parsed_sources"]!=25
        or c["source_any_native_q_supported"]!=13
        or c["source_complete_native_q_supported"]!=10
        or c["previous_48_signed_sunlight_external"]!=NEGATIVE
        or c["previous_ERA5_raw_SHA_12_HOLD"] is not True):
        raise ValueError("Original observation diagnostic contract changed")
    wave1=pinned(**{"path":c["original_16_wave_source_result"]["path"],
                    "sha":c["original_16_wave_source_result"]["git_blob_sha"]})
    all35=pinned(**{"path":c["original_35_wave_source_result"]["path"],
                    "sha":c["original_35_wave_source_result"]["git_blob_sha"]})
    original1=freeze_wave1()[1]
    original2=freeze_wave2()[1]
    rows=all35["combined_source_study_summary_rows"]
    if (len(original1)!=16 or len(original2)!=19
        or len(rows)!=35
        or len(wave1["all_first_wave_registered_native_group_effects"])!=27
        or sum(z["registered_native_count"] for z in rows)!=50
        or sum(z["source_MD5_confirmed"] for z in rows)!=28
        or sum(z["complete_source_native_q"] for z in rows)!=10
        or all35["human_location_source_no_access"] is not True):
        raise ValueError("Prior 35-source q or recorded source eligibility changed")
    return c,original1,original2,rows,wave1


def replay_ref(qold,qnow,name):
    if (qold is None)!=(qnow is None):
        raise ValueError("Original source-native "+name+" estimation classification altered")
    if qold is not None:
        if (qold["n_individuals"]!=qnow["n_individuals"]
            or not math.isclose(qold["mean"],qnow["mean"],rel_tol=0,abs_tol=1e-12)):
            raise ValueError("Previously frozen native "+name+" q value changed")


def one_source(wave,index):
    if wave not in ("wave1","wave2"):
        raise ValueError("Only two source waves, no source substitution")
    c,a,b,allrows,original1=frozen()
    count=16 if wave=="wave1" else 19
    if type(index)!=int or not 0<=index<count:
        raise ValueError("Only frozen original source indices allowed")
    src=(a if wave=="wave1" else b)[index]
    source_row=allrows[index if wave=="wave1" else 16+index]
    native=(src["registered_native_labels"] if wave=="wave1"
            else src["registered_native_taxa"])
    uid=(src["repository_study_id"] if wave=="wave1"
         else src["original_repository_study_id"])
    expected_label=(src["source_DOI_label"] if wave=="wave1"
                    else src["original_DOI_labels"][0])
    if (source_row["source_original_repository_study_id"]!=uid
        or source_row["source_DOI_labels"][0]!=expected_label
        or len(native)!=source_row["registered_native_count"]):
        raise ValueError("Prior original source roster study identity drift")
    report={
      "version":"RSDS_35_ORIGINAL_GPS_SOURCE_OBSERVATION_FUNNEL_V1",
      "wave":wave,"index":index,
      "source_original_DOI":expected_label,
      "source_repository_study_id":uid,
      "original_registered_native_taxa":native,
      "prior_first_source_q_status":source_row["original_source_status"],
      "original_first_q_study_complete":source_row["complete_source_native_q"],
      "original_first_daylength_study_complete":source_row["complete_source_daylength"],
      "original_MD5_replayed":False,
      "source_status":"ARCHIVED_HOLD",
      "original_native_funnel":[],
      "source_UTC_observation_cadence_replay_status":"NOT_EVALUATED",
      "original_raw_original_event_count":None,
      "source_hold_reason":source_row["source_hold_error_type"],
      "original_48_signed_sunlight_external_NONREPLICATION_UNCHANGED":NEGATIVE,
      "original_12_raw_ERA5_source_HOLD_UNCHANGED":True,
      "no_individual_ID_GPS_location_date_timestamp_or_q_values_exported":True
    }
    archived=source_row["original_source_status"]
    if archived in (
        "PRE_EVENT_ORIGINAL_FILE_METADATA_HOLD",
        "SOURCE_FILE_METADATA_HOLD",
        "HUMAN_LOCATION_PRIVACY_HOLD",
        "ORIGINAL_SOURCE_GPS_HOLD",
        "SOURCE_GPS_HOLD"):
        if "Homo sapiens" in native:
            report["source_status"]="HUMAN_LOCATION_PRIVACY_HOLD_NO_ACCESS"
            report["source_UTC_observation_cadence_replay_status"]="NO_HUMAN_GPS_ACCESS"
        else:
            report["source_status"]="ARCHIVED_ORIGINAL_SOURCE_HOLD_UNCHANGED"
            report["source_UTC_observation_cadence_replay_status"]="NO_ADMISSIBLE_ARCHIVED_NATIVE_SOURCE"
        return report
    if "Homo sapiens" in native:
        raise ValueError("Human source original privacy HOLD was altered")
    metadata=src["original_file"]
    if metadata is None:
        raise ValueError("Parsed original source missing first provider MD5 file")
    file_name=metadata["filename"]
    suffix=(".csv.gz" if file_name.lower().endswith(".csv.gz") else
            ".zip" if file_name.lower().endswith(".zip") else ".csv")
    result_status="SOURCE_ORIGINAL_GPS_REPLAY_HOLD"
    stage="SOURCE_ORIGINAL_PROVIDER_FILE_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_private_obs_funnel_") as td:
            local=Path(td)/("original"+suffix)
            fs=[x for x in original_files(src["package_uuid"])
                if x["bitstream_uuid"]==metadata["bitstream_uuid"]]
            if (len(fs)!=1 or fs[0].get("bytes")!=metadata["bytes"]
                or fs[0]["checksum"]["checkSumAlgorithm"].upper()!="MD5"
                or fs[0]["checksum"]["value"].lower()!=metadata["checksum"]["value"]):
                raise ValueError("Original provider source MD5 identity changed")
            if checksum_and_download(fs[0],local)["provider_checksum_verified"] is not True:
                raise ValueError("Original provider MD5 bytes mismatch")
            report["original_MD5_replayed"]=True
            stage="ORIGINAL_HOURLY_GPS_AND_SOURCE_Q_REPLAY"
            animals,coverage,source_metadata=load_source_taxa(local,file_name,native)
            oldrows=(original1["all_first_wave_registered_native_group_effects"]
                     if wave=="wave1" else None)
            profiles=[]
            native_q=0
            native_L=0
            native_qmeans=[]
            native_Lmeans=[]
            for name in native:
                one=profile_native(name,animals[name])
                d=one["replayed_archived_native_q_support"]
                if oldrows is not None:
                    source_match=[r for r in oldrows
                                  if r["source_study_DOI"]==expected_label
                                  and r["native_taxon"]==name]
                    if len(source_match)!=1:
                        raise ValueError("First original wave1 registered native row not found")
                    replay_ref(source_match[0]["q"], {
                        "n_individuals":d["n_tag_q"],
                        "mean":d["source_q_mean_for_private_replay_only"]
                    } if d["q_estimable"] else None,"q")
                    replay_ref(source_match[0]["daylength"], {
                        "n_individuals":d["n_tag_daylength"],
                        "mean":d["source_daylength_mean_for_private_replay_only"]
                    } if d["daylength_estimable"] else None,"daylength")
                native_q+=int(d["q_estimable"])
                native_L+=int(d["daylength_estimable"])
                if d["q_estimable"]:native_qmeans.append(d["source_q_mean_for_private_replay_only"])
                if d["daylength_estimable"]:native_Lmeans.append(d["source_daylength_mean_for_private_replay_only"])
                profiles.append(redact_native_replay(one))
            if (native_q!=source_row["native_q_supported"]
                or native_L!=source_row["native_daylength_supported"]
                or source_metadata["original_total_rows"]!=source_row["original_source_event_row_count"]):
                raise ValueError("Archived source support and GPS row count changed")
            if source_row["complete_source_native_q"]:
                if (native_q!=len(native) or
                    not math.isclose(sum(native_qmeans)/len(native),source_row["source_q_mean_equal_registered_native"],rel_tol=0,abs_tol=1e-12)):
                    raise ValueError("Frozen source-full q mean changed")
            if source_row["complete_source_daylength"]:
                if (native_L!=len(native) or not math.isclose(sum(native_Lmeans)/len(native),
                       source_row["source_daylength_mean_equal_registered_native"],
                       rel_tol=0,abs_tol=1e-12)):
                    raise ValueError("Frozen source-full q~daylength mean changed")
            report["original_native_funnel"]=profiles
            report["source_status"]="REPLAYED_ARCHIVED_NATIVE_GPS_WITH_OBSERVATION_FUNNEL"
            report["source_UTC_observation_cadence_replay_status"]="REPLAY_VERIFIED"
            report["original_raw_original_event_count"]=source_metadata["original_total_rows"]
            report["original_native_q_eligible_count"]=native_q
            report["original_native_daylength_eligible_count"]=native_L
            report["source_hold_reason"]=None
            animals.clear()
    except Exception as exc:
        # Stage failure is a meaningful source-replay discrepancy, not
        # permission to replace the sample or revise the archived outcome.
        report["source_status"]=result_status
        report["source_UTC_observation_cadence_replay_status"]="REPLAY_DIAGNOSTIC_HOLD"
        report["original_native_funnel"]=[]
        report["source_hold_reason"]="REPLAY_"+stage+"_"+type(exc).__name__
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--wave",choices=("wave1","wave2"),required=True)
    p.add_argument("--index",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=one_source(a.wave,a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                     encoding="utf-8")
    print("RSDS_ORIGINAL_GPS_OBSERVATION_FUNNEL",json.dumps({
       "wave":a.wave,"index":a.index,"status":result["source_status"],
       "md5":result["original_MD5_replayed"],
       "native":len(result["original_registered_native_taxa"]),
       "replayed_native":len(result["original_native_funnel"]),
       "first_complete_source_q":result["original_first_q_study_complete"],
       "hold":result["source_hold_reason"],
    },sort_keys=True),flush=True)
    if result["source_status"]=="SOURCE_ORIGINAL_GPS_REPLAY_HOLD":
        raise SystemExit(2)

if __name__=="__main__":
    main()
