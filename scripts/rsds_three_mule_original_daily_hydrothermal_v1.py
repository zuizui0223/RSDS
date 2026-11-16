"""Three original mule deer source-native GPS q × daily climate, source-first.

All 3 source DOI/native tags were frozen after prior same-species patterns,
BEFORE this new daily 1-degree ERA5 retrieval. Replays original q and
astronomical daylength n/mean before weather. Missing provider or model grids
remain explicit source/weather HOLD; no selection of favorable modelled dates.
Results are descriptive observational associations, not causal daylight or heat.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import tempfile

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned,git_blob
from scripts.rsds_five_source_native_tag_month_conditioned_v1 import frozen as old_frozen,replay_original_baseline
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_native_ecology_20261019 import original_files,checksum_and_download
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.movebank_wyoming_individual_era5_20261023 import source_admitted_daily
from scripts.movebank_wyoming_daily_climate_20261023 import daily_source_era5,adjusted_daily_effects

SPEC=HERE/"research/RSDS_THREE_MULE_ORIGINAL_DAILY_HYDROTHERMAL_CONTRACT_20261116.json"
VERSION="RSDS_THREE_MULE_ORIGINAL_DAILY_HYDROTHERMAL_V1"
SOURCE_VERSION="RSDS_THREE_MULE_ORIGINAL_SOURCE_DAILY_ERA5_GPS_FIRST_V1"
NEG="NOT_REPLICATED"
OLD_HOLD="HOLD"

def frozen():
    c=json.loads(SPEC.read_text(encoding="utf-8"))
    checks={
        "version":VERSION,
        "three_source_keys":["mule_old_utah","mule_old_wyoming","mule_new_sublette"],
        "three_source_dois":["10.5441/001/1.711","10.5441/001/1.703","10.5441/001/1.708"],
        "target_native":"Odocoileus hemionus",
        "max_original_GPS_bytes":200000000,
        "max_new_source_weather_cells":160,
        "maximum_source_weather_request_days":3653,
        "prior_original_five_source_manifest_blob":"7a25beccbb097b266d171edfaec8945e7f7d524f",
        "prior_original_month_result_blob":"27643d9c69ca6683edde5ee75945f6a4fb0eb413",
        "prior_original_matched_month_result_blob":"1a3b57bb928ecc6169e3563f83b777d2980d63ef",
        "prior_2023_wyoming_daily_era5_source_result_blob":"104b2d62e6aacf7ccfc41d5a6d866832caaadab6",
        "historical_48_source_signed_sun_azimuth_irradiance_replication":NEG,
        "historical_12_ERA5_source_raw_SHA":OLD_HOLD,
        "no_original_trajectory_tag_coordinates_dates_or_individual_slopes_exported":True,
    }
    if any(c.get(k)!=v for k,v in checks.items()):
        raise ValueError("Frozen 3-source q meteorology declaration modified")
    previous_manifest=pinned(c["prior_original_five_source_manifest_path"],
                             c["prior_original_five_source_manifest_blob"])
    previous_month=pinned(c["prior_original_month_result_path"],
                          c["prior_original_month_result_blob"])
    previous_match=pinned(c["prior_original_matched_month_result_path"],
                          c["prior_original_matched_month_result_blob"])
    prior_weather=pinned(c["prior_2023_wyoming_daily_era5_source_result_path"],
                         c["prior_2023_wyoming_daily_era5_source_result_blob"])
    for key in ("prior_wyoming_diel_climate_code","prior_source_admission_code"):
        b=(HERE/c[key+"_path"]).read_bytes()
        if git_blob(b)!=c[key+"_blob"]:raise ValueError("Source weather code version changed: "+key)
    assert previous_month["source_status_count"]=={"MONTH_NATIVE_ESTIMABLE":5}
    assert previous_match["source_status_counts"]=={"MATCHED_SUPPORT_PAIRED_SOURCE_ESTIMABLE":5}
    assert prior_weather["source_doi"]=="10.5441/001/1.703"
    if prior_weather["status"]["previous_GPS_source_replay"] is not True:
        raise ValueError("Historical Wyoming source q weather context absent")
    catalog=old_frozen()
    if catalog!=previous_manifest:
        raise ValueError("Previous original 5-study manifest has been edited")
    entries=[z for z in catalog["source_manifest"] if z["taxon"]==c["target_native"]]
    if ([z["source_key"] for z in entries]!=c["three_source_keys"] or
        [z["doi"] for z in entries]!=c["three_source_dois"] or len(entries)!=3 or
        any(z["original_file"]["bytes"]>c["max_original_GPS_bytes"] for z in entries)):
        raise ValueError("Previously selected 3 original source studies changed")
    previous_month_map={z["source_key"]:z for z in previous_month["all_five_frozen_original_source_native_analyses_including_HOLD"]}
    previous_match_map={z["source_key"]:z for z in previous_match["all_five_original_study_source_native_matched_month_results_including_HOLD"]}
    for z in entries:
        if (previous_month_map[z["source_key"]]["source_gps_status"]!="MONTH_NATIVE_ESTIMABLE"
            or previous_match_map[z["source_key"]]["source_status"]!="MATCHED_SUPPORT_PAIRED_SOURCE_ESTIMABLE"):
            raise ValueError("Previous source native q/month cohort eligibility changed")
    return c,entries


def source_days(tracks):
    admitted=source_admitted_daily(tracks)
    per_cell=defaultdict(set)
    for _tag,_year,rows in admitted:
        for d,L,q,cell,ym in rows:
            if d.year!=ym[0] or d.month!=ym[1]:
                raise ValueError("Original source day/year-month discordance")
            per_cell[cell].add(d)
    return admitted,per_cell


def physical_fingerprint(weather:dict):
    """Single SHA on all new source daily weather arrays; never export grids/dates."""
    h=hashlib.sha256()
    for cell,table in sorted(weather.items()):
        h.update(repr(cell).encode())
        for day,values in sorted(table.items()):
            h.update(day.isoformat().encode())
            if values is None:h.update(b"NULL")
            else:h.update("|".join(f"{v:.6f}" for v in values).encode())
    return h.hexdigest()


def source(index:int):
    c,entries=frozen()
    if type(index) is not int or not 0<=index<3:raise ValueError("Only exactly three frozen original mule deer source IDs")
    s=entries[index]
    md=s["original_file"]
    result={
        "version":SOURCE_VERSION,
        "source_index":index,
        "source_key":s["source_key"],
        "original_study_DOI":s["doi"],
        "original_study_id":str(s["study_id"]),
        "original_native_taxon":s["taxon"],
        "registered_original_native_taxa":s["registered_taxa"],
        "original_source_file_MD5_verified":False,
        "original_source_native_q_and_daylength_replayed_before_weather":False,
        "original_GPS_source_event_rows_private_only":0,
        "source_status":"SOURCE_HOLD",
        "source_hold_stage":None,
        "source_hold_type":None,
        "source_weather_context_version":"NEW_20261116_ERA5_DAILY_NOT_OLD_RAW_SHA_REPLAY",
        "source_ERA5_cell_count":0,
        "source_ERA5_grid_HOLD_count":0,
        "new_weather_physical_SHA256_not_coordinates":None,
        "source_weather_grid_results_differ_from_historical_not_compared":True,
        "original_q_eligible_tag_years":0,
        "all_source_native_weather_adjusted_effects":None,
        "original_48_sunlight_heading_external_NOT_REPLICATED":True,
        "original_12_ERA5_raw_SHA_HOLD_unchanged":True,
        "no_original_tag_ID_exact_dates_coordinates_GPS_weather_events_exported":True,
    }
    ext=(".zip" if md["filename"].lower().endswith(".zip") else
         ".csv.gz" if md["filename"].lower().endswith(".csv.gz") else ".csv")
    stage="PREVIOUS_PROVIDER_ORIGINAL_UUID_SIZE_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_3mule_env_private_") as directory:
            path=Path(directory)/("original"+ext)
            ff=[f for f in original_files(s["package_uuid"])
                if f["bitstream_uuid"]==md["uuid"]]
            if (len(ff)!=1 or ff[0]["bytes"]!=md["bytes"] or
                ff[0].get("checksum",{}).get("checkSumAlgorithm","").upper()!="MD5" or
                ff[0]["checksum"]["value"].lower()!=md["md5"].lower()):
                raise ValueError("Original provider file identity/MD5 mismatch")
            if not checksum_and_download(ff[0],path)["provider_checksum_verified"]:
                raise ValueError("Original provider file bytes MD5 verification failed")
            result["original_source_file_MD5_verified"]=True
            stage="SOURCE_ORIGINAL_GPS_AND_Q_BASELINE_REPLAY"
            if s["source_key"]=="mule_old_utah":
                tracks,counts,unknown=extract_time_only_original_source(path,s["registered_taxa"])
                nrows=sum(z["rows"] for z in counts.values())+unknown
                nonzero={t:counts[t]["rows"]>0 for t in s["registered_taxa"]}
            else:
                tracks,coverage,report=load_source_taxa(path,md["filename"],s["registered_taxa"])
                nrows=report["original_total_rows"]
                nonzero={t:coverage[t]["original_taxon_rows"]>0 for t in s["registered_taxa"]}
            result["original_GPS_source_event_rows_private_only"]=nrows
            if not all(nonzero.values()):raise ValueError("Original source registered native taxonomy missing")
            replay_original_baseline(s,tracks,catalog_tolerance(c))
            result["original_source_native_q_and_daylength_replayed_before_weather"]=True
            stage="ORIGINAL_GPS_DAYPART_AND_ERA5_CONTEXT_LINK"
            admitted,cells=source_days(tracks[s["taxon"]])
            result["original_q_eligible_tag_years"]=len(admitted)
            result["source_ERA5_cell_count"]=len(cells)
            if not 0<len(cells)<=c["max_new_source_weather_cells"]:
                raise ValueError("Original GPS q-dates require >160/source or zero weather cells; no geographic pruning")
            weather={};errors=0;excessive=0;total_daily=0
            for cell,dates in sorted(cells.items()):
                if (max(dates)-min(dates)).days+30>c["maximum_source_weather_request_days"]:
                    excessive+=1;errors+=1;continue
                try:
                    table,receipt=daily_source_era5(cell,dates)
                    weather[cell]=table
                    total_daily+=receipt["total_daily_rows"]
                except Exception:
                    errors+=1
            result["source_ERA5_grid_HOLD_count"]=errors
            result["source_ERA5_grid_REQUEST_DAYS_excessive"]=excessive
            result["source_ERA5_daily_grid_records_total_nonindependent"]=total_daily
            result["new_weather_physical_SHA256_not_coordinates"]=physical_fingerprint(weather) if weather else None
            # If any provider climate context failed, keep source result separately
            # labeled partial; NEVER call this complete all-source weather evidence.
            all_=adjusted_daily_effects(admitted,weather)
            all_["source_weather_grid_count_before_request"]=len(cells)
            all_["all_weather_grid_inputs_retrieved"]=errors==0
            result["all_source_native_weather_adjusted_effects"]=all_
            result["source_status"]=(
                "ORIGINAL_SOURCE_ERA5_FULL_GRID_COVERAGE_EVALUATED" if errors==0 else
                "SOURCE_PARTIAL_ERA5_GRID_HOLD_DESCRIPTIVE_ONLY"
            )
            tracks.clear()
    except Exception as e:
        result["source_status"]="SOURCE_HOLD"
        result["source_hold_stage"]=stage
        result["source_hold_type"]=type(e).__name__
        result["all_source_native_weather_adjusted_effects"]=None
    return result


def catalog_tolerance(contract):return 1e-10


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--index",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    r=source(a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    coeff=r.get("all_source_native_weather_adjusted_effects") or {}
    daylight=coeff.get("daylight_given_tmax_and_water") or {}
    print("RSDS_THREE_MULE_CLIMATE_SOURCE",json.dumps({
        "index":a.index,"status":r["source_status"],"full_q_replayed":r["original_source_native_q_and_daylength_replayed_before_weather"],
        "original_cells":r["source_ERA5_cell_count"],
        "weather_holds":r["source_ERA5_grid_HOLD_count"],
        "climate_daylength_joint_mean":daylight.get("mean"),
        "climate_daylength_joint_95":daylight.get("animal_bootstrap_95"),
        "n_original_tagged_animals":daylight.get("admitted_animals"),
        "hold_stage":r["source_hold_stage"],"hold_type":r["source_hold_type"],
    },sort_keys=True),flush=True)


if __name__=="__main__":main()
