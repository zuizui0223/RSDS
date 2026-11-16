"""All-six native ungulate q/DAY/NIGHT Tmax within original tag×year×month.

New original GPS+fresh contemporaneous ERA5 weather only after immutable
two-source/four-plus-two-native manifest freeze. Old q and daylight original
native-tag means/n must REPLAY before obtaining any new weather data.
The source DOI is the unit, and each source's species can have different
telemetry locations and dates: this does NOT estimate a species×Tmax effect.
Raw IDs/dates/GPS/ERA5 grid and individual coefficient draws stay ephemeral.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import tempfile

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned,git_blob
from scripts import rsds_three_mule_day_night_heat_components_v1 as daynight
from scripts import rsds_three_mule_within_month_joint_weather_v1 as month
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as climate
from scripts.movebank_native_ecology_20261019 import (
    original_files,checksum_and_download
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.movebank_wyoming_daily_climate_20261023 import (
    daily_source_era5, climate_at_gps_day
)

CONTRACT=HERE/"research/RSDS_TWO_SOURCE_SIX_UNGULATE_NATIVE_MONTH_HEAT_CONTRACT_20261117.json"
VERSION="RSDS_TWO_ORIGINAL_UNGULATE_STUDIES_ALL_SIX_NATIVE_WITHIN_MONTH_TMAX_V1"
SOURCE_VERSION="RSDS_TWO_SOURCE_ALL_NATIVE_MONTH_HEAT_NIGHT_DAY_FIRST_V1"
TOL=1e-10
NATIVE_KEYS=("month_day_tmax","month_night_tmax","month_q_tmax")

def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if (c.get("version")!=VERSION or c.get("source_studies")!=2
        or c.get("original_native_taxon_source_rows")!=6
        or c.get("source_weather_max_grids")!=160
        or c.get("source_weather_request_max_days")!=3653
        or c.get("source_original_max_bytes")!=200_000_000
        or c.get("previous_48_solar_heading_signed_external")!="NOT_REPLICATED"
        or c.get("previous_12_old_ERA5_raw_SHA")!="HOLD"):
        raise ValueError("Original two source six native species contract changed")
    pinnedfiles={}
    for z in c["immutable_prior_files"]:
        if z["path"].endswith(".py"):
            if git_blob((HERE/z["path"]).read_bytes())!=z["blob"]:
                raise ValueError("Original source/native code Git blob version drift")
        else:
            pinnedfiles[z["path"]]=pinned(z["path"],z["blob"])
    ut=pinnedfiles["research/MOVEBANK_UTAH_FOUR_TAXA_SOLAR_V0_20261019.json"]["source"]
    wy=pinnedfiles["research/MOVEBANK_WYOMING_TWO_UNGULATES_SOURCE_V0_20261022.json"]["source"]
    original=[
       {
         "doi":ut["doi"],"id":c["source_0"]["study_id"],"package":ut["package_uuid"],
         "uuid":ut["original_bitstream_uuid"],"size":ut["original_bytes"],"md5":ut["original_md5"],
         "filename":ut["original_file"],"taxa":ut["catalog_taxa"],"reader":"UTAH_ZIP"
       },
       {
         "doi":wy["doi"],"id":c["source_1"]["study_id"],"package":wy["package_uuid"],
         "uuid":wy["original_bitstream_uuid"],"size":wy["original_bytes"],"md5":wy["original_md5"],
         "filename":wy["original_filename"],"taxa":wy["published_species"],"reader":"SOURCE_CSV"
       }
    ]
    for i,s in enumerate(original):
        p=c["source_"+str(i)]
        if (s["doi"]!=p["doi"] or s["id"]!=p["study_id"] or
            s["package"]!=p["package_uuid"] or s["uuid"]!=p["provider_file_uuid"] or
            s["size"]!=p["provider_original_bytes"] or
            s["md5"]!=p["provider_MD5"] or s["filename"]!=p["original_filename"] or
            s["taxa"]!=p["native_taxa"] or s["size"]>c["source_original_max_bytes"]):
            raise ValueError("Original provider source identity or registered taxa drift")
    l=pinnedfiles["research/RSDS_12_STUDY_22_NATIVE_TEMPORAL_MACRO_FIRST_RESULT_20261031.json"]
    legacy=defaultdict(dict)
    for row in l["all_original_study_native_taxon_rows_including_HOLD"]:
        if row["original_study_doi"] not in (original[0]["doi"],original[1]["doi"]):continue
        legacy[row["original_study_doi"]][row["native_taxon"]]=row
    if (len(legacy)!=2 or
        any(set(legacy[s["doi"]])!=set(s["taxa"]) for s in original) or
        sum(len(v) for v in legacy.values())!=6):
        raise ValueError("Original full source all 6 native baseline identities missing")
    old_heat=pinnedfiles["research/RSDS_UTAH_WYOMING_SIX_UNGULATE_NATIVE_HEAT_Q_REFERENCE_20261116.json"]
    if (old_heat["source_count"]!=2 or old_heat["source_native_rows"]!=6
        or len(old_heat["original_tmax_q_native_rows"])!=6):
        raise ValueError("All previously examined historical native source q records altered")
    return c,original,legacy


def source_eligible_climate_or_HOLD(native_admitted,weather):
    total=missing=0
    for rows in native_admitted.values():
        for tag,year,days in rows:
            for d,L,q,cell,ym in days:
                total+=1
                if climate_at_gps_day(weather.get(cell,{}),d) is None:
                    missing+=1
    if missing:
        raise ValueError("Not all original native weather source/q dates have complete past-30d provider ERA5")
    return total


def original_native_month_result(original_q,parts,weather,seed):
    first=month.per_source_same_sample_models(original_q,weather,seed)
    # This is the SAME old original q/day observation opportunity and model
    # seed, but a NEW provider weather snapshot, so only this NEW run is replayed.
    checked=daynight.summarize_original_tag_components(
       parts,weather,{"new_month_same_tag_model":first},seed)
    if checked["source_status"]=="SOURCE_COMPONENTS_ESTIMABLE":
        orig=first["month_partial_joint_estimands"]["month_Tmax_given_L_water"]
        night_q=checked["components"]["month_q_tmax"]
        if (orig["original_tag_n"]!=night_q["original_tag_n"] or
            abs(orig["estimate"]-night_q["estimate"])>TOL or
            any(abs(a-b)>TOL for a,b in zip(orig["tag_bootstrap95"],
                                              night_q["tag_bootstrap95"]))):
            raise ValueError("Original q temperature/source native bootstrap not reconstructed")
        for name in ("raw","month"):
            x=checked["components"]
            if abs(x[name+"_night_tmax"]["estimate"]-
                   x[name+"_day_tmax"]["estimate"]-
                   x[name+"_q_tmax"]["estimate"])>TOL:
                raise ValueError("Native source q temperature identity != night-day")
    return first,checked


def source(index):
    c,sources,legacy=frozen()
    if type(index)!=int or index not in (0,1):
        raise ValueError("Only two fixed original source study DOI IDs accepted")
    s=sources[index]
    out={
        "version":SOURCE_VERSION,"source_index":index,
        "original_study_DOI":s["doi"],"original_repository_study_id":s["id"],
        "registered_native_taxa":s["taxa"],
        "source_status":"SOURCE_HOLD","hold_stage":None,"hold_type":None,
        "source_original_GPS_MD5_verified":False,
        "all_registered_native_q_and_solar_group_first_results_replayed":False,
        "original_q_NIGHT_minus_DAY_per_tag_date_reconstructed":False,
        "new_contemporaneous_source_native_ERA5_full_grid_coverage":False,
        "source_ERA5_daily_physical_SHA256":None,
        "source_ERA5_grid_count":0,"source_ERA5_grid_HOLDs":0,
        "source_native_group_results":[],
        "original_source_mule_and_bighorn_shared_source_grid_yearmonth_count":None,
        "original_source_mule_and_bighorn_shared_grid_yearmonth_with_3_tags_each":None,
        "source_new_weather_is_not_previous_20261023_or_20261116_physical_replay":True,
        "no_original_animal_ID_tracked_coordinates_exact_date_GPS_weather_grid_or_individual_effect_exported":True,
        "previous_48_source_irradiance_heading_direction_NOT_REPLICATED":True,
        "previous_12_source_ERA5_raw_SHA_HOLD":True,
    }
    ext=".zip" if s["filename"].lower().endswith(".zip") else ".csv"
    stage="ORIGINAL_PROVIDER_SOURCE_UUID_MD5_BYTES"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_2source_6native_private_") as tmp:
            dest=Path(tmp)/("original"+ext)
            source_files=[r for r in original_files(s["package"]) if r["bitstream_uuid"]==s["uuid"]]
            if (len(source_files)!=1 or
                source_files[0]["bytes"]!=s["size"] or
                source_files[0]["checksum"]["checkSumAlgorithm"].upper()!="MD5" or
                source_files[0]["checksum"]["value"].lower()!=s["md5"]):
                raise ValueError("Original published source file UUID/size/MD5 altered")
            if not checksum_and_download(source_files[0],dest)["provider_checksum_verified"]:
                raise ValueError("Original published source bytes did not pass MD5")
            out["source_original_GPS_MD5_verified"]=True
            stage="ORIGINAL_NATIVE_SIX_Q_DAYLENGTH_SOURCE_GPS_REPLAY"
            if s["reader"]=="UTAH_ZIP":
                tracks,counts,unmatched=extract_time_only_original_source(dest,s["taxa"])
                found={tax:counts[tax]["rows"]>0 for tax in s["taxa"]}
                out["source_total_original_GPS_rows_seen_in_private_runner"]=sum(
                    z["rows"] for z in counts.values())+unmatched
            else:
                tracks,coverage,metadata=load_source_taxa(dest,s["filename"],s["taxa"])
                found={tax:coverage[tax]["original_taxon_rows"]>0 for tax in s["taxa"]}
                out["source_total_original_GPS_rows_seen_in_private_runner"]=metadata["original_total_rows"]
            if not all(found.values()):raise ValueError("One native original source species missing from original GPS")
            admitted_by_native={}
            components_by_native={}
            climate_cells=defaultdict(set)
            shared=defaultdict(lambda:defaultdict(set))
            for native in s["taxa"]:
                baseline=legacy[s["doi"]][native]
                source_effect=one_taxon(native,tracks[native])[0]
                for f,oldkey,current in (
                  ("q","day_night_speed_contrast",source_effect["day_night_speed_contrast"]),
                  ("daylength","daylength_association",source_effect["daylength_association"]),
                ):
                    previous=baseline[f]
                    if (current is None or previous.get("mean") is None or
                        previous["source_tags"]!=current["n_individuals"] or
                        not math.isclose(previous["mean"],current["mean"],rel_tol=0,abs_tol=TOL)):
                        raise ValueError("Original source-native old "+native+" "+f+" q/daylight not replayed")
                admitted,cells=climate.source_days(tracks[native])
                components=daynight.original_source_components(tracks[native])
                daynight.verify_original_daily_q_reconstruction(admitted,components)
                admitted_by_native[native]=admitted
                components_by_native[native]=components
                for cell,dates in cells.items():
                    climate_cells[cell].update(dates)
                if native in ("Odocoileus hemionus","Ovis canadensis"):
                    for tag,year,days in admitted:
                        for d,L,q,cell,ym in days:
                            shared[native][(cell,ym)].add(tag)
            out["all_registered_native_q_and_solar_group_first_results_replayed"]=True
            out["original_q_NIGHT_minus_DAY_per_tag_date_reconstructed"]=True
            mule=shared["Odocoileus hemionus"]
            sheep=shared["Ovis canadensis"]
            common=set(mule)&set(sheep)
            out["original_source_mule_and_bighorn_shared_source_grid_yearmonth_count"]=len(common)
            out["original_source_mule_and_bighorn_shared_grid_yearmonth_with_3_tags_each"]=sum(
                len(mule[k])>=3 and len(sheep[k])>=3 for k in common)
            stage="ALL_NATIVE_ORIGINAL_Q_CLIMATE_FULL_GRID_ERA5"
            if not 0<len(climate_cells)<=c["source_weather_max_grids"]:
                raise ValueError("All native species jointly require too many/zero originally GPS-observed 1deg cells")
            weather={};errors=0
            out["source_ERA5_grid_count"]=len(climate_cells)
            for cell,dates in sorted(climate_cells.items()):
                if (max(dates)-min(dates)).days+30>c["source_weather_request_max_days"]:
                    errors+=1;continue
                try:
                    t,receipt=daily_source_era5(cell,dates)
                    weather[cell]=t
                except Exception:
                    errors+=1
            out["source_ERA5_grid_HOLDs"]=errors
            if errors or len(weather)!=len(climate_cells):
                raise ValueError("Some originally eligible native source q days have missing ERA5 grids")
            source_eligible_climate_or_HOLD(admitted_by_native,weather)
            out["new_contemporaneous_source_native_ERA5_full_grid_coverage"]=True
            out["source_ERA5_daily_physical_SHA256"]=climate.physical_fingerprint(weather)
            stage="SIX_NATIVE_ORIGINAL_SAME_YEAR_MONTH_THERMAL_SOURCE_RESULTS"
            source_results=[]
            for i,taxon in enumerate(s["taxa"]):
                old_source=legacy[s["doi"]][taxon]
                qframe,components=original_native_month_result(
                    admitted_by_native[taxon],components_by_native[taxon],
                    weather,seed=20261117+index*1000+i*10)
                eff=(components.get("components") or {})
                q_heat=qframe["month_partial_joint_estimands"]["month_Tmax_given_L_water"]
                if components["source_status"]=="SOURCE_COMPONENTS_ESTIMABLE":
                    status="SOURCE_NATIVE_MONTH_COMPONENTS_ESTIMABLE"
                else:
                    status="SOURCE_NATIVE_MONTH_INSUFFICIENT_TAGS"
                source_results.append({
                    "original_native_taxon":taxon,
                    "source_native_status":status,
                    "first_full_original_q":old_source["q"]["mean"],
                    "first_full_original_q_tagged_animals":old_source["q"]["source_tags"],
                    "first_original_daylength_slope":old_source["daylength"]["mean"],
                    "first_original_daylength_tagged_animals":old_source["daylength"]["source_tags"],
                    "month_joint_original_tag_n":qframe["eligible_original_tagged_animals_joint_all_three"],
                    "month_native_q_Tmax_given_daylength_antecedent_water":q_heat,
                    "month_native_DAY_Tmax_given_daylength_antecedent_water":eff.get("month_day_tmax"),
                    "month_native_NIGHT_Tmax_given_daylength_antecedent_water":eff.get("month_night_tmax"),
                    "original_native_q_NIGHT_minus_DAY_model_and_bootstrap_algebra_verified":
                        components.get("q_night_minus_day_algebra_exact_in_all_tag_years_and_bootstrap",False),
                    "first_source_native_q_Tmax_new_ERA5_tag_bootstrap_replayed":
                        components.get("source_q_temperature_exact_mean_and_95_CI_replayed",False),
                })
            out["source_native_group_results"]=source_results
            out["source_status"]="SOURCE_ALL_NATIVE_ERA5_MONTH_MODELS_EVALUATED"
            weather.clear();tracks.clear()
    except Exception as e:
        out["source_status"]="SOURCE_HOLD"
        out["hold_stage"]=stage
        out["hold_type"]=type(e).__name__
        out["source_native_group_results"]=[]
        # source-only counts before a HOLD cannot be promoted to result
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--index",required=True,type=int)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    x=source(a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("RSDS_TWO_SOURCE_ALL_NATIVE_MONTH_TEMP",json.dumps({
        "i":x["source_index"],"status":x["source_status"],
        "source_weather_grids":x["source_ERA5_grid_count"],
        "source_weather_holds":x["source_ERA5_grid_HOLDs"],
        "source_registered_native_taxa":x["registered_native_taxa"],
        "native_model_outcomes":[{
            "name":r["original_native_taxon"],
            "status":r["source_native_status"],
            "q_month_heat":r["month_native_q_Tmax_given_daylength_antecedent_water"],
            "DAY":r["month_native_DAY_Tmax_given_daylength_antecedent_water"],
            "NIGHT":r["month_native_NIGHT_Tmax_given_daylength_antecedent_water"],
        } for r in x["source_native_group_results"]],
        "shared_mule_bighorn_cellmonth":x["original_source_mule_and_bighorn_shared_source_grid_yearmonth_count"],
        "hold_stage":x["hold_stage"],"hold_type":x["hold_type"],
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
