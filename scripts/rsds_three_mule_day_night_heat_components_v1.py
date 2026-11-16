"""Original GPS DAY and NIGHT speed components behind three-study temperature~q.

This is post-outcome descriptive decomposition, NOT thermoregulation causation.
The key scientific protection is to reproduce the EXACT previously admitted
original tag/year/solar date q values and the exact previous q~Tmax month model,
then estimate NIGHT and DAY on that SAME original source+weather opportunity.
q is algebraically log1p(median NIGHT speed) minus log1p(median DAY speed).
No original private animal/tag, coordinates, exact dates or response details
leave the GitHub runner. The previously failed 48-source solar-direction
external replication and historical 12 ERA5 source raw-SHA holds are unchanged.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import random
from statistics import mean, median, pstdev
import tempfile

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned, git_blob
from scripts import rsds_three_mule_within_month_joint_weather_v1 as old
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as origin
from scripts.movebank_native_ecology_20261019 import (
    original_files, checksum_and_download, daylength,
)
from scripts.movebank_solar_movement_20261019 import valid_interval
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.rsds_five_source_native_tag_month_conditioned_v1 import replay_original_baseline
from scripts.movebank_wyoming_daily_climate_20261023 import (
    daily_source_era5, adjusted_daily_effects, climate_at_gps_day,
)
from scripts.movebank_utah_calendar_condition_20261022 import residuals

CONTRACT=HERE/"research/RSDS_THREE_MULE_DAY_NIGHT_HEAT_COMPONENTS_CONTRACT_20261116.json"
VERSION="RSDS_THREE_MULE_SOURCE_NATIVE_DAY_NIGHT_HEAT_COMPONENTS_V1"
SOURCE_VERSION="RSDS_THREE_MULE_ORIGINAL_SOURCE_DAY_NIGHT_THERMAL_COMPONENT_FIRST_V1"
KEYS=(
    "raw_day_tmax","raw_night_tmax","raw_q_tmax",
    "month_day_tmax","month_night_tmax","month_q_tmax",
    "delta_day_month_minus_raw","delta_night_month_minus_raw",
    "delta_q_month_minus_raw",
)
NBOOT=2000
TOL=1e-10


def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    expected={
        "version":VERSION,
        "original_old_three_source_month_climate_contract_blob":
            "4e5d31d06b157db329b50018dc87754f0bd7f311",
        "original_old_three_source_month_climate_first_result_blob":
            "3ea7481049cc01e6601e6f2d76e482124e503a2a",
        "original_old_three_source_month_climate_script_blob":
            "42cfb74f8a6542b46c4e3158908fa207d884ef3b",
        "original_daily_source_admission_script_blob":
            "d98b7c8aed8dd8cc603dd6516e44efc671008556",
        "original_original_GPS_step_interval_script_blob":
            "7e02edd1e806a9adcdbc60b6a3943bca721d4145",
        "target_species":"Odocoileus hemionus",
        "original_source_IDs":[
            "10.5441/001/1.711","10.5441/001/1.703","10.5441/001/1.708"],
        "original_first_month_source_individual_count_by_source":[1596,109,129],
        "original_month_Tmax_joint_nominal_direction_by_source":
            ["POSITIVE","POSITIVE","POSITIVE"],
        "original_month_Tmax_joint_means_by_source":[
            0.0010870861267833559,0.0005708082421229938,
            0.002045894391447765],
        "source_independence_verified":False,
        "thermal_causality_verified":False,
        "study_species_selection_made_after_prior_heat_findings":True,
        "same_source_weather_sha_required":True,
        "no_raw_GPS_trajectory_source_tag_ID_exact_dates_private_ERA5_cell_or_individual_component_results_exported":True,
        "previous_48_irradiance_heading_external_signed_result":"NOT_REPLICATED",
        "previous_12_ERA5_raw_SHA_source_status":"HOLD",
    }
    if any(c.get(k)!=value for k,value in expected.items()):
        raise ValueError("Source DAY/NIGHT thermal q component science route changed")
    previous_contract=pinned(
        c["original_old_three_source_month_climate_contract_path"],
        c["original_old_three_source_month_climate_contract_blob"])
    previous=pinned(c["original_old_three_source_month_climate_first_result_path"],
                    c["original_old_three_source_month_climate_first_result_blob"])
    for route in (
        "original_old_three_source_month_climate_script",
        "original_daily_source_admission_script",
        "original_original_GPS_step_interval_script",
    ):
        p=c[route+"_path"]
        if git_blob((HERE/p).read_bytes())!=c[route+"_blob"]:
            raise ValueError("Historical source q/GPS exact method code drift: "+route)
    p_contract,sources,old_records=old.freeze()
    if previous_contract!=p_contract or len(sources)!=3 or len(old_records)!=3:
        raise ValueError("Previously frozen original GPS 3-source scientific identity changed")
    if (previous["original_source_study_DOI_denominator"]!=3 or
        previous["source_status_counts"]!={"SOURCE_MONTH_JOINT_ESTIMABLE":3} or
        previous["nominal_within_month_Tmax_positive_source_count"]!=3):
        raise ValueError("Previous original q~Tmax month source result changed")
    for i,(s,prior) in enumerate(zip(sources,old_records)):
        current=previous["all_original_three_native_source_studies_including_HOLD"][i]
        if (s["doi"]!=c["original_source_IDs"][i] or
            prior["original_study_DOI"]!=s["doi"] or
            current["source_original_DOI"]!=s["doi"] or
            current["source_status"]!="SOURCE_MONTH_JOINT_ESTIMABLE" or
            current["new_month_same_tag_model"]["eligible_original_tagged_animals_joint_all_three"]!=c["original_first_month_source_individual_count_by_source"][i] or
            current["source_month_Tmax_nominal_CI"]!="POSITIVE" or
            abs(current["new_month_same_tag_model"]["month_partial_joint_estimands"]["month_Tmax_given_L_water"]["estimate"]-
                c["original_month_Tmax_joint_means_by_source"][i])>1e-12):
            raise ValueError("Original 3-source thermal respondent cohort changed")
    return c,sources,old_records,previous["all_original_three_native_source_studies_including_HOLD"]


def original_source_components(animals):
    """Reproduce original source_admitted_daily including non-response day IDs.

    Separate speed-components use the SAME original valid_interval, solar dates,
    GPS hour selection, speed bounds, per-date daypart >=3 and tag-year >=10.
    Return private original tag-year/date components only; NEVER serialize.
    """
    result=[]
    for animal, hour_events in animals.items():
        ordered=sorted(hour_events.values())
        daily=defaultdict(lambda:{"day":[],"night":[],"lat":[],"lon":[]})
        for a,b in zip(ordered,ordered[1:]):
            v=valid_interval(a[0],a[1],a[2],b[0],b[1],b[2])
            if v is None:continue
            d,part,velocity,hours,lat=v
            day=daily[d]
            day[part].append(velocity)
            day["lat"].append(lat)
            day["lon"].append((a[2]+b[2])/2)
        years=defaultdict(list)
        for d,record in daily.items():
            if len(record["day"])<3 or len(record["night"])<3:continue
            lat=mean(record["lat"])
            lon=mean(record["lon"])
            yday=math.log1p(median(record["day"]))
            ynight=math.log1p(median(record["night"]))
            q=ynight-yday
            years[d.year].append((
                d,daylength(lat,d),q,(math.floor(lat),math.floor(lon)),
                (d.year,d.month),yday,ynight,
            ))
        for year,days in years.items():
            if len(days)>=10:
                result.append((animal,year,days))
    return result


def verify_original_daily_q_reconstruction(original, components):
    """Enforce exact original source animal/year/solar-date cohort and daily q."""
    a={(tag,year):rows for tag,year,rows in original}
    b={(tag,year):rows for tag,year,rows in components}
    if len(a)!=len(original) or len(b)!=len(components) or set(a)!=set(b):
        raise ValueError("Original source q-eligible tagged animal-year identity changed")
    days=0
    max_abs=0.
    for key,prior in a.items():
        recon=b[key]
        old={z[0]:z for z in prior}
        fresh={z[0]:z for z in recon}
        if len(old)!=len(prior) or len(fresh)!=len(recon) or set(old)!=set(fresh):
            raise ValueError("Original source q-qualified GPS solar-date mismatch")
        for d,o in old.items():
            r=fresh[d]
            if o[3]!=r[3] or o[4]!=r[4]:
                raise ValueError("Original source q GPS weather-grid/year-month changed")
            if not math.isclose(o[1],r[1],rel_tol=0,abs_tol=TOL):
                raise ValueError("Original GPS local astronomical daylength changed")
            if not math.isclose(o[2],r[2],rel_tol=0,abs_tol=TOL):
                raise ValueError("Original source original day q speed contrast changed")
            if not math.isclose(r[6]-r[5],o[2],rel_tol=0,abs_tol=TOL):
                raise ValueError("Original DAY/NIGHT q speed-algebra identity changed")
            max_abs=max(max_abs,abs((r[6]-r[5])-o[2]))
            days+=1
    return {"original_q_reconstructed_eligible_tag_years":len(a),
            "original_q_reconstructed_tag_solar_dates_not_independent":days,
            "maximum_original_q_day_night_algebra_absolute_error":max_abs,
            "original_source_tag_day_grid_yearmonth_coverage_exactly_identical":True}


def heat_partial(records,response_index):
    predictors=[0,1,2]
    x=[r[0] for r in records]
    y=[r[response_index] for r in records]
    controls=[[r[j] for r in records] for j in (1,2)]
    rx=residuals(x,controls)
    ry=residuals(y,controls)
    if rx is None or ry is None or pstdev(rx)<.25:
        return None
    denom=sum(v*v for v in rx)
    return sum(a*b for a,b in zip(rx,ry))/denom if denom>1e-12 else None


def matched_year_component_effect(rows):
    """Original weather-complete entries d,tmax,water,L,q,day,night.

    Month admission on previously q/weather criteria, identical sample for
    q,D,N in BOTH unadjusted and within-tag-year-month models.
    """
    qrows=[(d,t,w,L,q) for d,t,w,L,q,day,night in rows]
    first=old.within_year_month(qrows)
    if not first["supported"]:return None
    month=defaultdict(list)
    for d,t,w,L,q,day,night in rows:
        if not all(math.isfinite(v) for v in (t,w,L,q,day,night)):
            raise ValueError("Nonfinite original source DAY/NIGHT or original climate")
        if abs((night-day)-q)>TOL:
            raise ValueError("Original date night-day q identity changed")
        month[d.month].append((float(t),float(L),float(w),
                               float(q),float(day),float(night)))
    groups=[]
    for m,records in sorted(month.items()):
        if len(records)<5 or max(z[1] for z in records)-min(z[1] for z in records)<.05:
            continue
        groups.append(records)
    if len(groups)<2 or sum(len(z) for z in groups)!=first["eligible_original_dates"]:
        raise ValueError("Original source year-month weather-eligible date pool changed")
    raw=[row for group in groups for row in group]
    conditional=[]
    for group in groups:
        average=[mean(z[i] for z in group) for i in range(6)]
        conditional.extend(tuple(row[i]-average[i] for i in range(6))
                           for row in group)
    out={}
    for prefix,records in (("raw",raw),("month",conditional)):
        day=heat_partial(records,4)
        night=heat_partial(records,5)
        q=heat_partial(records,3)
        if any(v is None or not math.isfinite(v) for v in (day,night,q)):
            raise ValueError("Original q weather rank accepted but component projection failed")
        if abs((night-day)-q)>TOL:
            raise ValueError("Original source per-year q~Tmax != night-day heat components")
        original_key=prefix+"_Tmax_given_L_water"
        if not math.isclose(q,first["estimates"][original_key],
                            rel_tol=0,abs_tol=TOL):
            raise ValueError("Previously fixed joint q source year/month heat coefficient changed")
        out[prefix+"_day_tmax"]=day
        out[prefix+"_night_tmax"]=night
        out[prefix+"_q_tmax"]=q
    out["delta_day_month_minus_raw"]=out["month_day_tmax"]-out["raw_day_tmax"]
    out["delta_night_month_minus_raw"]=out["month_night_tmax"]-out["raw_night_tmax"]
    out["delta_q_month_minus_raw"]=out["month_q_tmax"]-out["raw_q_tmax"]
    if abs((out["delta_night_month_minus_raw"]-
            out["delta_day_month_minus_raw"])-out["delta_q_month_minus_raw"])>TOL:
        raise ValueError("Original source month-minus-raw q component identity changed")
    return out


def summarize_original_tag_components(components,weather,prior,seed):
    by_tag=defaultdict(dict)
    weather_matches=0
    for tag,year,days in components:
        if year in by_tag[tag]:
            raise ValueError("Duplicate original source tagged animal-year")
        ready=[]
        for d,L,q,cell,ym,day,night in days:
            cl=climate_at_gps_day(weather.get(cell,{}),d)
            if cl is None:continue
            tm,tmax,water=cl
            ready.append((d,tmax,water,L,q,day,night))
            weather_matches+=1
        by_tag[tag][year]=ready
    eligible={}
    tagyears=0;original_dates=0
    for tag,yearmap in by_tag.items():
        estimates=[]
        for year,records in yearmap.items():
            fit=matched_year_component_effect(records)
            if fit is not None:
                estimates.append(fit)
                tagyears+=1
                # Counts are provenance/observation opportunities, not sample size
                original_dates+=len(records)
        if estimates:
            eligible[tag]={key:mean(z[key] for z in estimates) for key in KEYS}
    n=len(eligible)
    if n!=prior["new_month_same_tag_model"]["eligible_original_tagged_animals_joint_all_three"]:
        raise ValueError("Original source same q-eligible independent tag cohort changed")
    if n<5:
        return {"source_status":"SOURCE_COMPONENT_UNAVAILABLE",
                "eligible_original_tagged_animals":n,
                "original_component_tag_years":tagyears,
                "original_weather_qualified_candidate_dates":weather_matches,
                "original_component_tag_dates_not_independent":original_dates,
                "components":{key:{"status":"UNAVAILABLE","estimate":None,
                                  "tag_bootstrap95":None,"nominal_sign":"UNAVAILABLE",
                                  "original_tag_n":n} for key in KEYS},
                "all_three_endpoints_identical_tag_date_weather_support":True}
    names=sorted(eligible)
    result={key:[] for key in KEYS}
    rng=random.Random(seed)
    # Exactly the previous q-source bootstrap original-tag draw stream for
    # the shared source; never independently resample DAY/NIGHT.
    for i in range(NBOOT):
        ids=[names[rng.randrange(n)] for j in range(n)]
        for key in KEYS:
            result[key].append(mean(eligible[z][key] for z in ids))
        for p in ("raw","month"):
            j=(result[p+"_night_tmax"][-1]-result[p+"_day_tmax"][-1])
            if abs(j-result[p+"_q_tmax"][-1])>TOL:
                raise ValueError("Bootstrap draw original q != NIGHT minus DAY")
    out={}
    for key in KEYS:
        xs=sorted(result[key])
        ci=[xs[49],xs[1949]]
        out[key]={"status":"ESTIMABLE",
                  "estimate":mean(eligible[t][key] for t in names),
                  "tag_bootstrap95":ci,
                  "nominal_sign":"POSITIVE" if ci[0]>0 else
                    "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO",
                  "original_tag_n":n}
    # Verify source first month effect EXACT including original tag-bootstrap
    for key,oldkey in (("raw_q_tmax","raw_Tmax_given_L_water"),
                       ("month_q_tmax","month_Tmax_given_L_water"),
                       ("delta_q_month_minus_raw","delta_month_minus_raw_Tmax")):
        baseline=prior["new_month_same_tag_model"]["month_partial_joint_estimands"][oldkey]
        new=out[key]
        if new["original_tag_n"]!=baseline["original_tag_n"] or \
           abs(new["estimate"]-baseline["estimate"])>TOL or \
           any(abs(a-b)>TOL for a,b in zip(new["tag_bootstrap95"],baseline["tag_bootstrap95"])):
            raise ValueError("Original q~Tmax same-source first month bootstrap is not replayed")
    for prefix in ("raw","month"):
        if abs(out[prefix+"_night_tmax"]["estimate"]-
               out[prefix+"_day_tmax"]["estimate"]-
               out[prefix+"_q_tmax"]["estimate"])>TOL:
            raise ValueError("Original source mean q != NIGHT-DAY temperature slope")
    return {
        "source_status":"SOURCE_COMPONENTS_ESTIMABLE",
        "eligible_original_tagged_animals":n,
        "original_component_tag_years":tagyears,
        "original_weather_qualified_candidate_dates":weather_matches,
        "original_component_tag_dates_not_independent":original_dates,
        "original_tag_bootstrap_draws_joint":NBOOT,
        "all_three_endpoints_identical_tag_date_weather_support":True,
        "source_q_temperature_exact_mean_and_95_CI_replayed":True,
        "q_night_minus_day_algebra_exact_in_all_tag_years_and_bootstrap":True,
        "components":out,
    }


def process(index):
    c,sources,old_records,prior_records=frozen()
    if type(index) is not int or index not in (0,1,2):
        raise ValueError("Exactly three earlier source IDs, never replaced")
    source=sources[index]
    prior=prior_records[index]
    md=source["original_file"]
    out={"version":SOURCE_VERSION,"source_index":index,
         "original_study_DOI":source["doi"],"source_key":source["source_key"],
         "original_study_id":str(source["study_id"]),
         "original_registered_taxa":source["registered_taxa"],
         "status":"SOURCE_HOLD",
         "hold_stage":None,"hold_type":None,
         "source_original_GPS_MD5_bytes_verified":False,
         "original_full_native_q_daylength_replayed":False,
         "original_source_daily_q_with_day_night_components_reconstructed":False,
         "same_20261116_daily_ERA5_physical_SHA_replayed":False,
         "original_six_weather_source_model_means_replayed":False,
         "original_q_Tmax_month_first_bootstrap_replayed":False,
         "prior_three_study_ID_and_source_region_UNCHANGED":True,
         "no_original_GPS_tag_ids_positions_exact_dates_climate_grid_component_coefficients_exported":True,
         "old_48_irradiance_heading_direction_signed_test_NOT_REPLICATED":True,
         "old_12_ERA5_raw_SHA_HOLD":True,
         "original_GPS_source_event_rows_private_only":0,
         "source_weather_grid_count":0,"source_weather_grid_HOLD_count":0,
         "source_component_summary":None}
    ext=(".zip" if md["filename"].lower().endswith(".zip") else
         ".csv.gz" if md["filename"].lower().endswith(".csv.gz") else ".csv")
    stage="ORIGINAL_PROVIDER_MD5_UUID_BYTES"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_daynight_private_") as directory:
            path=Path(directory)/("original"+ext)
            ff=[z for z in original_files(source["package_uuid"])
                if z["bitstream_uuid"]==md["uuid"]]
            if (len(ff)!=1 or ff[0]["bytes"]!=md["bytes"] or
                ff[0].get("checksum",{}).get("checkSumAlgorithm","").upper()!="MD5" or
                ff[0]["checksum"]["value"].lower()!=md["md5"].lower()):
                raise ValueError("Previously frozen original animal GPS source file modified")
            if not checksum_and_download(ff[0],path)["provider_checksum_verified"]:
                raise ValueError("Original provider GPS bytes checksum not verified")
            out["source_original_GPS_MD5_bytes_verified"]=True
            stage="ORIGINAL_GPS_Q_AND_DAY_NIGHT_ALGEBRA"
            if source["source_key"]=="mule_old_utah":
                tracks,counts,nonmatch=extract_time_only_original_source(
                    path,source["registered_taxa"])
                out["original_GPS_source_event_rows_private_only"]=nonmatch+sum(
                    z["rows"] for z in counts.values())
            else:
                tracks,coverage,report=load_source_taxa(
                    path,md["filename"],source["registered_taxa"])
                out["original_GPS_source_event_rows_private_only"]=report["original_total_rows"]
            replay_original_baseline(source,tracks,1e-10)
            out["original_full_native_q_daylength_replayed"]=True
            admitted,cells=origin.source_days(tracks[source["taxon"]])
            component=original_source_components(tracks[source["taxon"]])
            quality=verify_original_daily_q_reconstruction(admitted,component)
            out["original_source_daily_q_with_day_night_components_reconstructed"]=True
            out["source_original_q_reconstruction_quality"]=quality
            stage="OLD_VERIFIED_ERA5_SOURCE_PHYSICAL_VERSION"
            if len(cells)!=old_records[index]["new_weather_source_grid_count"]:
                raise ValueError("Original source GPS q year/cell weather opportunity changed")
            out["source_weather_grid_count"]=len(cells)
            weather={};errors=0
            for cell,dates in sorted(cells.items()):
                try:
                    table,receipt=daily_source_era5(cell,dates)
                    weather[cell]=table
                except Exception:
                    errors+=1
            out["source_weather_grid_HOLD_count"]=errors
            if errors or len(weather)!=len(cells) or \
               origin.physical_fingerprint(weather)!=old_records[index]["current_provider_daily_physical_fingerprint"]:
                raise ValueError("ERA5 source missing/different from 2026-11-16 physical snapshot")
            out["same_20261116_daily_ERA5_physical_SHA_replayed"]=True
            stage="ORIGINAL_SIX_SOURCE_WEATHER_AND_MONTH_Q_REPLAY"
            original_models=adjusted_daily_effects(admitted,weather)
            for key,past in old_records[index]["source_weather_adjusted_endpoints"].items():
                now=original_models[key]
                if (now["admitted_animals"]!=past["original_source_tagged_animals"]
                    or now["admitted_animal_years"]!=past["source_original_animal_years"]
                    or ((now["mean"] is None)!=(past["mean"] is None)) or
                    (now["mean"] is not None and
                     not math.isclose(now["mean"],past["mean"],
                                      rel_tol=0,abs_tol=TOL))):
                    raise ValueError("Original six daily weather source group results drifted")
            out["original_six_weather_source_model_means_replayed"]=True
            result=summarize_original_tag_components(
                component,weather,prior,seed=20261116+1000*index)
            out["original_q_Tmax_month_first_bootstrap_replayed"]=(
                result.get("source_q_temperature_exact_mean_and_95_CI_replayed") is True)
            out["source_component_summary"]=result
            out["status"]=result["source_status"]
            weather.clear();tracks.clear()
    except Exception as error:
        out["status"]="SOURCE_HOLD"
        out["hold_stage"]=stage
        out["hold_type"]=type(error).__name__
        out["source_component_summary"]=None
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--index",required=True,type=int)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    result=process(a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    c=(result["source_component_summary"] or {}).get("components") or {}
    print("RSDS_THREE_MULE_DAY_NIGHT_HEAT_SOURCE",json.dumps({
        "i":a.index,"status":result["status"],"hold_stage":result["hold_stage"],
        "hold_type":result["hold_type"],
        "source_q_bootstrap_replayed":result["original_q_Tmax_month_first_bootstrap_replayed"],
        "tag_n":(result["source_component_summary"] or {}).get("eligible_original_tagged_animals"),
        "month_DAY_heat":c.get("month_day_tmax"),
        "month_NIGHT_heat":c.get("month_night_tmax"),
        "month_Q_heat":c.get("month_q_tmax"),
    },sort_keys=True),flush=True)


if __name__=="__main__":main()
