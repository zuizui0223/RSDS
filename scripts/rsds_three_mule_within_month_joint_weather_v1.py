"""Within-animal/year/month heat, daylight, water GPS q associations in 3 studies.

Post-outcome observational sensitivity. Source native originally tracked GPS, no
raw outputs. The prior 2026-11-16 verified daily ERA5 source physical fingerprint
must replay *exactly* before a new month partial slope is calculated. Original
provider GPS UUID/bytes/MD5, original q and six old daily weather source means
are also replayed. No alternative weather grid, model or source is admitted.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import random
from statistics import mean,pstdev
import tempfile

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned,git_blob
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as orig
from scripts.movebank_native_ecology_20261019 import original_files,checksum_and_download
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.rsds_five_source_native_tag_month_conditioned_v1 import replay_original_baseline
from scripts.movebank_wyoming_daily_climate_20261023 import (
    daily_source_era5,climate_at_gps_day,adjusted_daily_effects
)
from scripts.movebank_utah_calendar_condition_20261022 import residuals

CONTRACT=HERE/"research/RSDS_THREE_MULE_WITHIN_MONTH_HEAT_CONTRACT_20261116.json"
VERSION="RSDS_THREE_MULE_WITHIN_TAG_YEAR_MONTH_JOINT_DAILY_CLIMATE_V1"
SOURCE_VERSION="RSDS_THREE_MULE_FIRST_ORIGINAL_Q_WEATHER_MATCHED_MONTH_V1"
TARGETS=("Tmax_given_L_water","L_given_Tmax_water","water_given_Tmax_L")
RAW_KEYS=("raw_Tmax_given_L_water","raw_L_given_Tmax_water","raw_water_given_Tmax_L")
MONTH_KEYS=("month_Tmax_given_L_water","month_L_given_Tmax_water","month_water_given_Tmax_L")
DELTA_KEYS=("delta_month_minus_raw_Tmax","delta_month_minus_raw_L","delta_month_minus_raw_water")
CLIMATE_KEYS=("tmax_given_daylight_and_water","daylight_given_tmax_and_water","water_given_daylight_and_tmax")
LIMITS=(.25,.05,.05)
NBOOT=2000

def freeze():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    expected={
        "version":VERSION,
        "original_daily_weather_source_contract_blob":"6a102152df8d11b27007c3fbecf45f48b20d2c4a",
        "original_first_three_source_result_blob":"25db513c9e6804a282bcc8110565e033273f03b6",
        "original_daily_source_script_blob":"c1afa02d63129a256f2a2307dab1093ba59ac2d8",
        "original_weather_group_model_script_blob":"d9d886cb469fe8f69cb61ec7df7f4d00b03b08c0",
        "original_within_month_q_contract_blob":"7a25beccbb097b266d171edfaec8945e7f7d524f",
        "prior_three_source_complete_grid_count":3,
        "prior_positive_Tmax_nominal_source_count":3,
        "prior_daylength_joint_nominal_signs":["NEGATIVE","POSITIVE","NEGATIVE"],
        "source_IDs":["10.5441/001/1.711","10.5441/001/1.703","10.5441/001/1.708"],
        "min_tag_month_q_dates":5,
        "min_tag_year_qualifying_months":2,
        "min_tag_year_qualified_days":15,
        "min_tag_month_photoperiod_span_h":.05,
        "test_window_Daylength_Sxx":.05,
        "min_residual_temperature_SD_C":.25,
        "min_residual_daylength_SD_h":.05,
        "min_residual_water_SD_per100mm":.05,
        "bootstrap_original_tag_resamples":NBOOT,
        "eligibility_min_original_tags_per_source":5,
        "source_independence_not_proven":True,
        "species_selection_after_prior_positive_temperature_findings":True,
        "original_tag_trajectories_IDs_coordinate_dates_and_individual_slopes_not_exported":True,
        "historic_48_source_solar_heading_direction_signed_external_result":"NOT_REPLICATED",
        "historic_12_ERA5_raw_SHA_source_status":"HOLD",
    }
    if any(c.get(k)!=v for k,v in expected.items()):
        raise ValueError("Previously frozen after-outcome three-mule month-weather contract changed")
    d=pinned(c["original_first_three_source_result_path"],
             c["original_first_three_source_result_blob"])
    source_contract=pinned(c["original_daily_weather_source_contract_path"],
                           c["original_daily_weather_source_contract_blob"])
    pinned(c["original_within_month_q_contract_path"],
           c["original_within_month_q_contract_blob"])
    for prefix in ("original_daily_source_script","original_weather_group_model_script"):
        if git_blob((HERE/c[prefix+"_path"]).read_bytes())!=c[prefix+"_blob"]:
            raise ValueError("Previously frozen GPS/weather SOURCE code drift: "+prefix)
    if (d.get("source_study_DOI_denominator")!=3 or
        d.get("original_source_status_counts")!={"ORIGINAL_SOURCE_ERA5_FULL_GRID_COVERAGE_EVALUATED":3} or
        d.get("full_weather_coverage_source_studies")!=3 or
        len(d["all_three_original_source_rows_including_HOLD"])!=3 or
        d["primary_daylength_joint_adjusted_nominal_CI_signs_by_fully_covered_source"]!={"mule_old_utah":"NEGATIVE","mule_old_wyoming":"POSITIVE","mule_new_sublette":"NEGATIVE"} or
        d["Tmax_given_daylength_water_nominal_CI_signs_by_fully_covered_source"]!={"mule_old_utah":"POSITIVE","mule_old_wyoming":"POSITIVE","mule_new_sublette":"POSITIVE"}):
        raise ValueError("Three original 2026 source weather results or signs changed")
    original_contract,studies=orig.frozen()
    if original_contract!=source_contract or len(studies)!=3:
        raise ValueError("Original day-level source contract/studies mutated")
    first=d["all_three_original_source_rows_including_HOLD"]
    for a,b in zip(studies,first):
        if (a["source_key"]!=b["original_source_key"] or a["doi"]!=b["original_study_DOI"] or
            b["native_GPS_source_MD5_replayed"] is not True or
            b["original_native_q_first_result_replayed"] is not True or
            b["new_weather_source_grid_HOLD_count"]!=0 or
            len(b["current_provider_daily_physical_fingerprint"])!=64):
            raise ValueError("Missing verified previous original source+weather fingerprint")
    return c,studies,first


def within_year_month(rows):
    """Rows: original (tmax, water/100 mm, original local daylight h, GPS q).

    A/B calculations use exactly the same date/weather source rows and same
    qualifying years per tag. Month centering de-confounds between months but
    cannot identify a causal sunlight/temperature treatment.
    """
    months=defaultdict(list)
    for day,tmax,water,L,q in rows:
        if not all(math.isfinite(x) for x in (tmax,water,L,q)):
            raise ValueError("Nonfinite original source qualified q/weather row")
        months[day.month].append((float(tmax),float(L),float(water),float(q)))
    qualified=[]
    for month,records in sorted(months.items()):
        if (len(records)<5 or
            max(v[1] for v in records)-min(v[1] for v in records)<.05):
            continue
        qualified.append(records)
    total=sum(len(v) for v in qualified)
    if len(qualified)<2 or total<15:
        return {"supported":False,"eligible_month_blocks":len(qualified),
                "eligible_original_dates":total,"estimates":None}
    sxx=sum(sum((r[1]-mean(k[1] for k in block))**2 for r in block)
            for block in qualified)
    if sxx<.05:
        return {"supported":False,"eligible_month_blocks":len(qualified),
                "eligible_original_dates":total,"estimates":None}
    raw=[v for block in qualified for v in block]
    demeaned=[]
    for block in qualified:
        mu=[mean(v[i] for v in block) for i in range(4)]
        demeaned.extend(tuple(v[i]-mu[i] for i in range(4)) for v in block)
    if len(raw)!=len(demeaned):
        raise ValueError("Original month weather conditioning changed date support")
    # The same FULL three-predictor joint rank/variance restrictions in both
    # models means the paired result compares the identical source tag/year
    # date opportunity, never endpoints estimated on different tag samples.
    def partial(records,i):
        other=[j for j in range(3) if j!=i]
        rx=residuals([x[i] for x in records],[[x[j] for x in records] for j in other])
        ry=residuals([x[3] for x in records],[[x[j] for x in records] for j in other])
        if rx is None or ry is None or pstdev(rx)<LIMITS[i]:return None
        denom=sum(z*z for z in rx)
        return sum(x*y for x,y in zip(rx,ry))/denom if denom>1e-12 else None
    rawcoef=[partial(raw,i) for i in range(3)]
    within=[partial(demeaned,i) for i in range(3)]
    if any(v is None or not math.isfinite(v) for v in (*rawcoef,*within)):
        return {"supported":False,"eligible_month_blocks":len(qualified),
                "eligible_original_dates":total,"estimates":None}
    return {"supported":True,"eligible_month_blocks":len(qualified),
            "eligible_original_dates":total,
            "estimates":{
                **{k:float(z) for k,z in zip(RAW_KEYS,rawcoef)},
                **{k:float(z) for k,z in zip(MONTH_KEYS,within)},
                **{k:float(m-r) for k,m,r in zip(DELTA_KEYS,within,rawcoef)},
            }}


def per_source_same_sample_models(admitted,weather,seed):
    source_tags=defaultdict(lambda:defaultdict(list))
    for tag,year,days in admitted:
        key=(tag,year)
        if key in source_tags[tag] or year in source_tags[tag]:
            raise ValueError("Repeated original tag-year source calendar admission")
        observed=[]
        for d,L,q,cell,ym in days:
            climate=climate_at_gps_day(weather.get(cell,{}),d)
            if climate is None:continue
            tm,tmax,water=climate
            observed.append((d,tmax,water,L,q))
        if observed:source_tags[tag][year]=observed
    keys=RAW_KEYS+MONTH_KEYS+DELTA_KEYS
    eligible={}
    qualified_years=0
    qualified_days=0
    qualified_months=0
    for tag,years in source_tags.items():
        vals=[]
        for year,rows in years.items():
            r=within_year_month(rows)
            if r["supported"]:
                qualified_years+=1;qualified_days+=r["eligible_original_dates"]
                qualified_months+=r["eligible_month_blocks"]
                vals.append(r["estimates"])
        if vals:
            eligible[tag]={k:mean(v[k] for v in vals) for k in keys}
    tag_n=len(eligible)
    result={
        "original_tag_qualified_month_joint_years":qualified_years,
        "original_joint_date_count_not_independent":qualified_days,
        "original_joint_tag_month_blocks":qualified_months,
        "eligible_original_tagged_animals_joint_all_three":tag_n,
        "all_nine_estimands_share_identical_original_tag_year_month_date_support":True,
        "month_partial_joint_estimands":{},
        "bootstrap_tags_resampled_jointly_all_nine":True,
    }
    if tag_n<5:
        for key in keys:
            result["month_partial_joint_estimands"][key]={
                "status":"INSUFFICIENT_ORIGINAL_TAGS",
                "original_tag_n":tag_n,
                "estimate":None,
                "tag_bootstrap95":None,
                "nominal_CI_sign":"UNAVAILABLE"
            }
        return result
    names=sorted(eligible)
    rng=random.Random(seed)
    n=tag_n
    draws=[[names[rng.randrange(n)] for _ in range(n)] for z in range(NBOOT)]
    for k in keys:
        v=[eligible[t][k] for t in names]
        sample=sorted(mean(eligible[tag][k] for tag in ids) for ids in draws)
        ci=[sample[49],sample[1949]]
        result["month_partial_joint_estimands"][k]={
            "status":"ESTIMABLE",
            "original_tag_n":tag_n,
            "estimate":mean(v),
            "tag_bootstrap95":ci,
            "nominal_CI_sign":"POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
        }
    for a,b,d in zip(RAW_KEYS,MONTH_KEYS,DELTA_KEYS):
        lhs=(result["month_partial_joint_estimands"][b]["estimate"]-
             result["month_partial_joint_estimands"][a]["estimate"])
        if not math.isclose(lhs,result["month_partial_joint_estimands"][d]["estimate"],
                            abs_tol=1e-12,rel_tol=0):
            raise ValueError("Matched-tag/date month delta source mean failed identity")
    return result


def run(index):
    c,source,prior=freeze()
    if type(index)!=int or index not in range(3):
        raise ValueError("Only three original source-native fixed study slots")
    s=source[index];p=prior[index];md=s["original_file"]
    output={
        "version":SOURCE_VERSION,"source_index":index,
        "source_key":s["source_key"],"source_original_DOI":s["doi"],
        "source_original_study_id":str(s["study_id"]),
        "original_source_registered_native_taxa":s["registered_taxa"],
        "status":"SOURCE_HOLD",
        "hold_stage":None,"hold_type":None,
        "provider_original_MD5_verified":False,
        "original_native_q_and_daylength_mean_n_replayed_before_ERA5":False,
        "original_weather_version_exact_20261116_physical_sha_replayed":False,
        "prior_six_climate_model_mean_n_replayed":False,
        "original_GPS_rows_private_only":0,
        "weather_grid_count":0,
        "weather_grid_HOLD_count":0,
        "source_month_joint":None,
        "original_48_solar_heading_signed_external_not_replicated":True,
        "original_12_ERA5_old_raw_sha_hold":True,
        "GPS_animal_ID_coordinates_dates_weather_grid_key_and_individual_slopes_not_exported":True,
    }
    ext=(".zip" if md["filename"].lower().endswith(".zip") else
         ".csv.gz" if md["filename"].lower().endswith(".csv.gz") else ".csv")
    step="ORIGINAL_SOURCE_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_3mule_wm_private_") as tmp:
            dest=Path(tmp)/("original"+ext)
            files=[x for x in original_files(s["package_uuid"]) if x["bitstream_uuid"]==md["uuid"]]
            if (len(files)!=1 or files[0]["bytes"]!=md["bytes"] or
                files[0].get("checksum",{}).get("checkSumAlgorithm","").upper()!="MD5" or
                files[0]["checksum"]["value"].lower()!=md["md5"].lower()):
                raise ValueError("Original GPS bitstream MD5/size changed")
            if not checksum_and_download(files[0],dest)["provider_checksum_verified"]:
                raise ValueError("Original provider GPS MD5 download not verified")
            output["provider_original_MD5_verified"]=True
            step="ORIGINAL_GPS_NATIVE_Q_AND_DAYLIGHT_REPLAY"
            if s["source_key"]=="mule_old_utah":
                tracks,counts,nonmatch=extract_time_only_original_source(dest,s["registered_taxa"])
                output["original_GPS_rows_private_only"]=nonmatch+sum(x["rows"] for x in counts.values())
            else:
                tracks,coverage,info=load_source_taxa(dest,md["filename"],s["registered_taxa"])
                output["original_GPS_rows_private_only"]=info["original_total_rows"]
            replay_original_baseline(s,tracks,1e-10)
            output["original_native_q_and_daylength_mean_n_replayed_before_ERA5"]=True
            admitted,cells=orig.source_days(tracks[s["taxon"]])
            if not cells or len(cells)!=p["new_weather_source_grid_count"]:
                raise ValueError("Original source q/calendar weather cell support changed")
            output["weather_grid_count"]=len(cells)
            step="SAME_VERSION_ALL_SOURCE_ERA5_DAILY_PHYSICAL_SHA"
            weather={};errors=0
            for cell,dates in sorted(cells.items()):
                try:
                    records,receipt=daily_source_era5(cell,dates)
                    weather[cell]=records
                except Exception:errors+=1
            output["weather_grid_HOLD_count"]=errors
            if errors or len(weather)!=len(cells):
                raise ValueError("Old complete weather source has new missing grids; no partial month inference")
            current_sha=orig.physical_fingerprint(weather)
            if current_sha!=p["current_provider_daily_physical_fingerprint"]:
                raise ValueError("ERA5 physical snapshot not identical to previous original q+weather")
            output["original_weather_version_exact_20261116_physical_sha_replayed"]=True
            step="PRIOR_SIX_WEATHER_ENDPOINTS_N_MEAN_REPLAY"
            old=adjusted_daily_effects(admitted,weather)
            for key in CLIMATE_KEYS+("unadjusted_daylight_on_climate_days","daylight_given_tmax","daylight_given_antecedent_water"):
                v=old[key]
                oldv=p["source_weather_adjusted_endpoints"][key]
                if (v["admitted_animals"]!=oldv["original_source_tagged_animals"] or
                    v["admitted_animal_years"]!=oldv["source_original_animal_years"] or
                    (v["mean"] is None)!=(oldv["mean"] is None) or
                    (v["mean"] is not None and
                     not math.isclose(v["mean"],oldv["mean"],abs_tol=1e-10,rel_tol=0))):
                    raise ValueError("Old first daily climate source native original n/mean drifted")
            output["prior_six_climate_model_mean_n_replayed"]=True
            step="FIRST_ORIGINAL_TAG_YEAR_MONTH_JOINT_MODELS"
            new=per_source_same_sample_models(admitted,weather,seed=20261116+1000*index)
            output["source_month_joint"]=new
            output["status"]="SOURCE_MONTH_JOINT_ESTIMABLE" if new["eligible_original_tagged_animals_joint_all_three"]>=5 else "SOURCE_MONTH_JOINT_UNAVAILABLE"
            tracks.clear();weather.clear()
    except Exception as e:
        output["hold_stage"]=step
        output["hold_type"]=type(e).__name__
        output["source_month_joint"]=None
        output["status"]="SOURCE_HOLD"
    return output


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--index",type=int,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    r=run(args.index)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    e=(r["source_month_joint"] or {}).get("month_partial_joint_estimands",{})
    print("RSDS_THREE_MULE_MONTH_WEATHER",json.dumps({
        "index":args.index,"source_status":r["status"],"hold_stage":r["hold_stage"],
        "hold_type":r["hold_type"],"old_models_replayed":r["prior_six_climate_model_mean_n_replayed"],
        "qualified_independent_tags":(r["source_month_joint"] or {}).get("eligible_original_tagged_animals_joint_all_three"),
        "Tmax_month_CI":e.get("month_Tmax_given_L_water"),
        "daylength_month_CI":e.get("month_L_given_Tmax_water"),
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
