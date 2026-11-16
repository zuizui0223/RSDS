"""Original mule GPS equal-three DAY/NIGHT interval-opportunity robustness.

No GPS tag/day observations or ERA5 grid/dated weather leave the ephemeral
source Action. Exact previously frozen q, original GPS tag/date, first ERA5
physical SHA, prior six climate effects and first q-Tmax 2000 tag bootstraps
must replay on the same cohort before an equal-three-step diagnostic counts.
This is post-discovery sampling sensitivity, not a causal heat test.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import random
from statistics import mean,median
import tempfile

from scripts.rsds_437_systematic_source_frame_v1 import HERE,pinned,git_blob
from scripts import rsds_three_mule_day_night_heat_components_v1 as original
from scripts import rsds_three_mule_within_month_joint_weather_v1 as old
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as weather_origin
from scripts.movebank_native_ecology_20261019 import (
    daylength,original_files,checksum_and_download
)
from scripts.movebank_solar_movement_20261019 import valid_interval
from scripts.global_solar_azimuth_12_new_original_events_v1 import load_source_taxa
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.rsds_five_source_native_tag_month_conditioned_v1 import replay_original_baseline
from scripts.movebank_wyoming_daily_climate_20261023 import (
    daily_source_era5,adjusted_daily_effects,climate_at_gps_day
)

CONTRACT=HERE/"research/RSDS_THREE_MULE_BALANCED_DAY_NIGHT_OPPORTUNITY_CONTRACT_20261116.json"
VERSION="RSDS_THREE_MULE_SAME_TAG_DATE_EQUAL_GPS_DAY_NIGHT_OPPORTUNITY_V1"
SOURCE_VERSION="RSDS_THREE_MULE_ORIGINAL_SOURCE_EQUAL_THREE_DAY_NIGHT_GPS_FIRST_V1"
NBOOT=2000
TOL=1e-10
KEYS=(
"old_month_DAY_Tmax","old_month_NIGHT_Tmax","old_month_Q_Tmax",
"balanced_month_DAY_Tmax","balanced_month_NIGHT_Tmax","balanced_month_Q_Tmax",
"balanced_minus_old_month_DAY_Tmax","balanced_minus_old_month_NIGHT_Tmax",
"balanced_minus_old_month_Q_Tmax")

def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    expected={
        "version":VERSION,
        "prior_three_source_day_night_first_result_blob":"804bcbfcb3ab4ddf8fc1ea2e3d932ee00dba96fd",
        "prior_three_source_day_night_contract_blob":"a604bcac05e64b14647a7568b73764f2a923db32",
        "prior_three_source_day_night_code_blob":"1d64c213302ed2914ab374f2b8cd31eedfd95882",
        "prior_three_source_month_joint_code_blob":"42cfb74f8a6542b46c4e3158908fa207d884ef3b",
        "original_DOIs":["10.5441/001/1.711","10.5441/001/1.703","10.5441/001/1.708"],
        "original_species":"Odocoileus hemionus",
        "previous_original_month_cohort_tags":[1596,109,129],
        "previous_month_Tmax_nominal_signs":["POSITIVE","POSITIVE","POSITIVE"],
        "selected_day_night_step_count_per_tag_date":3,
        "source_independence_verified":False if "source_independence_verified" in c else None,
        "same_source_weather_sha_required":True if "same_source_weather_sha_required" in c else None,
        "old_48_source_irradiance_heading_external_direction_signed_result":"NOT_REPLICATED",
        "old_12_ERA5_raw_SHA_hold":"HOLD",
    }
    for k,v in expected.items():
        if k in ("source_independence_verified","same_source_weather_sha_required") and v is None:
            continue
        if c.get(k)!=v:
            raise ValueError("Frozen original source equal-opportunity scientific contract changed: "+k)
    p=pinned(c["prior_three_source_day_night_first_result_path"],
             c["prior_three_source_day_night_first_result_blob"])
    pinned(c["prior_three_source_day_night_contract_path"],
           c["prior_three_source_day_night_contract_blob"])
    for k in ("prior_three_source_day_night_code","prior_three_source_month_joint_code"):
        if git_blob((HERE/c[k+"_path"]).read_bytes())!=c[k+"_blob"]:
            raise ValueError("Original previous DAY/NIGHT source code drift")
    spec,sources,old_source,previous=original.frozen()
    if (p["registered_original_source_DOI_studies"]!=3 or
        p["original_source_status_counts"]!={"SOURCE_COMPONENTS_ESTIMABLE":3} or
        len(p["all_three_original_source_studies_including_SOURCE_HOLD"])!=3):
        raise ValueError("Original three-source before-balanced result changed")
    for i,(s,oldfirst,prev) in enumerate(zip(sources,old_source,
             p["all_three_original_source_studies_including_SOURCE_HOLD"])):
        if (s["doi"]!=c["original_DOIs"][i] or
            prev["source_original_DOI"]!=s["doi"] or
            prev["source_status"]!="SOURCE_COMPONENTS_ESTIMABLE" or
            prev["source_original_individual_tag_n"]!=c["previous_original_month_cohort_tags"][i] or
            prev["same_original_individual_tag_day_month_Q_Tmax"]["nominal_sign"]!="POSITIVE" or
            prev["source_original_q_daylength_model_replayed"] is not True or
            prev["first_2026_11_16_ERA5_physical_fingerprint_replayed"] is not True):
            raise ValueError("Original three source DAY/NIGHT q thermal pattern changed")
    return c,sources,old_source,p["all_three_original_source_studies_including_SOURCE_HOLD"]


def choose_three_by_original_interval_order(speeds):
    """Deterministic chronological 0/middle/last indices, NEVER sorted by speed."""
    if len(speeds)<3:raise ValueError("Original solar-date DAY/NIGHT >=3 not met")
    return [speeds[0],speeds[(len(speeds)-1)//2],speeds[-1]]


def original_balanced_date_rows(animals):
    """Private all original admitted day/night medians and equal 3-step medians."""
    output=[]
    for tag,source in animals.items():
        points=sorted(source.values())
        day_data=defaultdict(lambda:{"day":[],"night":[],"lat":[],"lon":[]})
        for a,b in zip(points,points[1:]):
            interval=valid_interval(a[0],a[1],a[2],b[0],b[1],b[2])
            if interval is None:continue
            d,part,kmph,hours,lat=interval
            rec=day_data[d]
            rec[part].append(kmph);rec["lat"].append(lat)
            rec["lon"].append((a[2]+b[2])/2)
        years=defaultdict(list)
        for day,z in day_data.items():
            if len(z["day"])<3 or len(z["night"])<3:continue
            lat=mean(z["lat"]);lon=mean(z["lon"])
            day_old=math.log1p(median(z["day"]))
            night_old=math.log1p(median(z["night"]))
            day_three=math.log1p(median(choose_three_by_original_interval_order(z["day"])))
            night_three=math.log1p(median(choose_three_by_original_interval_order(z["night"])))
            years[day.year].append((
                day,daylength(lat,day),night_old-day_old,
                (math.floor(lat),math.floor(lon)),(day.year,day.month),
                day_old,night_old,
                day_three,night_three,len(z["day"]),len(z["night"])
            ))
        for year,days in years.items():
            if len(days)>=10:output.append((tag,year,days))
    return output


def balanced_month_model(rows):
    """Previously q-qualified source models determine admission, not qbalanced."""
    baseline=original.matched_year_component_effect([
        (d,t,w,L,q,day,night)
        for d,t,w,L,q,day,night,bd,bn,nDay,nNight in rows
    ])
    if baseline is None:return None
    month=defaultdict(list)
    for d,t,w,L,q,day,night,bd,bn,nDay,nNight in rows:
        if abs((night-day)-q)>TOL:raise ValueError("Original day q replay failed")
        if nDay<3 or nNight<3:raise ValueError("Original solar daypart minimum segment support changed")
        month[d.month].append((
            float(t),float(L),float(w),
            float(q),float(day),float(night),
            float(bn-bd),float(bd),float(bn),nDay,nNight))
    blocks=[]
    for mo,g in sorted(month.items()):
        if len(g)<5 or max(x[1] for x in g)-min(x[1] for x in g)<.05:continue
        blocks.append(g)
    if len(blocks)<2:
        raise ValueError("Previously q-admitted source tag-year month group missing")
    def model(records,idx):
        # The original old month q estimator uses same positive SOURCE model
        # all three predictors Tmax,L,P-ET0/100 and FWL residual checks.
        return original.heat_partial(records,idx)
    centered=[]
    for group in blocks:
        bar=[mean(z[i] for z in group) for i in range(9)]
        centered.extend(tuple(z[i]-bar[i] for i in range(9)) for z in group)
    for i,attr in ((3,"month_q_tmax"),(4,"month_day_tmax"),(5,"month_night_tmax")):
        z=model(centered,i)
        if z is None or abs(z-baseline[attr])>TOL:
            raise ValueError("Original q, DAY, NIGHT source month FWL model did not replay")
    balanced={}
    for idx,key in ((6,"Q"),(7,"DAY"),(8,"NIGHT")):
        z=model(centered,idx)
        if z is None or not math.isfinite(z):
            raise ValueError("Equal-three q/part new slope not finite")
        balanced[key]=z
    if abs(balanced["NIGHT"]-balanced["DAY"]-balanced["Q"])>TOL:
        raise ValueError("Equal-three per-tag-year Q != NIGHT - DAY")
    values={
        "old_month_DAY_Tmax":baseline["month_day_tmax"],
        "old_month_NIGHT_Tmax":baseline["month_night_tmax"],
        "old_month_Q_Tmax":baseline["month_q_tmax"],
        "balanced_month_DAY_Tmax":balanced["DAY"],
        "balanced_month_NIGHT_Tmax":balanced["NIGHT"],
        "balanced_month_Q_Tmax":balanced["Q"]
    }
    values.update({
        "balanced_minus_old_month_DAY_Tmax":
            values["balanced_month_DAY_Tmax"]-values["old_month_DAY_Tmax"],
        "balanced_minus_old_month_NIGHT_Tmax":
            values["balanced_month_NIGHT_Tmax"]-values["old_month_NIGHT_Tmax"],
        "balanced_minus_old_month_Q_Tmax":
            values["balanced_month_Q_Tmax"]-values["old_month_Q_Tmax"]
    })
    day_records=sum(len(b) for b in blocks)
    unequal=sum(z[9]!=z[10] for block in blocks for z in block)
    return values,day_records,unequal


def source_tag_bootstrap(balanced_rows,weather,prior,index):
    by_tag=defaultdict(dict)
    total_candidate_days=0
    for tag,year,dates in balanced_rows:
        if year in by_tag[tag]:raise ValueError("Repeated original tag/year")
        ready=[]
        for d,L,q,cell,ym,day,night,bd,bn,nd,nn in dates:
            climate=climate_at_gps_day(weather.get(cell,{}),d)
            if climate is None:continue
            Tmean,tmax,water=climate
            ready.append((d,tmax,water,L,q,day,night,bd,bn,nd,nn))
        total_candidate_days+=len(ready)
        by_tag[tag][year]=ready
    cohort={}
    eligible_years=eligible_days=unbalanced_days=0
    for tag,years in by_tag.items():
        included=[]
        for yr,rows in years.items():
            r=balanced_month_model(rows)
            if r is not None:
                values,days,unequal=r
                eligible_years+=1
                eligible_days+=days
                unbalanced_days+=unequal
                included.append(values)
        if included:
            cohort[tag]={k:mean(r[k] for r in included) for k in KEYS}
    n=len(cohort)
    previous=prior["same_original_individual_tag_day_month_component_regression_means_and_paired_95"]
    if n!=prior["source_original_individual_tag_n"] or n!=previous["eligible_original_tagged_animals"]:
        raise ValueError("New equal-three diagnostic changed eligible original q/Tmax animal cohort")
    if n<5:raise ValueError("Previously supported source tagged animal q cohort vanished")
    # Replay exact previous q source-tag bootstrap, with identical original
    # random seed and sorted original tag labels.
    tags=sorted(cohort)
    rng=random.Random(20261116+1000*index)
    draws=[[tags[rng.randrange(n)] for _ in range(n)] for _ in range(NBOOT)]
    computed={}
    for key in KEYS:
        values=[cohort[tag][key] for tag in tags]
        samples=sorted(mean(cohort[t][key] for t in idxs) for idxs in draws)
        ci=[samples[49],samples[1949]]
        computed[key]={
            "status":"ESTIMABLE",
            "source_tag_n":n,
            "mean":mean(values),
            "tag_bootstrap95":ci,
            "nominal_sign":
                "POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO"
        }
    for key,previd in (("old_month_DAY_Tmax","month_day_tmax"),
                       ("old_month_NIGHT_Tmax","month_night_tmax"),
                       ("old_month_Q_Tmax","month_q_tmax")):
        old=previous["components"][previd]
        now=computed[key]
        if (old["original_tag_n"]!=n or
            abs(old["estimate"]-now["mean"])>TOL or
            any(abs(a-b)>TOL for a,b in zip(old["tag_bootstrap95"],now["tag_bootstrap95"]))):
            raise ValueError("Previous original Q/DAY/NIGHT source tag bootstrap not reproduced")
    for label in ("old_month","balanced_month"):
        vals=computed
        if abs(vals[label+"_NIGHT_Tmax"]["mean"]-
               vals[label+"_DAY_Tmax"]["mean"]-
               vals[label+"_Q_Tmax"]["mean"])>TOL:
            raise ValueError("Original/balanced study source Q != NIGHT - DAY")
    return {
        "original_q_weather_month_tagged_animals_identical_to_first":n,
        "source_q_model_qualified_original_tag_years":eligible_years,
        "original_tag_solar_dates_in_admitted_models_not_independent":eligible_days,
        "original_tag_solar_dates_with_unequal_day_night_valid_GPS_step_counts":unbalanced_days,
        "original_tag_solar_dates_with_unequal_count_fraction":
            unbalanced_days/eligible_days if eligible_days else None,
        "balanced_rule_three_accepted_intervals_per_daypart_per_tag_solar_date":True,
        "original_source_daily_q_and_first_heat_bootstrap_exact_replay":True,
        "equal_three_new_q_NIGHT_minus_DAY_algebra_exact":True,
        "bootstrap_same_original_tag_draws_all_endpoints":True,
        "source_effects":computed,
    }


def process(index):
    c,sources,oldsource,prior=frozen()
    if type(index)!=int or index not in (0,1,2):
        raise ValueError("Only three original frozen source studies")
    s=sources[index];p=prior[index];md=s["original_file"]
    out={
        "version":SOURCE_VERSION,
        "source_index":index,"source_key":s["source_key"],
        "original_study_DOI":s["doi"],"original_study_ID":str(s["study_id"]),
        "source_native":"Odocoileus hemionus",
        "original_registered_native_taxa":s["registered_taxa"],
        "status":"SOURCE_HOLD","hold_stage":None,"hold_type":None,
        "provider_source_MD5_verified":False,
        "source_native_full_q_and_daylength_replayed":False,
        "original_source_q_qualified_tag_solar_date_set_reconstructed":False,
        "exact_first_ERA5_daily_physical_SHA_verified":False,
        "previous_original_source_six_weather_model_means_replayed":False,
        "original_first_three_source_DAY_NIGHT_Q_tag_bootstraps_replayed":False,
        "source_weather_grids":0,"source_weather_grids_HOLD":0,
        "result":None,
        "source_48_irradiance_heading_direction_external_sign_NOT_REPLICATED":True,
        "source_12_old_ERA5_raw_SHA_HOLD":True,
        "no_original_tag_ID_GPS_exact_date_climate_cell_or_individual_speeds_exported":True,
    }
    ext=".zip" if md["filename"].lower().endswith(".zip") else (
         ".csv.gz" if md["filename"].lower().endswith(".csv.gz") else ".csv")
    stage="SOURCE_ORIGINAL_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_original_equal_step_private_") as tmp:
            dest=Path(tmp)/("original"+ext)
            rows=[r for r in original_files(s["package_uuid"])
                  if r["bitstream_uuid"]==md["uuid"]]
            if (len(rows)!=1 or rows[0]["bytes"]!=md["bytes"] or
                rows[0]["checksum"]["checkSumAlgorithm"].upper()!="MD5" or
                rows[0]["checksum"]["value"].lower()!=md["md5"].lower()):
                raise ValueError("Original source UUID/size/MD5 mismatch")
            if not checksum_and_download(rows[0],dest)["provider_checksum_verified"]:
                raise ValueError("Original source MD5 check failed")
            out["provider_source_MD5_verified"]=True
            stage="ORIGINAL_GPS_Q_REPLAY_WITH_EQUAL_INTERVAL_ROSTER"
            if s["source_key"]=="mule_old_utah":
                tracks,counts,nonmatching=extract_time_only_original_source(
                    dest,s["registered_taxa"])
            else:
                tracks,coverage,records=load_source_taxa(
                    dest,md["filename"],s["registered_taxa"])
            replay_original_baseline(s,tracks,1e-10)
            out["source_native_full_q_and_daylength_replayed"]=True
            original_admitted,cells=weather_origin.source_days(tracks[s["taxon"]])
            balanced=original_balanced_date_rows(tracks[s["taxon"]])
            audit=original.verify_original_daily_q_reconstruction(
                original_admitted,[(tag,year,[r[:7] for r in rs])
                                  for tag,year,rs in balanced])
            out["original_source_q_qualified_tag_solar_date_set_reconstructed"]=True
            out["original_q_algebra_source_diagnostics"]=audit
            stage="IDENTICAL_PREVIOUS_ERA5_SOURCE_PHYSCIAL_SHA"
            if len(cells)!=oldsource[index]["new_weather_source_grid_count"]:
                raise ValueError("Source eligible climate grid count changed")
            out["source_weather_grids"]=len(cells)
            weather={};holds=0
            for cell,dates in sorted(cells.items()):
                try:
                    x,receipt=daily_source_era5(cell,dates);weather[cell]=x
                except Exception:
                    holds+=1
            out["source_weather_grids_HOLD"]=holds
            if (holds or len(weather)!=len(cells) or
                weather_origin.physical_fingerprint(weather)!=
                    oldsource[index]["current_provider_daily_physical_fingerprint"]):
                raise ValueError("Old ERA5 canonical physical content source not exact")
            out["exact_first_ERA5_daily_physical_SHA_verified"]=True
            stage="PRIOR_SIX_SOURCE_WEATHER_MODELS"
            first=adjusted_daily_effects(original_admitted,weather)
            for name,oldresult in oldsource[index]["source_weather_adjusted_endpoints"].items():
                current=first[name]
                if (current["admitted_animals"]!=oldresult["original_source_tagged_animals"] or
                    current["admitted_animal_years"]!=oldresult["source_original_animal_years"] or
                    (current["mean"] is None)!=(oldresult["mean"] is None) or
                    (current["mean"] is not None and
                     not math.isclose(current["mean"],oldresult["mean"],rel_tol=0,abs_tol=TOL))):
                    raise ValueError("Old six original q-source climate models changed")
            out["previous_original_source_six_weather_model_means_replayed"]=True
            stage="NEW_EQUAL_THREE_GPS_VALID_STEP_OPPORTUNITY"
            new=source_tag_bootstrap(balanced,weather,p,index)
            out["original_first_three_source_DAY_NIGHT_Q_tag_bootstraps_replayed"]=True
            out["result"]=new
            out["status"]="SOURCE_BALANCED_3_EACH_COMPONENT_ESTIMABLE"
            tracks.clear();weather.clear()
    except Exception as exc:
        out["status"]="SOURCE_HOLD"
        out["hold_stage"]=stage
        out["hold_type"]=type(exc).__name__
        out["result"]=None
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--index",required=True,type=int)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    r=process(a.index)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                     encoding="utf-8")
    output=(r["result"] or {}).get("source_effects") or {}
    print("RSDS_THREE_MULE_EQUAL_DAY_NIGHT_GPS",json.dumps({
        "source":a.index,"status":r["status"],
        "hold_stage":r["hold_stage"],"hold_type":r["hold_type"],
        "original_q_replayed":r["original_first_three_source_DAY_NIGHT_Q_tag_bootstraps_replayed"],
        "original_source_tags":(r["result"] or {}).get(
            "original_q_weather_month_tagged_animals_identical_to_first"),
        "proportion_unequal_day_night_interval_counts":(r["result"] or {}).get(
            "original_tag_solar_dates_with_unequal_count_fraction"),
        "balanced_heat_q":output.get("balanced_month_Q_Tmax"),
        "balanced_heat_day":output.get("balanced_month_DAY_Tmax"),
        "balanced_heat_night":output.get("balanced_month_NIGHT_Tmax"),
    },sort_keys=True),flush=True)


if __name__=="__main__":main()
