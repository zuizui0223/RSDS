"""Same-place AND same-calendar-date mule deer vs bighorn DAY/NIGHT heat response.

Source data: ONE original Utah 4-native Movebank GPS archive, SHA/MD5 frozen.
After original q/daylength and all four native source q/day GPS source
rosters replay, check exact 2026-11-17 daily ERA5 physical SOURCE SHA.

Pairs live on the SAME original 1-degree grid, calendar year-month and
calendar DATE. Weather per cell/date is physically the same. Estimates are
cell-equal block slopes with cell-cluster bootstrap, not tag/day pseudo
replicates. Post-outcome ecological sensitivity: no species adaptation claim.
No raw tag IDs, coordinates, exact dates or GPS records exported.
"""
from __future__ import annotations

import argparse
from collections import defaultdict,Counter
from datetime import date
import json
import math
from pathlib import Path
import random
from statistics import mean,pstdev
import tempfile

from scripts.rsds_437_systematic_source_frame_v1 import HERE, pinned, git_blob
from scripts import rsds_two_source_all_native_ungulate_month_heat_v1 as originals
from scripts import rsds_three_mule_original_daily_hydrothermal_v1 as climate
from scripts import rsds_three_mule_day_night_heat_components_v1 as components
from scripts.movebank_native_ecology_20261019 import (
    original_files,checksum_and_download,daylength
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.movebank_wyoming_daily_climate_20261023 import (
    daily_source_era5,climate_at_gps_day
)
from scripts.movebank_utah_calendar_condition_20261022 import residuals

CONTRACT=HERE/"research/RSDS_UTAH_SHARED_CELL_DATE_MULE_BIGHORN_THERMAL_COMPONENT_CONTRACT_20261118.json"
VERSION="RSDS_UTAH_SHARED_ORIGINAL_GPS_CELL_DAY_MULE_BIGHORN_THERMAL_COMPONENTS_V1"
OUT_VERSION="RSDS_UTAH_FIRST_ORIGINAL_SAME_CELL_SAME_DAY_MULE_BIGHORN_DAY_NIGHT_HEAT_V1"
MULE="Odocoileus hemionus"
SHEEP="Ovis canadensis"
RESPONSES=("q","day","night")
TOL=1e-10
NBOOT=2000


def frozen():
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    required={
        "version":VERSION,
        "source_DOI":"10.5441/001/1.711",
        "source_study_ID":"7364502758",
        "original_all_four_registered_native_taxa":[
            "Antilocapra americana","Cervus elaphus",
            "Odocoileus hemionus","Ovis canadensis"],
        "previous_six_native_source_result_blob":"9bf37249223b7312ee86c800cb2170adc422f996",
        "source_first_physical_ERA5_daily_SHA256":"0d1992f19c79aacaaf3535182dd0b7ea6ca805a88e2b67caf6721955e6fa4f69",
        "source_first_full_ERA5_grid_count":28,
        "prior_source_shared_grid_yearmonth_pairs_any_gps_tag":118,
        "prior_source_shared_grid_yearmonth_with_at_least_3_distinct_tags_each_taxon":95,
        "minimum_both_species_shared_original_dates_within_cell_yearmonth":5,
        "minimum_tagged_animals_per_taxon_in_each_candidate_grid_yearmonth":3,
        "within_calendar_cell_yearmonth_daylength_range_h_min":.05,
        "minimum_qualified_cell_months_per_one_degree_cell":2,
        "minimum_shared_source_days_per_qualified_cell":15,
        "minimum_centered_L_sxx_h2_per_cell":.05,
        "required_after_L_and_water_adjustment_Tmax_residual_sd_C":.25,
        "source_cell_min_admissible_independent_geographic_blocks_for_bootstrap":5,
        "original_source_first_q_GPS_MD5":"08200c9dbdf7989fb5a35467213b7e5a",
        "original_source_GPS_provider_bitstream_uuid":"9c4aaf76-7476-4305-8b14-9a6d8886e39a",
        "original_source_GPS_bytes":63369773,
        "original_20261022_q_mule_minus_bighorn_Tmax_with_L_water_was_NOT_nominally_supported":True,
        "selected_after_prior_source_positive_temperature_outcomes":True,
        "previous_48_signed_solar_heading_external_verdict":"NOT_REPLICATED",
        "previous_12_ERA5_raw_SHA_source_status":"HOLD",
    }
    if any(c.get(k)!=v for k,v in required.items()):
        raise ValueError("Original Utah same-grid/same-date science contract changed")
    for stem in ("original_source_GPS_native_q_first_code",
                 "original_utah_source_zip_reader_code",
                 "original_source_gps_cell_ERA5_code",
                 "original_source_q_date_admission_code",
                 "original_source_DAY_NIGHT_components_code",
                 "original_source_weather_code"):
        if stem=="original_source_GPS_native_q_first_code":
            codepath=c["original_source_GPS_native_q_first_code_path"]
            blob=c["original_source_GPS_native_q_first_code_blob"]
        else:
            codepath=c[stem+"_path"];blob=c[stem+"_blob"]
        if git_blob((HERE/codepath).read_bytes())!=blob:
            raise ValueError("Original source GPS/ERA5/solar code version drift: "+stem)
    previous=pinned(c["previous_six_native_source_result_path"],
                    c["previous_six_native_source_result_blob"])
    pinned(c["previous_six_native_first_original_contract_path"],
           c["previous_six_native_first_original_contract_blob"])
    pinned(c["original_utah_source_path"],c["original_utah_source_blob"])
    legacy=pinned(c["original_previous_legacy_12_native_result_path"],
                  c["original_previous_legacy_12_native_result_blob"])
    sourcecontract,sources,old=originals.frozen()
    ut=sources[0]
    if (ut["doi"]!=c["source_DOI"] or ut["id"]!=c["source_study_ID"]
        or ut["taxa"]!=c["original_all_four_registered_native_taxa"]
        or ut["uuid"]!=c["original_source_GPS_provider_bitstream_uuid"]
        or ut["md5"]!=c["original_source_first_q_GPS_MD5"]
        or ut["size"]!=c["original_source_GPS_bytes"]):
        raise ValueError("Source original Utah MD5/study ID or registered taxa changed")
    assert len(legacy["all_original_study_native_taxon_rows_including_HOLD"])==22
    if (previous["original_source_studies"]!=2 or
        previous["original_native_taxon_study_rows"]!=6 or
        previous["original_source_complete_tag_month_native_rows"]!=6):
        raise ValueError("Original six source native temperature results incomplete")
    uta=previous["original_source_statuses"][0]
    if (uta["source"]!=ut["doi"] or uta["status"]!="SOURCE_ALL_NATIVE_ERA5_MONTH_MODELS_EVALUATED"
        or uta["original_source_weather_grids"]!=28
        or uta["weather_HOLDs"]!=0 or
        uta["new_source_ERA5_fingerprint"]!=c["source_first_physical_ERA5_daily_SHA256"]):
        raise ValueError("Original Utah all-native physical ERA5 source snapshot differs")
    prior=[x for x in previous["all_original_six_source_native_rows_including_HOLD"]
           if x["source_original_DOI"]==ut["doi"]]
    if len(prior)!=4 or [x["native_taxon"] for x in prior]!=ut["taxa"]:
        raise ValueError("Original Utah source four native support changed")
    return c,ut,old[ut["doi"]],prior


def source_native_cell_calendar(source_rows):
    """Only original eligible source tag q DAY/NIGHT, never private IDs public."""
    record=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for native,rows in source_rows.items():
        for tag,year,days in rows:
            for d,L,q,cell,ym,D,N in days:
                if d.year!=year or ym!=(d.year,d.month):
                    raise ValueError("Original tagged source q-date/year-month mismatch")
                if abs((N-D)-q)>TOL:
                    raise ValueError("Original native q != NIGHT minus DAY on source tag/day")
                key=(cell,ym)
                rec=record[native][key][d]
                rec.append((tag,q,D,N))
    return record


def source_shared_day_blocks(src,weather,contract):
    mule=src[MULE];sheep=src[SHEEP]
    shared=set(mule)&set(sheep)
    if len(shared)!=contract["prior_source_shared_grid_yearmonth_pairs_any_gps_tag"]:
        raise ValueError("Previously frozen 118 Utah source shared grid/year-month opportunities changed")
    ntag=int(contract["minimum_tagged_animals_per_taxon_in_each_candidate_grid_yearmonth"])
    def tags_in_month(days):
        return {v[0] for records in days.values() for v in records}
    eligible=[key for key in sorted(shared)
              if len(tags_in_month(mule[key]))>=ntag and
                 len(tags_in_month(sheep[key]))>=ntag]
    if len(eligible)!=contract["prior_source_shared_grid_yearmonth_with_at_least_3_distinct_tags_each_taxon"]:
        raise ValueError("Previously frozen 95 both-native >=3 tagged animal grid/month frame changed")
    cellblocks=defaultdict(list)
    exact_days=0;month_support=0;missing_weather=0;ge3_same_dates=0
    for cell,ym in eligible:
        both=sorted(set(mule[(cell,ym)]) & set(sheep[(cell,ym)]))
        exact_days+=len(both)
        if len(both)<contract["minimum_both_species_shared_original_dates_within_cell_yearmonth"]:
            continue
        rows=[]
        for d in both:
            cl=climate_at_gps_day(weather.get(cell,{}),d)
            if cl is None:
                missing_weather+=1;continue
            if (len(mule[(cell,ym)][d])>=3 and len(sheep[(cell,ym)][d])>=3):
                ge3_same_dates+=1
            # Equal original-tag weights within each species on the SAME
            # calendar day, with precisely the same weather exposure.
            A=mule[(cell,ym)][d];B=sheep[(cell,ym)][d]
            diff=[mean(x[i] for x in A)-mean(x[i] for x in B) for i in (1,2,3)]
            L=daylength(cell[0]+.5,d)
            rows.append((float(cl[1]),float(L),float(cl[2]),*diff))
        if missing_weather:
            raise ValueError("Original exact same-date 2026 ERA5 provider complete climate lost")
        if (len(rows)<contract["minimum_both_species_shared_original_dates_within_cell_yearmonth"]
            or max(r[1] for r in rows)-min(r[1] for r in rows)
               <contract["within_calendar_cell_yearmonth_daylength_range_h_min"]):
            continue
        month_support+=1
        cellblocks[cell].append(rows)
    return {
        "source_both_native_original_shared_cell_yearmonth_total":len(shared),
        "source_both_native_original_at_least_three_tags_each_shared_cell_yearmonth":len(eligible),
        "original_exact_shared_tag_calendar_dates_across_all_candidate_cell_months":exact_days,
        "original_exact_shared_dates_with_at_least_three_tags_each_species":ge3_same_dates,
        "source_cell_yearmonth_with_min_same_original_GPS_dates":month_support,
        "source_exact_shared_date_weather_missing":missing_weather,
    },cellblocks


def cell_month_demeaned_temperature_slopes(blocks,contract):
    """Returns per-cell DAY/NIGHT/Q Tmax slope for same-day mule−bighorn means.

    A cell is the geographic cluster. Month FE is applied *within cell*; no
    GPS tagged/day pseudo-replication. All three responses use IDENTICAL dates.
    """
    if len(blocks)<contract["minimum_qualified_cell_months_per_one_degree_cell"]:
        return None
    selected=[row for block in blocks for row in block]
    if len(selected)<contract["minimum_shared_source_days_per_qualified_cell"]:
        return None
    centered=[]
    lxx=0.
    for block in blocks:
        mus=[mean(r[i] for r in block) for i in range(6)]
        centered.extend(tuple(r[i]-mus[i] for i in range(6)) for r in block)
        lxx+=sum((r[1]-mus[1])**2 for r in block)
    if lxx<contract["minimum_centered_L_sxx_h2_per_cell"]:
        return None
    # Tmax residualized against astronomical daylength and 30d climatic water
    # using previously frozen pure 2covariate-FWL helper.
    x=[r[0] for r in centered]
    co=[[r[j] for r in centered] for j in (1,2)]
    rx=residuals(x,co)
    if rx is None or pstdev(rx)<contract["required_after_L_and_water_adjustment_Tmax_residual_sd_C"]:
        return None
    denom=sum(v*v for v in rx)
    if denom<=1e-12:return None
    out={}
    for key,j in (("q",3),("day",4),("night",5)):
        ry=residuals([r[j] for r in centered],co)
        if ry is None:return None
        v=sum(a*b for a,b in zip(rx,ry))/denom
        if not math.isfinite(v):return None
        out[key]=v
    if abs(out["night"]-out["day"]-out["q"])>TOL:
        raise ValueError("Within Utah original source cell month species Q!=NIGHT-DAY")
    return out,len(selected),len(blocks)


def original_source_equal_cell_bootstrap(cellblocks,contract):
    supported={}
    eligible_block_months=eligible_dates=0
    for cell,blocks in cellblocks.items():
        x=cell_month_demeaned_temperature_slopes(blocks,contract)
        if x is None:continue
        supported[cell]=x[0]
        eligible_dates+=x[1];eligible_block_months+=x[2]
    n=len(supported)
    total=len(cellblocks)
    result={
        "original_geographic_cells_with_any_shared_month_and_enough_daily_records":total,
        "original_geographic_cells_with_identifiable_joint_Tmax_L_water":n,
        "original_common_cellmonth_blocks_in_final_cell_level_regression":eligible_block_months,
        "original_common_solar_dates_in_final_regressions_not_independent":eligible_dates,
        "source_shared_cell_day_species_components_estimates":{},
    }
    if n<contract["source_cell_min_admissible_independent_geographic_blocks_for_bootstrap"]:
        result["source_status"]="INSUFFICIENT_SHARED_CELL_MONTH_GEOSPATIAL_SUPPORT"
        for name in RESPONSES:
            result["source_shared_cell_day_species_components_estimates"][name]={
                "status":"UNAVAILABLE","mean":None,
                "original_one_degree_cell_cluster_bootstrap95":None,
                "source_cell_n":n,"nominal_sign":"UNAVAILABLE",
            }
        return result
    keys=sorted(supported)
    rng=random.Random(20261118)
    draws=[[keys[rng.randrange(n)] for _ in range(n)] for z in range(NBOOT)]
    first={}
    for name in RESPONSES:
        v=[supported[k][name] for k in keys]
        samples=sorted(mean(supported[k][name] for k in drawn) for drawn in draws)
        ci=[samples[49],samples[1949]]
        first[name]={
            "status":"ESTIMABLE",
            "mean":mean(v),
            "original_one_degree_cell_cluster_bootstrap95":ci,
            "source_cell_n":n,
            "nominal_sign":"POSITIVE" if ci[0]>0 else
                           "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO",
            "positive_cell_count":sum(x>0 for x in v),
            "negative_cell_count":sum(x<0 for x in v)
        }
    if abs(first["night"]["mean"]-first["day"]["mean"]-first["q"]["mean"])>TOL:
        raise ValueError("Source cell-equal original native NIGHT-DAY source regression identity failed")
    for sample in draws:
        delta=mean(supported[k]["night"]-supported[k]["day"]-supported[k]["q"] for k in sample)
        if abs(delta)>TOL:
            raise ValueError("Spatial cell bootstrap species Q=NIGHT−DAY exact draw failed")
    result["source_status"]="SHARED_ORIGINAL_CELL_DAY_MULE_BIGHORN_EFFECT_ESTIMABLE"
    result["source_shared_cell_day_species_components_estimates"]=first
    result["same_geographic_cell_bootstrap_draws_for_q_DAY_NIGHT"]=True
    result["original_effect_study_unit_geographic_cells_not_animals_or_days"]=True
    return result


def process():
    c,s,legacy,first=frozen()
    result={
        "version":OUT_VERSION,
        "source_original_DOI":s["doi"],
        "source_native_pair":[MULE,SHEEP],
        "source_all_registered_native_taxa":s["taxa"],
        "status":"SOURCE_HOLD",
        "source_hold_stage":None,"source_hold_type":None,
        "original_published_source_GPS_MD5_verified":False,
        "original_four_native_full_q_daylight_baseline_replayed":False,
        "original_q_equals_original_NIGHT_minus_DAY_daily_replayed":False,
        "original_20261117_all_native_ERA5_physical_SHA_verified":False,
        "source_old_full_weather_grids":0,"source_weather_missing_grids":0,
        "original_source_shared_cell_day_support":None,
        "original_source_shared_cell_day_block_model":None,
        "Wyoming_source_original_same_grid_month_candidate_count":0,
        "old_48_source_irradiance_heading_external_NOT_REPLICATED":True,
        "old_12_ERA5_raw_SHA_HOLD":True,
        "no_original_GPS_tags_precise_dates_coordinates_original_grid_or_individual_day_values_exported":True,
        "source_study_selected_POSTOUTCOME":True,
        "original_20261022_common_cellmonth_Tmax_species_contrast_nominally_includes_zero":True,
    }
    stage="ORIGINAL_UTAH_PROVIDER_UUID_BYTE_MD5"
    try:
        with tempfile.TemporaryDirectory(prefix="rsds_utah_shared_cell_day_private_") as directory:
            path=Path(directory)/"original.zip"
            matches=[x for x in original_files(s["package"]) if x["bitstream_uuid"]==s["uuid"]]
            if (len(matches)!=1 or matches[0]["bytes"]!=s["size"]
                or matches[0]["checksum"]["checkSumAlgorithm"].upper()!="MD5"
                or matches[0]["checksum"]["value"].lower()!=s["md5"]):
                raise ValueError("Utah source original GPS provider UUID/MD5/bytes mismatch")
            if not checksum_and_download(matches[0],path)["provider_checksum_verified"]:
                raise ValueError("Utah original GPS provider MD5 not verified")
            result["original_published_source_GPS_MD5_verified"]=True
            stage="ALL_FOUR_NATIVE_ORIGINAL_Q_DAYLIGHT_FIRST_SOURCE_REPLAY"
            tracks,counts,nonmatching=extract_time_only_original_source(path,s["taxa"])
            if not all(counts[taxon]["rows"]>0 for taxon in s["taxa"]):
                raise ValueError("Original Utah all four native source roster altered")
            admitted={};parts={};weather_cells=defaultdict(set)
            for taxon in s["taxa"]:
                source_effect=one_taxon(taxon,tracks[taxon])[0]
                original=legacy[taxon]
                for key,firstkey in (("q","day_night_speed_contrast"),
                                     ("daylength","daylength_association")):
                    old=original[key];new=source_effect[firstkey]
                    if new is None or old["source_tags"]!=new["n_individuals"] or \
                       not math.isclose(old["mean"],new["mean"],rel_tol=0,abs_tol=TOL):
                        raise ValueError("Utah original source native first q/daylength baseline not replayed")
                native_admitted,cells=climate.source_days(tracks[taxon])
                for cell,dates in cells.items():
                    weather_cells[cell].update(dates)
                admitted[taxon]=native_admitted
                if taxon in (MULE,SHEEP):
                    native_parts=components.original_source_components(tracks[taxon])
                    components.verify_original_daily_q_reconstruction(native_admitted,native_parts)
                    parts[taxon]=native_parts
            result["original_four_native_full_q_daylight_baseline_replayed"]=True
            result["original_q_equals_original_NIGHT_minus_DAY_daily_replayed"]=True
            stage="REVALIDATE_ORIGINAL_ALL_NATIVE_ERA5_VERSION"
            if len(weather_cells)!=c["source_first_full_ERA5_grid_count"]:
                raise ValueError("Original all Utah four-native q weather grids changed")
            result["source_old_full_weather_grids"]=len(weather_cells)
            weather={};errors=0
            for cell,dates in sorted(weather_cells.items()):
                try:
                    table,receipt=daily_source_era5(cell,dates)
                    weather[cell]=table
                except Exception:
                    errors+=1
            result["source_weather_missing_grids"]=errors
            if (errors or len(weather)!=len(weather_cells) or
                climate.physical_fingerprint(weather)!=c["source_first_physical_ERA5_daily_SHA256"]):
                raise ValueError("Original Utah six-native first ERA5 physical content changed")
            result["original_20261117_all_native_ERA5_physical_SHA_verified"]=True
            stage="FIRST_COMMON_SOURCE_CELL_AND_EXACT_CALENDAR_DAY_SUPPORT"
            cohort=source_native_cell_calendar(parts)
            support,cellblocks=source_shared_day_blocks(cohort,weather,c)
            result["original_source_shared_cell_day_support"]=support
            stage="SAME_ORIGINAL_DAY_GEO_CELL_MONTH_SPECIES_TMAX_COMPONENT_CONTRAST"
            pair=original_source_equal_cell_bootstrap(cellblocks,c)
            result["original_source_shared_cell_day_block_model"]=pair
            result["status"]=pair["source_status"]
            tracks.clear();weather.clear()
    except Exception as e:
        result["status"]="SOURCE_HOLD"
        result["source_hold_stage"]=stage
        result["source_hold_type"]=type(e).__name__
        result["original_source_shared_cell_day_block_model"]=None
        # Keep source identity/provenance flags only; no partial effect values.
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    result=process()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print("RSDS_UTAH_SAME_SOURCE_CELL_DAY",json.dumps({
        "source_status":result["status"],
        "hold_stage":result["source_hold_stage"],
        "hold_type":result["source_hold_type"],
        "original_first_ERA5_physical_replayed":
            result["original_20261117_all_native_ERA5_physical_SHA_verified"],
        "candidate_month_support":result["original_source_shared_cell_day_support"],
        "matched_cell_model":result["original_source_shared_cell_day_block_model"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
