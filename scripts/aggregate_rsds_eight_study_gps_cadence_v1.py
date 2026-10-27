"""Combine ALL eight provider-locked GPS observation-cadence original source results.

Original source DOI is one field study; within-study native species are NOT
independent ecosystems. This is a protocol/measurement sensitivity, not
ecological change over time or a causal daylight/radiation hypothesis.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    HERE,PAIRS,verify_frame,NEGATIVE
)

OUTPUT=HERE/"research/RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_FIRST_RESULT_20261028.json"


def aggregate(sources):
    old,new,fixed,metadata=verify_frame()
    if len(sources)!=8 or set(sources)!=set(range(8)):
        raise ValueError("Original eight field study source frame was silently reduced")
    original={
        0:("10.5441/001/1.711",set(x["taxon"] for x in old["sources"][0]["taxon_results"])),
        1:("10.5441/001/1.703",set(x["taxon"] for x in old["sources"][1]["taxon_results"])),
    }
    for i,x in enumerate(fixed,2):
        original[i]=(x["doi"],{x["taxon"]})
    study_status=Counter()
    support={p:{"matched_q":Counter(),"matched_daylength":Counter()} for p in PAIRS}
    flattened=[]
    source_records=[]
    source_equal={}
    for i in range(8):
        row=sources[i]
        expected_doi,taxa=original[i]
        if (row.get("source_index")!=i or row.get("original_study_doi")!=expected_doi
                or row.get("version")!="RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_V1_SOURCE"
                or set(row.get("source_native_taxa",[]))!=taxa
                or row.get("unchanged_48_source_signed_radiation_external_NONREPLICATION")!=NEGATIVE
                or row.get("unchanged_old_ERA5_12_RAW_SHA_HOLDS") is not True):
            raise ValueError("Original source identity/frame/previous negative result changed")
        state=row["source_status"]
        study_status[state]+=1
        entries=row["native_taxon_results"]
        if len(entries)!=len(taxa) or {t["taxon"] for t in entries}!=taxa:
            raise ValueError("One original native source taxon was dropped or duplicated")
        if state=="SOURCE_CADENCE_SENSITIVITY_COMPLETED" and (
                not row["original_GPS_provider_MD5_verified"]
                or not row["baseline_1h_q_and_daylength_source_effects_replayed"]):
            raise ValueError("Original q or source MD5 not verified before cadence response")
        entry_summaries=[]
        for t in entries:
            phase=t.get("source_original_even_odd_2h_phase_sensitivity")
            if state=="SOURCE_CADENCE_SENSITIVITY_COMPLETED" and (
                    not isinstance(phase,dict) or set(phase)!=set(PAIRS)):
                raise ValueError("Missing original even/odd UTC 2h phase for an eligible source")
            if not isinstance(phase,dict):
                phase=None
            if phase:
                for mode in PAIRS:
                    effects=phase[mode]
                    p=effects["matched_q_velocity_proxy_difference"]
                    l=effects["matched_q_daylength_slope_difference_per_hour"]
                    if p["n_original_tagged_animals_paired_between_cadences"]<5:
                        if p["paired_tag_mean_2h_minus_same_day_1h"] is not None:
                            raise ValueError("Unestimated matched original tag q became a numeric zero")
                    for axis,z in (("matched_q",p),("matched_daylength",l)):
                        if z["status"] not in (
                                "SOURCE_TAGGED_ANIMAL_DESCRIPTIVE","INSUFFICIENT_MATCHED_TAGS"):
                            raise ValueError("Unknown original source parity support status")
                        if z["status"]=="SOURCE_TAGGED_ANIMAL_DESCRIPTIVE":
                            if (z["n_original_tagged_animals_paired_between_cadences"]<5
                                    or not isinstance(z["original_tag_resample95"],list)
                                    or len(z["original_tag_resample95"])!=2):
                                raise ValueError("Cadence effect without required independent tagged animals")
                            low,high=z["original_tag_resample95"]
                            m=z["paired_tag_mean_2h_minus_same_day_1h"]
                            if not all(math.isfinite(v) for v in (low,high,m)) or not low<=m<=high:
                                raise ValueError("Original source matched-tag interval invalid")
                        support[mode][axis][z["nominal_interval_sign"]]+=1
                if t["source_original_full_hour_baseline_support"]["tags_with_some_q_observable"]<0:
                    raise ValueError("Negative native GPS observed-tag count")
            entry_summaries.append({"taxon":t["taxon"],"original_source_status":state,
                                    "original_hourly_GPS_data_support":
                                      t.get("source_original_full_hour_baseline_support"),
                                    "even_odd_cadence_effects":phase})
            flattened.append(entry_summaries[-1])
        source_records.append({
            "original_study_doi":expected_doi,
            "fixed_original_study_index":i,
            "source_status":state,
            "native_taxon_count":len(taxa),
            "source_native_taxon_results":entry_summaries,
            "original_MD5_accepted":row["original_GPS_provider_MD5_verified"],
            "fullhour_q_daylength_source_replayed":
                row["baseline_1h_q_and_daylength_source_effects_replayed"],
            "source_hold_reason":row.get("source_hold_error_type") or row.get("hold_reason")
        })
        source_equal[expected_doi]={}
        for mode in PAIRS:
            valid=[t["even_odd_cadence_effects"][mode]["matched_q_velocity_proxy_difference"][
                       "paired_tag_mean_2h_minus_same_day_1h"]
                   for t in entry_summaries if t["even_odd_cadence_effects"] is not None
                   and t["even_odd_cadence_effects"][mode]["matched_q_velocity_proxy_difference"][
                       "status"]=="SOURCE_TAGGED_ANIMAL_DESCRIPTIVE"]
            source_equal[expected_doi][mode]=mean(valid) if valid else None
    if len(flattened)!=12 or len(source_records)!=8:
        raise ValueError("Original eight source DOI/twelve native taxon denominator was lost")
    return {
        "version":"RSDS_EIGHT_STUDY_UTC_GPS_1H_VS_2H_MATCHED_DAY_V1_FIRST_RESULT",
        "date":"2026-10-28",
        "analysis_stage":"POST_ORIGINAL_OUTCOMES_SAME_TAG_SAME_DAY_TELEMETRY_CADENCE_SENSITIVITY",
        "original_field_studies":8,
        "original_source_taxon_rows":12,
        "status_by_original_source_study":dict(study_status),
        "original_q_paired_cadence_support_counts":{
            phase:{k:dict(val) for k,val in support[phase].items()} for phase in PAIRS
        },
        "all_original_source_study_and_native_taxa":source_records,
        "one_original_DOI_one_study_phase_delta_mean":source_equal,
        "historic_ERA5_old_12_source_raw_SHA_HOLDS_UNCHANGED":True,
        "historical_48_source_negative_radiation_replication_UNCHANGED":NEGATIVE,
        "limits":[
            "A 2h movement straight-line speed difference from 1h is inherently sensitive to travel tortuosity, measured time intervals and diurnal edges; do not interpret as a true activity budget change or animal response to sampling.",
            "2h phases are even/odd UTC hour bins, not morning/afternoon, true solar period exposure, or independently recruited individuals.",
            "Within each source tagged animal/year the original dates eligible under 1h and 2h are intersected. Numbers of matching dates, individuals, slopes and source holds are reported and no missing effect replaced by zero.",
            "The baseline 1h q and 2h source effects are different observation estimands even on the same animal/date because displacement chord length varies with interval length.",
            "Native taxa within Utah and Wyoming share their underlying original studies. Study-equal source summaries cannot be treated as independent native taxa, random world ecosystem samples or causal within-study environmental effects.",
            "All tag-resampled 95% intervals are nominal descriptive source-conditional uncertainties uncorrected for multiplicity.",
            "No ERA5 requested. Earlier original GPS-linked ERA5 12/12 raw SHA version holds and external 48-study directional shortwave negative-replication remain unchanged.",
            "No original animal tag IDs, raw GPS location/time, private source grids, per-animal slope or restricted records are redistributed."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,default=OUTPUT)
    a=p.parse_args()
    inputs={}
    for i in range(8):
        path=a.input/("original_study_"+str(i)+"_cadence_v1.json")
        if not path.is_file():
            raise ValueError("One fixed original study is absent; do not replace: "+str(i))
        inputs[i]=json.loads(path.read_text(encoding="utf-8"))
    d=aggregate(inputs)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                     encoding="utf-8")
    print("RSDS_EIGHT_STUDY_CADENCE",json.dumps({
        "n_original_studies":d["original_field_studies"],
        "n_native_study_taxon":d["original_source_taxon_rows"],
        "source_states":d["status_by_original_source_study"],
        "source_effect_support":d["original_q_paired_cadence_support_counts"],
        "historic_source_HOLDS_unchanged":d["historic_ERA5_old_12_source_raw_SHA_HOLDS_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
