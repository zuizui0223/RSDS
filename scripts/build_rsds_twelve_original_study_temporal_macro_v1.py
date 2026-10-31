"""RSDS original twelve-field-study temporal q and daylight macro evidence.

All numerical source inputs are PREVIOUSLY MD5-validated first original-data
study-group results, pinned to Git blob identities. No GPS event, tag ID,
private location or hourly weather is re-fetched or released here.

Within a DOI average available native taxa equally, then weight each of
the NINE supported original field studies once. Preserve all TWELVE study
DOIs and TWENTY-TWO source×native-taxon rows, including every original HOLD.
This is descriptive, purposively selected, post-first-result ecology, NOT
a global random-effects estimator or independent causal photoperiod test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"research/RSDS_12_STUDY_ORIGINAL_GPS_TEMPORAL_MACRO_V1_CONTRACT.json"
OUTPUT=ROOT/"research/RSDS_12_STUDY_22_NATIVE_TEMPORAL_MACRO_FIRST_RESULT_20261031.json"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
FIELDS=("q","daylength")
BOOTSTRAP_SAMPLES=5000
RANDOM_SEED=20261031


def pin_file(path,sha):
    raw=(ROOT/path).read_bytes()
    real=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
    if real!=sha:
        raise ValueError("Frozen original group scientific result changed: "+path)
    return json.loads(raw.decode("utf-8"))


def input_contract():
    spec=json.loads(CONTRACT.read_text(encoding="utf-8"))
    if spec["version"]!="RSDS_ORIGINAL_TWELVE_STUDY_TEMPORAL_MACRO_V1":
        raise ValueError("Unexpected post-first-result original study macro route")
    x={p["path"]:pin_file(p["path"],p["git_blob_sha"]) for p in spec["source_blobs"]}
    def get(substring):
        matches=[z for p,z in x.items() if substring in p]
        if len(matches)!=1:
            raise ValueError("Missing or ambiguous frozen original study-group source "+substring)
        return matches[0]
    return spec,get("8_STUDY_TIME_SOLAR"),get("FOUR_NEW_MULTINATIVE_GPS_SHARED"),\
        get("FOUR_MULTINATIVE_ORIGINAL_GPS_ECOLOGICAL_VALIDATION"),\
        get("FOUR_NEW_MULTINATIVE_SOURCE_METADATA")


def validated(n,point,interval,source_status):
    if point is None:
        if interval is not None:
            raise ValueError("Unavailable original source effect was given an interval")
        return {"status":"UNAVAILABLE","mean":None,"source_tags":n,"tag_bootstrap95":None,
                "origin_source_status":source_status}
    if not isinstance(n,int) or n<5 or interval is None or len(interval)!=2:
        raise ValueError("Original non-HOLD native effect lacks five tags and a valid interval")
    if any(type(z) not in (float,int) or not math.isfinite(z)
           for z in (point,interval[0],interval[1])):
        raise ValueError("Nonfinite original source q or photoperiod")
    low,high=interval
    if not low<=point<=high:
        raise ValueError("Native original source mean conflicts with its source-tag interval")
    sign=("POSITIVE" if low>0 else "NEGATIVE" if high<0 else "INCLUDES_ZERO")
    return {"status":sign,"mean":point,"source_tags":n,"tag_bootstrap95":[low,high],
            "origin_source_status":source_status}


def original_eight(old):
    if (old["original_study_denominator"]!=8 or old["source_taxon_rows"]!=12
            or old["source_taxa_with_q"]!=11 or old["source_taxa_with_solar_raw"]!=10
            or old["old_48_radiation_study_external_negative_direction_test_UNCHANGED"]!=NEGATIVE
            or old["old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED"] is not True):
        raise ValueError("Original eight field-study science frame changed")
    out=[]
    for r in old["all_sources_including_egret_original_metadata_HOLD"]:
        out.append({
            "original_study_doi":r["source_study_doi"],
            "native_taxon":r["taxon"],
            "original_source_stratum":"PREVIOUS_EIGHT_STUDIES",
            "original_source_status":r["source_status"],
            "q":validated(r["q"]["source_animals"],r["q"]["mean"],r["q"]["animal_bootstrap95"],
                          r["source_status"]),
            "daylength":validated(r["daylength_slope"]["source_animals"],
                                  r["daylength_slope"]["mean"],
                                  r["daylength_slope"]["animal_bootstrap95"],r["source_status"])
        })
    return out


def original_four(new,validation,metadata):
    if (new["fixed_original_published_studies"]!=4 or
            new["registered_original_study_native_taxon_rows"]!=10 or
            new["prior_48_source_shortwave_direction_NEGATIVE_EXTERNAL_REPLICATION_UNCHANGED"]!=NEGATIVE or
            not new["prior_12_ERA5_version_source_HOLDS_UNCHANGED"] or
            validation["source_original_MD5_verified"]!=4 or
            validation["new_geese_original_same_source_calendar_days"]!=120 or
            validation["new_geese_source_year_month_blocks"]!=5 or
            validation["previous_48_source_signed_shortwave_direction_EXTERNAL_NONREPLICATED_UNCHANGED"]!=NEGATIVE):
        raise ValueError("New original four study source MD5/shared calendar/legacy denial changed")
    if len(new["all_frozen_source_study_original_results"])!=4:
        raise ValueError("All four pre-source-registered field studies must remain")
    m=metadata.get("source_entries") or metadata.get("original_source_entries")
    if m is None:
        # Metadata pin + scientific validation are still independently checked;
        # schema and native names checked against first original source result.
        if metadata.get("version") is None:
            raise ValueError("Missing original pre-event provider metadata snapshot")
    out=[]
    for source in new["all_frozen_source_study_original_results"]:
        doi=source["published_DOI"]
        native=source["original_native_source_q_and_daylength"]
        if sorted(x["taxon"] for x in native)!=sorted(source["source_registered_native_taxa"]):
            raise ValueError("Registered native taxonomy and opened original source rows diverged")
        for x in native:
            group={}
            for k,p in (("q","original_q_full_source_tag_effect"),
                        ("daylength","original_q_daylength_full_source_tag_effect")):
                z=x[p]
                n=x["original_source_animals_with_diel_q" if k=="q" else "original_source_animals_with_q_daylength"]
                if z is None:
                    group[k]=validated(n,None,None,source["source_status"])
                else:
                    if z["n_individuals"]!=n:
                        raise ValueError("Original source native effect tag support and archive disagree")
                    group[k]=validated(n,z["mean"],z["animal_bootstrap_95_percentile"],
                                       source["source_status"])
            out.append({
                "original_study_doi":doi,"native_taxon":x["taxon"],
                "original_source_stratum":"FOUR_NEW_SOURCE_FROZEN_STUDIES",
                "original_source_status":source["source_status"],**group,
            })
    return out


def percentiles(v):
    v=sorted(v)
    return [v[int(.025*len(v))],v[min(len(v)-1,int(.975*len(v))-1)]]


def study_level(rows,field,seed):
    by=defaultdict(list)
    for x in rows:
        by[x["original_study_doi"]].append(x)
    # Keep source DOI order explicit, regardless of result or group size.
    original=[]
    val=[]
    for doi in sorted(by):
        group=by[doi]
        native=[r[field]["mean"] for r in group if r[field]["mean"] is not None]
        m=mean(native) if native else None
        original.append({
            "source_study_doi":doi,
            "total_native_taxon_rows":len(group),
            "estimable_native_taxa":len(native),
            "equal_available_native_taxon_mean":m,
            "source_mean_status":"UNAVAILABLE" if m is None else
                                 "POSITIVE_POINT" if m>0 else
                                 "NEGATIVE_POINT" if m<0 else "EXACT_ZERO_POINT"
        })
        if native: val.append(m)
    native_all=[r[field]["mean"] for r in rows if r[field]["mean"] is not None]
    if len(val)<3:
        raise ValueError("Too few eligible original field studies")
    rng=random.Random(seed)
    n=len(val)
    draws=[mean(val[rng.randrange(n)] for _ in range(n)) for _ in range(BOOTSTRAP_SAMPLES)]
    loo=[mean(val[:j]+val[j+1:]) for j in range(n)]
    return {
        "original_field_study_DOIs_including_HOLD":len(by),
        "original_studies_with_estimable_native_source_mean":len(val),
        "native_source_taxon_rows_with_estimate":len(native_all),
        "source_DOI_equal_mean":mean(val),
        "native_taxon_row_equal_mean_OVERWEIGHTS_MULTI_TAXON_STUDIES":mean(native_all),
        "source_mean_minus_native_row_mean":mean(val)-mean(native_all),
        "source_DOI_bootstrap95_DESCRIPTIVE_NONPROBABILITY":percentiles(draws),
        "leave_one_original_study_DOI_out_mean_range":[min(loo),max(loo)],
        "original_DOI_source_mean_point_signs":{
            "POSITIVE":sum(x>0 for x in val),"NEGATIVE":sum(x<0 for x in val),
            "EXACT_ZERO":sum(x==0 for x in val)
        },
        "all_original_DOI_source_means_including_unavailable":original,
        "not_global_CI":True
    }


def main_analysis():
    spec,prior,new,validation,metadata=input_contract()
    priorrows=original_eight(prior)
    newrows=original_four(new,validation,metadata)
    allrows=priorrows+newrows
    target=spec["expected_denominators"]
    doi={x["original_study_doi"] for x in allrows}
    taxa={x["native_taxon"] for x in allrows}
    olddoi={x["original_study_doi"] for x in priorrows}
    newdoi={x["original_study_doi"] for x in newrows}
    if (olddoi&newdoi or len(olddoi)!=8 or len(newdoi)!=4 or
        len(doi)!=target["published_studies"] or
        len(allrows)!=target["native_source_taxon_rows"] or
        len(taxa)!=target["distinct_species_labels"] or
        sum(x["q"]["mean"] is not None for x in allrows)!=target["q_supported_rows"] or
        sum(x["daylength"]["mean"] is not None for x in allrows)!=target["photoperiod_supported_rows"]):
        raise ValueError("Original twelve DOI/22 native study×taxon denominator or q/daylight support changed")
    if len({(x["original_study_doi"],x["native_taxon"]) for x in allrows})!=len(allrows):
        raise ValueError("Duplicated native source study taxon falsely increased evidence")
    q=study_level(allrows,"q",RANDOM_SEED)
    day=study_level(allrows,"daylength",RANDOM_SEED+1)
    old_q=study_level(priorrows,"q",RANDOM_SEED+2)
    old_day=study_level(priorrows,"daylength",RANDOM_SEED+3)
    if (q["original_studies_with_estimable_native_source_mean"]!=target["q_supported_studies"] or
        day["original_studies_with_estimable_native_source_mean"]!=target["photoperiod_supported_studies"] or
        len(validation["goose_native_q_and_study_reference_daylength"])!=3 or
        len(validation["all_three_frozen_original_goose_source_native_pairs"])!=3):
        raise ValueError("New original source/three goose all-species calendar controls missing")
    goose=[{
        "native_taxon":x["native_taxon"],
        "shared_120_date_q":x["observed_NIGHT_minus_DAY_GPS_q"],
        "shared_120_date_q_source_tags":x["original_same_calendar_q_tag_count"],
        "shared_120_date_reference_daylength_slope_per_hour":x["q_per_hour_study_reference_daylength_slope"],
        "nominal_original_tag_bootstrap_q":x["q_original_tag_bootstrap95"],
        "nominal_original_tag_bootstrap_daylength":x["daylength_original_tag_bootstrap95"],
    } for x in validation["goose_native_q_and_study_reference_daylength"]]
    return {
        "version":"RSDS_TWELVE_ORIGINAL_STUDY_22_NATIVE_TEMPORAL_MACRO_V1_FIRST_RESULT",
        "date":"2026-10-31",
        "research_stage":"POST_FIRST_SOURCE_OUTCOMES_ORIGINAL_FIELD_STUDY_EQUAL_MACRO_COMPARISON",
        "all_original_studies":len(doi),
        "all_source_native_taxon_rows":len(allrows),
        "distinct_native_species_labels":len(taxa),
        "new_original_studies":len(newdoi),
        "temporal_q_original_study_equal_inference":q,
        "astronomical_daylength_original_study_equal_inference":day,
        "first_eight_study_original_q_DOI_equal_reference":old_q,
        "first_eight_study_original_daylength_DOI_equal_reference":old_day,
        "new_source_native_goose_original_shared_calendar_days":validation["new_geese_original_same_source_calendar_days"],
        "new_source_native_goose_original_year_month_blocks":validation["new_geese_source_year_month_blocks"],
        "one_source_coobserved_three_goose_native_q_and_reference_daylength":goose,
        "one_source_three_goose_pair_q_and_month_reference_slopes":
            validation["all_three_frozen_original_goose_source_native_pairs"],
        "all_original_study_native_taxon_rows_including_HOLD":allrows,
        "native_source_q_nominal_animal_CI_sign_counts":dict(Counter(x["q"]["status"] for x in allrows)),
        "native_source_daylength_nominal_animal_CI_sign_counts":
            dict(Counter(x["daylength"]["status"] for x in allrows)),
        "original_source_git_blob_receipts":{
            z["path"]:z["git_blob_sha"] for z in spec["source_blobs"]
        },
        "external_48_source_negative_shortwave_direction_nonreplication_UNCHANGED":NEGATIVE,
        "original_12_ERA5_raw_JSON_source_HOLDS_UNCHANGED":True,
        "interpretation_ceiling":[
            "Twelve source DOIs are independent publications not randomly selected representative world ecosystems; twenty-two native source-taxon rows are NOT twenty-two independent studies.",
            "The original GPS q is NIGHT−DAY log1p observed straight-line step displacement speed (km/h), not true latent activity time or fractions of daylight used.",
            "Eight-study and twelve-study source-equal comparisons differ in which studies exist; the change is not a time trend or biological change.",
            "Study-equal DOI bootstrap percentiles are resampling sensitivity conditional on purposively selected original field studies, NOT general global population uncertainty.",
            "Native within-study q and astronomical daylength effects may use different original source animal samples and year/day availability.",
            "The goose three-native shared 120-date study-reference daylength is a SEPARATE estimand from original full-source native local daylength slopes, and three goose source native taxa are one field-study unit.",
            "Some source failures reflect label/schema mismatch, GPS sampling cadence or file size rather than ecological absence; all original HOLD rows remain null.",
            "Original native source animal bootstrap nominal intervals are unadjusted for multiple comparisons, no causal photoperiod or evolutionary trait effect identified.",
            "Historic 48-source shortwave signed directional negative external replication and original 12 ERA5 raw JSON version holds remain unchanged.",
            "No private GPS tags, original locations, precise dates, individual effects or raw telemetry included."
        ],
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,default=OUTPUT)
    args=ap.parse_args()
    result=main_analysis()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_TWELVE_ORIGINAL_FIELD_STUDY_TEMPORAL_MACRO",json.dumps({
        "studies":result["all_original_studies"],
        "study_taxon":result["all_source_native_taxon_rows"],
        "source_q":result["temporal_q_original_study_equal_inference"]["source_DOI_equal_mean"],
        "native_q":result["temporal_q_original_study_equal_inference"]["native_taxon_row_equal_mean_OVERWEIGHTS_MULTI_TAXON_STUDIES"],
        "old_8_study_source_q":result["first_eight_study_original_q_DOI_equal_reference"]["source_DOI_equal_mean"],
        "source_q_n":result["temporal_q_original_study_equal_inference"]["original_studies_with_estimable_native_source_mean"],
        "goose_days":result["new_source_native_goose_original_shared_calendar_days"],
        "old_external_nonrep":result["external_48_source_negative_shortwave_direction_nonreplication_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
