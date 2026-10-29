"""Within the SAME published Utah Movebank study, do opposite photoperiod slopes
persist after the two native taxa are restricted to COOBSERVED CALENDAR DATES?

Biological/observational targets are the 2 original taxa (mule deer, bighorn).
The Utah taxa were fixed after seeing Wyoming opposite native signs:
an explicit date-support sensitivity, never a new independent replication.

Private original source provider GPS/animal IDs and exact dates are loaded only
inside an ephemeral CI runner, never exported. The only public source output is
DOI/taxon-level sample coverage and original-tag/date-block descriptive effects.
ERA5 is NOT used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download, daylength
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_utah_diel_photoperiod_20261019 import slope
from scripts.run_rsds_eight_study_gps_cadence_v1 import (
    PAIRS, NEGATIVE, git_blob, per_tag_year_days,
    source_baseline_check, source_reference, verify_frame
)
from scripts.run_rsds_triple_matched_ecological_q_v1 import (
    analyze_taxon as replay_first_strict_three
)

SPEC=HERE/"research/RSDS_UTAH_COOCCURRENT_DATES_PHOTOPERIOD_V1_CONTRACT.json"
PRIOR=HERE/"research/RSDS_TRIPLE_DATE_ECOLOGICAL_Q_DAYLENGTH_INVARIANCE_V1_FIRST_RESULT_20261028.json"
GPS_PRIOR=HERE/"research/RSDS_GPS_THREE_PHASE_EXACT_TAG_DATE_V2_FIRST_RESULT_20261028.json"
TAXA=("Odocoileus hemionus","Ovis canadensis")
SCHEDULES=("original","even","odd")
SEED=20261029
NBOOT=2000
REFERENCE_LATITUDE=43.0


def frozen():
    """Old 8-study source original biological level and triple-delta files SHA-pinned."""
    spec=json.loads(SPEC.read_text(encoding="utf-8"))
    if (spec["version"]!="RSDS_UTAH_COOCCURRENT_DATES_TWO_TAXON_PHOTOPERIOD_V1" and
        spec["version"]!="RSDS_UTAH_COOCCURRENT_DATES_TWO_TAXON_PHOTOPERIOD_V1"):
        raise ValueError("Utah source pre-registered cross-taxon shared-calendar contract changed")
    if (spec["original_study"]["doi"]!="10.5441/001/1.711"
        or spec["original_study"]["study_index"]!=0
        or spec["original_study"]["preselected_two_taxon_comparison"]!=list(TAXA)
        or len(spec["original_study"]["native_taxa_in_full_original_source"])!=4
        or spec["prespecified_support_selection"]["reference_astronomy"].find("43.0")<0):
        raise ValueError("Original Utah source taxa, study DOI or common 43N reference changed")
    prior={}
    for item in spec["immutable_source_results"]:
        path=HERE/item["path"]
        raw=path.read_bytes()
        if git_blob(raw)!=item["git_blob_sha"]:
            raise ValueError("Frozen full 8-study or Wyoming source result modified: "+item["path"])
        prior[item["path"]]=json.loads(raw.decode("utf-8"))
    absolute=prior[str(PRIOR.relative_to(HERE))]
    delta=prior[str(GPS_PRIOR.relative_to(HERE))]
    wyoming=prior["research/RSDS_WYOMING_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json"]
    if (absolute["original_source_studies"]!=8 or
        absolute["original_study_taxon_rows"]!=12 or
        delta["fixed_original_source_studies"]!=8 or
        delta["native_source_study_taxon_rows"]!=12 or
        wyoming["status"]!="STUDY_COMMON_DATE_EVALUATED" or
        wyoming["study_Doi"]!="10.5441/001/1.703" or
        absolute["historical_48_source_signed_shortwave_external_NONREPLICATED_UNCHANGED"]!=NEGATIVE or
        not absolute["historical_original_12_ERA5_source_version_HOLDS_UNCHANGED"]):
        raise ValueError("Original Utah/Wyoming science source denominator/old HOLD changed")
    study=next(x for x in absolute["all_fixed_original_studies_and_native_taxa"]
               if x["original_study_DOI"]=="10.5441/001/1.711")
    triple=next(x for x in delta["all_original_source_studies_and_native_tag_triple_q"]
                if x["original_study_doi"]=="10.5441/001/1.711")
    if ({x["taxon"] for x in study["native_source_results"]}!=
          set(spec["original_study"]["native_taxa_in_full_original_source"]) or
        {x["taxon"] for x in triple["taxa"]}!=
          set(spec["original_study"]["native_taxa_in_full_original_source"])):
        raise ValueError("Original Utah four native taxa source was selectively reduced")
    return spec,study,triple

def as_three_original_calendar(tracks):
    """PRIVATE per-source tag -> year -> astronomy-date -> three q values."""
    original,_=per_tag_year_days(tracks,"ORIGINAL_1H")
    even,_=per_tag_year_days(tracks,PAIRS[0])
    odd,_=per_tag_year_days(tracks,PAIRS[1])
    pool={}
    for tag in sorted(set(original)&set(even)&set(odd)):
        yearmap={}
        for year in sorted(set(original[tag])&set(even[tag])&set(odd[tag])):
            common=set(original[tag][year])&set(even[tag][year])&set(odd[tag][year])
            dates=sorted(common)
            if len(dates)<10:
                continue
            # Use a FIXED astronomy reference; never use q or effect signs to
            # select dates, tags or years. Same x vector for both taxa.
            x=[daylength(REFERENCE_LATITUDE,d) for d in dates]
            if max(x)-min(x)<.5:
                continue
            yearmap[year]={
                d:{"original":original[tag][year][d][0],
                   "even":even[tag][year][d][0],
                   "odd":odd[tag][year][d][0]}
                for d in dates
            }
        if yearmap:
            pool[tag]=yearmap
    return pool


def matched_two_taxon_calendar(pool):
    """Iteratively drop poorly covered tag-years, then intersect source dates.

    Only date support can enter this filter; q RESPONSE VALUES NEVER ENTER.
    """
    if set(pool)!=set(TAXA):
        raise ValueError("Exactly two frozen original native study taxa required")
    valid={sp:{t:{yr:dict(q) for yr,q in years.items()}
               for t,years in pool[sp].items()} for sp in TAXA}
    loops=0
    while True:
        loops+=1
        if loops>100:
            raise ValueError("Shared calendar date support selection did not converge")
        temporal=[]
        for taxon in TAXA:
            temporal.append({(yr,d) for track in valid[taxon].values()
                             for yr,obs in track.items() for d in obs})
        shared=temporal[0]&temporal[1]
        next_valid={}
        for taxon in TAXA:
            next_valid[taxon]={}
            for tag,years in valid[taxon].items():
                chosen={}
                for yr,days in years.items():
                    subset={d:v for d,v in days.items() if (yr,d) in shared}
                    if len(subset)>=10:
                        lights=[daylength(REFERENCE_LATITUDE,d) for d in subset]
                        if max(lights)-min(lights)>=.5:
                            chosen[yr]=subset
                if chosen:
                    next_valid[taxon][tag]=chosen
        if next_valid==valid:
            break
        valid=next_valid
    dates={taxon:{(yr,d) for track in valid[taxon].values()
                  for yr,ds in track.items() for d in ds} for taxon in TAXA}
    common=dates[TAXA[0]]&dates[TAXA[1]]
    if dates[TAXA[0]]!=dates[TAXA[1]]:
        # This should follow from the fixed-point iteration; fail-closed on
        # an accidental source-day support mismatch.
        raise ValueError("Native taxa do not share exactly the study-level calendar date pool")
    return valid,common,loops


def per_tag_effects(tagyears):
    """Equal-year within-tag q and q~common astronomy date daylength."""
    level={}
    gradients={}
    for tag,years in tagyears.items():
        qyear=[]
        lyr=[]
        for yr,records in years.items():
            dates=sorted(records)
            x=[daylength(REFERENCE_LATITUDE,d) for d in dates]
            y={schedule:[records[d][schedule] for d in dates]
               for schedule in SCHEDULES}
            if len(dates)<10 or max(x)-min(x)<.5:
                raise ValueError("Post-qualification original tag/year lost common date support")
            qyear.append({s:mean(y[s]) for s in SCHEDULES})
            derived={s:slope(list(zip(x,y[s]))) for s in SCHEDULES}
            if all(v is not None and math.isfinite(v) for v in derived.values()):
                lyr.append(derived)
        if qyear:
            level[tag]={s:mean(x[s] for x in qyear) for s in SCHEDULES}
        if lyr:
            gradients[tag]={s:mean(x[s] for x in lyr) for s in SCHEDULES}
    return level,gradients


def bootstrap_source(rows,seed):
    """Source animal mean and conditional nominal tag bootstrap across schedules."""
    names=sorted(rows)
    n=len(names)
    out={"original_tag_count":n,"status":"UNAVAILABLE",
         "schedules":{}}
    draws={s:[] for s in SCHEDULES}
    if n>=5:
        rng=random.Random(seed)
        for _ in range(NBOOT):
            sample=[names[rng.randrange(n)] for _ in range(n)]
            for sch in SCHEDULES:
                draws[sch].append(mean(rows[tag][sch] for tag in sample))
    for sch in SCHEDULES:
        values=[rows[tag][sch] for tag in names]
        if n>=5:
            ord_=sorted(draws[sch])
            ci=[ord_[49],ord_[1949]]
            mean_=mean(values)
            sign=("POSITIVE" if ci[0]>0 else
                  "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO")
        else:
            mean_=None
            ci=None
            sign="UNAVAILABLE"
        out["schedules"][sch]={"mean":mean_,"animal_bootstrap95":ci,
                               "nominal_sign":sign,
                               "positive_source_tags":sum(v>0 for v in values),
                               "negative_source_tags":sum(v<0 for v in values)}
    if n>=5:
        out["status"]="ESTIMABLE"
    out["same_tags_and_source_calendar_date_support_for_all_three_schedules"]=True
    return out


def between_taxa_effect(a,b,seed):
    """Resample distinct original tags separately in the two native taxa."""
    names_a=sorted(a)
    names_b=sorted(b)
    n1,n2=len(names_a),len(names_b)
    out={"native_taxon_tag_counts":{TAXA[0]:n1,TAXA[1]:n2},
         "schedules":{}}
    rng=random.Random(seed)
    draws={s:[] for s in SCHEDULES}
    if min(n1,n2)>=5:
        for _ in range(NBOOT):
            aa=[names_a[rng.randrange(n1)] for _ in range(n1)]
            bb=[names_b[rng.randrange(n2)] for _ in range(n2)]
            for sch in SCHEDULES:
                draws[sch].append(mean(a[k][sch] for k in aa)-
                                  mean(b[k][sch] for k in bb))
    for sch in SCHEDULES:
        if min(n1,n2)>=5:
            mean_=mean(a[k][sch] for k in names_a)-mean(b[k][sch] for k in names_b)
            ordered=sorted(draws[sch])
            ci=[ordered[49],ordered[1949]]
            status=("NOMINAL_POSITIVE" if ci[0]>0 else
                    "NOMINAL_NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO")
        else:
            mean_=None
            ci=None
            status="INSUFFICIENT_COMMON_DATE_TAGGED_ANIMALS"
        out["schedules"][sch]={"mule_minus_bighorn_slope":mean_,
                               "descriptive_original_tag_bootstrap95":ci,
                               "status":status}
    return out


def same_calendar_date_difference(pool,shared_dates):
    """PRIVATE date rows for each taxon, no date/time exported."""
    qmean={sp:{} for sp in TAXA}
    for sp in TAXA:
        datewise=defaultdict(lambda:{s:[] for s in SCHEDULES})
        for tagged in pool[sp].values():
            for year,dates in tagged.items():
                for d,v in dates.items():
                    for schedule in SCHEDULES:
                        datewise[(year,d)][schedule].append(v[schedule])
        for key,row in datewise.items():
            if min(len(row[s]) for s in SCHEDULES)<1:
                continue
            qmean[sp][key]={s:mean(row[s]) for s in SCHEDULES}
    if set(qmean[TAXA[0]])!=set(qmean[TAXA[1]]) or set(qmean[TAXA[0]])!=shared_dates:
        raise ValueError("Study daily native taxon common calendar support not equal")
    return {date: {s:qmean[TAXA[0]][date][s]-qmean[TAXA[1]][date][s]
                   for s in SCHEDULES}
            for date in sorted(shared_dates)}


def monthly_block_calendar_contrast(daydiff,seed):
    """Date-paired mule-bighorn q-slope, calendar month-year block resampling.

    Each source date contributes one difference per schedule (NOT number of
    GPS fixes/tags); month-year blocks preserve short-term day correlation.
    """
    rows=sorted(daydiff)
    groups=defaultdict(list)
    for y,d in rows:
        block=(y,d.month)
        groups[block].append((daylength(REFERENCE_LATITUDE,d),daydiff[(y,d)]))
    blocks=sorted(groups)
    out={"coobserved_original_source_dates":len(rows),
         "shared_year_month_blocks":len(blocks),
         "mean_source_qualified_tags_per_species_per_date_NOT_EXPORTED":True,
         "slopes":{}}
    sufficient=(len(rows)>=20 and len(blocks)>=5 and
                max(daylength(REFERENCE_LATITUDE,d) for _,d in rows)-
                min(daylength(REFERENCE_LATITUDE,d) for _,d in rows)>=.5)
    rng=random.Random(seed)
    for sch in SCHEDULES:
        value=slope([(daylength(REFERENCE_LATITUDE,d),daydiff[(y,d)][sch])
                     for y,d in rows]) if sufficient else None
        draws=[]
        if sufficient and value is not None:
            for _ in range(NBOOT):
                sampled=[blocks[rng.randrange(len(blocks))] for _ in blocks]
                flat=[pair for b in sampled for pair in groups[b]]
                yy=slope([(x,v[sch]) for x,v in flat])
                if yy is not None:
                    draws.append(yy)
        if sufficient and value is not None and len(draws)>=1600:
            z=sorted(draws)
            ci=[z[int(.025*len(z))],z[max(0,int(.975*len(z))-1)]]
            state=("NOMINAL_POSITIVE" if ci[0]>0 else
                   "NOMINAL_NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO")
        else:
            ci=None
            state="INSUFFICIENT_CALENDAR_BLOCK_SUPPORT" if not sufficient else "BLOCK_RESAMPLING_UNSTABLE"
        out["slopes"][sch]={"mule_q_minus_bighorn_q_vs_shared_reference_daylength":value,
                            "year_month_block_descriptive95":ci,
                            "status":state}
    return out


def replay_original_triple(tracks,prior_absolute,prior_delta):
    """Use original 3-schedule DELTA input for replay, compare with absolute levels."""
    for sp in TAXA:
        previous=next(x for x in prior_absolute["native_source_results"] if x["taxon"]==sp)
        archived_delta=next(x for x in prior_delta["taxa"] if x["taxon"]==sp)
        got=replay_first_strict_three(sp,tracks[sp],archived_delta)
        for newer,archived in (
            ("source_q_sign_recurrence","source_q_three_schedule_nominal_sign_stability"),
            ("source_daylength_slope_sign_recurrence","source_daylength_three_schedule_nominal_sign_stability")
        ):
            if got[newer]!=previous[archived]:
                raise ValueError("Old Wyoming original GPS q/photoperiod sign changed")
        for dst,old_key in (("q_original_tag_means_on_exact_same_dates",
                              "original_same_date_q_source_estimates"),
                            ("q_daylength_per_hour_on_exact_same_dates",
                              "original_same_date_daylength_slope_estimates")):
            for ax,new in got[dst].items():
                old=previous[old_key][ax]
                if (new["n_source_tags"]!=old["n_source_tags"]
                        or new["status"]!=old["status"] or
                        (new["group_original_tag_mean"] is None)!=
                        (old["group_original_tag_mean"] is None)):
                    raise ValueError("Old original tagged animal per taxon q/daylight source support changed")
                if new["group_original_tag_mean"] is not None and abs(
                        new["group_original_tag_mean"]-old["group_original_tag_mean"])>1e-12:
                    raise ValueError("Old original q/daylength coefficient changed")
    return True


def run(out):
    spec,previous,previous_delta=frozen()
    old,new,_,_=verify_frame()
    region="Utah"
    src,taxa,reference,_=source_reference(region,old,new)
    if (src["doi"]!=spec["original_study"]["doi"] or set(taxa)!=set(spec["original_study"]["native_taxa_in_full_original_source"])):
        raise ValueError("Utah original published source DOI/two species identity changed")
    result={
      "version":"RSDS_UTAH_COOCCURRENT_CALENDAR_DAY_PHOTOPERIOD_V1_ORIGINAL_RESULT",
      "date":"2026-10-29",
      "study_Doi":src["doi"],
      "native_taxa":list(TAXA),
      "status":"ORIGINAL_STUDY_HOLD",
      "original_provider_GPS_MD5_verified":False,
      "both_native_taxa_original_all_date_and_original_triple_source_estimates_replayed":False,
      "shared_support":None,
      "source_taxon_results":[],
      "same_study_between_native_taxon":None,
      "calendar_date_paired_source_difference":None,
      "historical_12_ERA5_raw_JSON_HOLDS_UNCHANGED":True,
      "historical_48_source_signed_radiation_external_NONREPLICATED_UNCHANGED":NEGATIVE
    }
    out.mkdir(parents=True,exist_ok=True)
    private=out/"private_utah_same_study_coobserved_calendar_gps.zip"
    stage="ORIGINAL_PROVIDER_METADATA"
    try:
        original=[x for x in original_files(src["package_uuid"])
                  if x["bitstream_uuid"]==src["original_bitstream_uuid"]]
        if len(original)!=1 or original[0]["bytes"]!=src["original_bytes"] or \
           original[0]["checksum"]["value"].lower()!=src["original_md5"].lower():
            raise ValueError("Original Utah original public source UUID/file bytes/MD5 changed")
        stage="ORIGINAL_PROVIDER_FILE_MD5"
        verify=checksum_and_download(original[0],private)
        if not verify["provider_checksum_verified"]:
            raise ValueError("Original Utah source provider MD5 did not match")
        result["original_provider_GPS_MD5_verified"]=True
        stage="ORIGINAL_SOURCE_GPS_PARSING"
        tracks,counts,unknown=extract_time_only_original_source(private,taxa)
        stage="ORIGINAL_NATIVE_FULL_DATE_Q_DAYLENGTH_REPLAY"
        source_baseline_check(tracks,reference)
        stage="PREVIOUS_TRIPLE_DATE_Q_DAYLENGTH_REPLAY"
        replay_original_triple(tracks,previous,previous_delta)
        result["both_native_taxa_original_all_date_and_original_triple_source_estimates_replayed"]=True
        print("RSDS_UTAH_SHARED_DATE_PREVIOUS_NATIVE_Q_REPLAY",len(taxa),flush=True)

        stage="STUDY_SHARED_CALENDAR_SUPPORT_AND_NEW_ORIGINAL_Q"
        private_pool={sp:as_three_original_calendar(tracks[sp]) for sp in TAXA}
        shared_pool,shared_dates,iterations=matched_two_taxon_calendar(private_pool)
        initial={sp:len(private_pool[sp]) for sp in TAXA}
        final={sp:len(shared_pool[sp]) for sp in TAXA}
        common_seasons={d.month for y,d in shared_dates}
        years={y for y,d in shared_dates}
        out_support={
            "original_tags_with_3_schedule_10_day_local_eligible_season":
                {sp:initial[sp] for sp in TAXA},
            "original_tags_with_10_days_on_both_taxa_shared_calendar":
                {sp:final[sp] for sp in TAXA},
            "original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules":
                len(shared_dates),
            "distinct_original_calendar_years_with_coobserved_dates":len(years),
            "distinct_calendar_month_numbers_in_coobserved_dates":len(common_seasons),
            "date_support_fixed_point_iterations":iterations,
            "common_daylight_standard_reference_latitude_degrees":REFERENCE_LATITUDE,
            "no_original_calendar_years_or_dates_or_tag_IDs_exported":True,
        }
        result["shared_support"]=out_support
        level={}
        slope_by={}
        for k,sp in enumerate(TAXA):
            q,s=per_tag_effects(shared_pool[sp])
            level[sp]=q
            slope_by[sp]=s
            item={
                "taxon":sp,
                "original_study_coobserved_date_tagged_animals_q":len(q),
                "original_study_coobserved_date_tagged_animals_daylength":len(s),
                "common_date_temporal_q_by_schedule":bootstrap_source(q,SEED+100*k),
                "common_date_photoperiod_q_slope_by_schedule":
                    bootstrap_source(s,SEED+200*k),
            }
            result["source_taxon_results"].append(item)
        result["same_study_between_native_taxon"]={
            "q_mule_minus_bighorn_by_schedule":
                between_taxa_effect(level[TAXA[0]],level[TAXA[1]],SEED+300),
            "q_daylength_slope_mule_minus_bighorn_by_schedule":
                between_taxa_effect(slope_by[TAXA[0]],slope_by[TAXA[1]],SEED+500),
        }
        if shared_dates:
            datewise=same_calendar_date_difference(shared_pool,shared_dates)
            result["calendar_date_paired_source_difference"] = \
                monthly_block_calendar_contrast(datewise,SEED+700)
        else:
            result["calendar_date_paired_source_difference"]={
                "coobserved_original_source_dates":0,"shared_year_month_blocks":0,
                "slopes":{sch:{"status":"INSUFFICIENT_CALENDAR_BLOCK_SUPPORT",
                              "mule_q_minus_bighorn_q_vs_shared_reference_daylength":None,
                              "year_month_block_descriptive95":None}
                          for sch in SCHEDULES}
            }
        result["status"]="STUDY_COMMON_DATE_EVALUATED" if min(final.values())>=5 else \
                         "INSUFFICIENT_COOCCURRENT_SOURCE_DATE_SUPPORT"
        for sp in TAXA:
            tracks[sp].clear()
    except Exception as e:
        result["status"]="ORIGINAL_UTAH_SOURCE_OR_SHARED_DATE_HOLD"
        result["hold_stage"]=stage
        result["hold_error_type"]=type(e).__name__
        print("RSDS_UTAH_SHARED_DATE_SOURCE_HOLD",json.dumps({
            "hold_stage":stage,"exception_type":type(e).__name__,
            "original_md5_passed":result["original_provider_GPS_MD5_verified"]
        },sort_keys=True),flush=True)
    finally:
        private.unlink(missing_ok=True)
    (out/"utah_shared_calendar_two_taxa_result_v1.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8"
    )
    print("RSDS_UTAH_TWO_TAXON_SHARED_DATES",json.dumps({
        "status":result["status"],
        "original_MD5":result["original_provider_GPS_MD5_verified"],
        "native_q_replayed":result["both_native_taxa_original_all_date_and_original_triple_source_estimates_replayed"],
        "source_taxa":len(result["source_taxon_results"]),
        "calendar_common_dates":(result["shared_support"] or {}).get(
            "original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"),
        "common_tag_n":(result["shared_support"] or {}).get(
            "original_tags_with_10_days_on_both_taxa_shared_calendar"),
        "old_holds_unchanged":True
    },sort_keys=True),flush=True)
    return result


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,required=True)
    a=cli.parse_args()
    run(a.out)


if __name__=="__main__":
    main()
