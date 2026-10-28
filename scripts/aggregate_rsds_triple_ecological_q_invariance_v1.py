"""Original eight-study temporal ecology sign stability across 3 matched GPS clocks.

All eight original Movebank source DOIs and 12 study-native taxa retained,
including original great egret file-size HOLD, original annual/day support and
new unavailable source states. Scientific unit = ORIGINAL source study, never
i.i.d. GPS fix or study×taxon pseudoreplication. Only descriptive original
tag-level bootstrap bands conditional on exactly triple matched original days.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean

from scripts.run_rsds_triple_matched_ecological_q_v1 import (
    HERE,freeze,Q_KEYS,SLOPE_KEYS,SCHEDULES,NEGATIVE,classify_stability
)

OUT=HERE/"research/RSDS_TRIPLE_DATE_ECOLOGICAL_Q_DAYLENGTH_INVARIANCE_V1_FIRST_RESULT_20261028.json"
ALLOWED=("POSITIVE","NEGATIVE","INCLUDES_ZERO","UNAVAILABLE")
ALL_K=Q_KEYS+SLOPE_KEYS


def consistency(source):
    spec,previous=freeze()
    if len(source)!=8 or set(source)!=set(range(8)):
        raise ValueError("Original eight source study DOI denominator incomplete")
    rows=[]
    studies=[]
    signs={axis:Counter() for axis in ALL_K}
    stability={"q":Counter(),"daylength":Counter()}
    state=Counter()
    original_study_q={}
    for i,old in enumerate(previous["all_original_source_studies_and_native_tag_triple_q"]):
        s=source[i]
        fixed={x["taxon"] for x in old["taxa"]}
        new=s.get("native_taxon_results",[])
        if (s.get("version")!="RSDS_8_STUDY_TRIPLE_DATE_ECOLOGICAL_Q_AND_DAYLENGTH_V1_SOURCE"
                or s.get("index")!=i
                or s.get("original_study_doi")!=old["original_study_doi"]
                or set(s.get("native_taxa",[]))!=fixed
                or len(new)!=len(fixed) or {t["taxon"] for t in new}!=fixed
                or not s.get("previous_12_ERA5_raw_version_HOLDS_UNCHANGED")
                or s.get("previous_48_source_signed_shortwave_external_NONREPLICATED_UNCHANGED")!=NEGATIVE):
            raise ValueError("Original source/provider/native study taxonomy or historic outcome was altered")
        if old["source_status"]=="PREDECLARED_FILE_METADATA_HOLD":
            if s["status"]!="PREDECLARED_FILE_METADATA_HOLD":
                raise ValueError("Original source egret >200MB file HOLD cannot change")
        elif s["status"]=="SOURCE_TRIPLE_ECOLOGICAL_Q_COMPLETED":
            if not s["original_provider_bitstream_MD5_confirmed"] or not s["original_source_q_and_daylength_replayed_before_new_q_levels"]:
                raise ValueError("Original GPS biological q/daylight levels were computed before source MD5/q replay")
        state[s["status"]]+=1
        entry={
            "original_study_DOI":s["original_study_doi"],
            "original_study_index":i,
            "study_native_taxa_count":len(new),
            "status":s["status"],
            "old_source_q_daylength_baseline_reproduced":s["original_source_q_and_daylength_replayed_before_new_q_levels"],
            "native_source_results":[]
        }
        q_native=[]
        for group in new:
            old_group=next(x for x in old["taxa"] if x["taxon"]==group["taxon"])
            oldsupport=old_group.get("triple_original_animal_day_support")
            q=group.get("q_original_tag_means_on_exact_same_dates")
            grad=group.get("q_daylength_per_hour_on_exact_same_dates")
            n=group.get("original_three_schedule_source_support")
            statecode=group.get("source_status","SOURCE_HOLD")
            if s["status"]=="SOURCE_TRIPLE_ECOLOGICAL_Q_COMPLETED":
                if not isinstance(q,dict) or not isinstance(grad,dict) or set(q)!=set(Q_KEYS) or set(grad)!=set(SLOPE_KEYS):
                    raise ValueError("Original source complete but three schedule absolute q/daylength levels missing")
                if n is None or oldsupport is None:
                    raise ValueError("Original source strict triple support missing")
                if (n["original_tags_with_at_least_ten_triple_qualified_days"]!=oldsupport["source_original_tags_with_at_least_ten_triple_days"] or
                    n["matched_original_tagged_animal_calendar_days"]!=oldsupport["triply_qualified_original_animal_days"] or
                    n["original_tags_with_common_photoperiod_slope_support"]!=oldsupport["source_original_tags_with_daylength_slope_triple_support"]):
                    raise ValueError("Strict original tagged animal/date sample changed")
                if statecode=="THREE_SCHEDULE_ECOLOGICAL_Q_ESTIMABLE" and n["original_tags_with_at_least_ten_triple_qualified_days"]<5:
                    raise ValueError("Original source biological q signed claim from fewer than five tagged animals")
                if statecode=="INSUFFICIENT_ORIGINAL_TRIPLE_TAG_DATE_SUPPORT" and n["original_tags_with_at_least_ten_triple_qualified_days"]>=5:
                    raise ValueError("Original valid source biological q was silently dropped")
                for k in ALL_K:
                    obj=(q if k in Q_KEYS else grad)[k]
                    if obj["status"] not in ALLOWED:
                        raise ValueError("Unknown source temporal ecological bootstrap sign")
                    nc=obj["n_source_tags"]
                    expected_n=(n["original_tags_with_at_least_ten_triple_qualified_days"]
                                if k in Q_KEYS else
                                n["original_tags_with_common_photoperiod_slope_support"])
                    if nc!=expected_n:
                        raise ValueError("Original native source tag number differs for q or photoperiod")
                    if nc<5:
                        if obj["status"]!="UNAVAILABLE" or obj["group_original_tag_mean"] is not None or obj["nominal_source_tag_bootstrap95"] is not None:
                            raise ValueError("Source insufficient original tags became fake biological zero")
                    else:
                        ci=obj["nominal_source_tag_bootstrap95"]
                        m=obj["group_original_tag_mean"]
                        if not isinstance(ci,list) or len(ci)!=2 or not all(math.isfinite(x) for x in (m,ci[0],ci[1])) or not ci[0]<=m<=ci[1]:
                            raise ValueError("Original tagged animal q or photoperiod source bootstrap invalid")
                        cls=("POSITIVE" if ci[0]>0 else "NEGATIVE" if ci[1]<0 else "INCLUDES_ZERO")
                        if obj["status"]!=cls:
                            raise ValueError("Original native source sign not consistent with its bootstrap95")
                    signs[k][obj["status"]]+=1
                qstable=group["source_q_sign_recurrence"]
                gstable=group["source_daylength_slope_sign_recurrence"]
                qset=[q[k]["status"] for k in Q_KEYS]
                gset=[grad[k]["status"] for k in SLOPE_KEYS]
                if qstable!=classify_stability(q,Q_KEYS):
                    raise ValueError("Original source full/even/odd q sign stability label disagrees with three actual intervals")
                if gstable!=classify_stability(grad,SLOPE_KEYS):
                    raise ValueError("Original source full/even/odd daylength sign stability label disagrees with three actual intervals")
                stability["q"][qstable]+=1
                stability["daylength"][gstable]+=1
                if all(q[k]["group_original_tag_mean"] is not None for k in Q_KEYS):
                    q_native.append({k:q[k]["group_original_tag_mean"] for k in Q_KEYS})
            entry["native_source_results"].append({
                "taxon":group["taxon"],
                "status":statecode,
                "strict_identical_original_tag_date_support":n,
                "original_same_date_q_source_estimates":q,
                "original_same_date_daylength_slope_estimates":grad,
                "source_q_three_schedule_nominal_sign_stability":group.get("source_q_sign_recurrence"),
                "source_daylength_three_schedule_nominal_sign_stability":group.get("source_daylength_slope_sign_recurrence")
            })
            rows.append({"source_study_DOI":s["original_study_doi"],
                         "taxon":group["taxon"],
                         "source_status":statecode,
                         "source_q_three_schedule_nominal_sign_stability":group.get("source_q_sign_recurrence"),
                         "source_daylength_three_schedule_nominal_sign_stability":group.get("source_daylength_slope_sign_recurrence")})
        original_study_q[s["original_study_doi"]]={
            "native_original_taxa_with_5_all_schedule_q_tags":len(q_native),
            "source_study_equal_native_taxon_means_on_triple_dates":{
                k:mean(x[k] for x in q_native) if q_native else None
                for k in Q_KEYS
            },
            "one_study_DOI_not_multiple_world_replicates":True,
        }
        studies.append(entry)
    if len(studies)!=8 or len(rows)!=12:
        raise ValueError("Eight published original GPS source DOI / twelve native taxa lost")
    return {
        "version":"RSDS_8_STUDY_TRIPLE_DATE_TEMPORAL_ECOLOGY_SIGN_SUPPORT_V1_FIRST",
        "date":"2026-10-28",
        "research_stage":"POST_PREVIOUS_ECOLOGICAL_RESULTS_EXACT_ORIGINAL_TAG_DATE_CADENCE_ROBUSTNESS",
        "original_source_studies":8,
        "original_study_taxon_rows":12,
        "study_status_counts":dict(state),
        "q_source_tag_bootstrap_interval_signs_by_fixed_schedule":{
            k:dict(signs[k]) for k in Q_KEYS
        },
        "q_daylength_source_tag_bootstrap_interval_signs_by_fixed_schedule":{
            k:dict(signs[k]) for k in SLOPE_KEYS
        },
        "native_taxa_q_three_schedule_nominal_sign_invariance":dict(stability["q"]),
        "native_taxa_daylength_three_schedule_nominal_sign_invariance":dict(stability["daylength"]),
        "original_study_DOI_equal_source_q_levels":original_study_q,
        "all_fixed_original_studies_and_native_taxa":studies,
        "historical_original_12_ERA5_source_version_HOLDS_UNCHANGED":True,
        "historical_48_source_signed_shortwave_external_NONREPLICATED_UNCHANGED":NEGATIVE,
        "science_ceiling":[
            "These are direct absolute original within-tag GPS q and within-tag-year astronomical daylength slopes on EXACT matched animal/date support in all three schedules, not the prior separately matched original first study group means.",
            "Source q is NIGHT−DAY log1p observed straight-line GPS displacement speed and not true latent biological activity duration. Source q from original/hourbin, even UTC and odd UTC use distinct measurement chords and sometimes vastly different fix availability.",
            "A stable nominal sign under all three schedules is a source sample/schedule conditional descriptive result; simultaneous confidence and world-population generality are NOT established, many taxon/axis comparisons were made and the analysis was posthoc.",
            "Some source study×taxon q and daylength effects have INCLUDES_ZERO intervals; this is not zero animal biology or biological equivalence.",
            "One original source DOI is the study unit; four Utah taxa and two Wyoming taxa are correlated native source cohorts, not independent field-study replications.",
            "Do not translate an apparent UTC sampling phase effect to photoperiod manipulation, animal clock entrainment, solar navigation, heat response or causal radiation.",
            "The old 12 source ERA5 provider serialized raw JSON SHA holds and 48-source independent negative shortwave direction external nonreplication remain unchanged.",
            "No individual original GPS tags, exact original position/date/UTC observation time or per-individual q values published."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    sources={}
    for i in range(8):
        src=args.input/("source_"+str(i)+"_triple_q_ecology.json")
        if not src.is_file():
            raise ValueError("Missing frozen original eight-study native GPS source: "+str(i))
        sources[i]=json.loads(src.read_text(encoding="utf-8"))
    result=consistency(sources)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("RSDS_TRIPLE_SAME_DATE_TIME_NICHE_ECOLOGY",json.dumps({
        "original_studies":result["original_source_studies"],
        "native_taxa":result["original_study_taxon_rows"],
        "study_states":result["study_status_counts"],
        "q_stability":result["native_taxa_q_three_schedule_nominal_sign_invariance"],
        "photoperiod_stability":result["native_taxa_daylength_three_schedule_nominal_sign_invariance"],
        "old_holds":result["historical_original_12_ERA5_source_version_HOLDS_UNCHANGED"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
