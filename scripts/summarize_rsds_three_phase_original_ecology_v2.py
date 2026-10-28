"""Independent audit of frozen original 8-study exact triple tag/date GPS q effects.

No new animal GPS is touched: original eight-study original-source source
first results and previous 2-phase findings are verified by their original
Git blob SHA. Summarizes *observation-clock-phase dependence*, not a new
solar/photoperiod mechanism, a representative world ecology effect or a
predeclared phase-response hypothesis test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE="research/RSDS_GPS_THREE_PHASE_EXACT_TAG_DATE_V2_FIRST_RESULT_20261028.json"
OLDER="research/RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_FIRST_RESULT_20261028.json"
PIN={
    SOURCE:"da679f8ca8d973de7da55463e5643037df6a0579",
    OLDER:"8b93fe361eccbdc13647e465ddfaa24ad10c7c98",
}
NEG="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
AXES=("q_even_minus_original","q_odd_minus_original","q_even_minus_odd",
      "daylength_even_minus_original","daylength_odd_minus_original",
      "daylength_even_minus_odd")
DEST="research/RSDS_GPS_TRIPLE_DATE_CADENCE_ECOLOGICAL_INTERPRETATION_V2_20261028.json"


def pin(root,path):
    raw=(root/path).read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
    if sha!=PIN[path]:
        raise ValueError("A frozen original GPS source group first result changed: "+path)
    return json.loads(raw.decode("utf-8"))


def run(root=ROOT):
    root=Path(root)
    current=pin(root,SOURCE)
    before=pin(root,OLDER)
    if (current["fixed_original_source_studies"]!=8 or
            current["native_source_study_taxon_rows"]!=12 or
            before["original_field_studies"]!=8 or
            current["old_48_original_source_negative_shortwave_direction_external_replication_UNCHANGED"]!=NEG
            or not current["old_12_ERA5_source_RAW_SHA_HOLDS_UNCHANGED"]):
        raise ValueError("Source original eight-study native taxon/negative outcome not intact")
    rows=[]
    summary={axis:Counter() for axis in AXES}
    study_counts=Counter()
    for s in current["all_original_source_studies_and_native_tag_triple_q"]:
        status=s["source_status"]
        study_counts[status]+=1
        for t in s["taxa"]:
            q=t.get("new_strict_triple_original_tag_effects")
            effects={}
            if q is not None:
                if set(q)!=set(AXES):
                    raise ValueError("One original GPS triple-match source endpoint missing")
                for ax in AXES:
                    v=q[ax]
                    sign=v["status"]
                    if sign not in ("NOMINAL_POSITIVE","NOMINAL_NEGATIVE","INCLUDES_ZERO","UNAVAILABLE"):
                        raise ValueError("Unknown original study native-tag bootstrap sign")
                    summary[ax][sign]+=1
                    ci=v["descriptive_source_tag_bootstrap95"]
                    n=v["source_original_same_tags_and_triple_dates"]
                    x=v["group_source_mean_delta"]
                    if sign=="UNAVAILABLE":
                        if n>=5 or x is not None or ci is not None:
                            raise ValueError("Unqualified original source produced numeric confidence")
                    else:
                        if n<5 or len(ci)!=2 or not ci[0]<=x<=ci[1]:
                            raise ValueError("Unreliable original tagged animal bootstrap interval")
                    effects[ax]={
                        "same_original_tags_and_calendar_dates_with_estimable_axis":n,
                        "original_tag_group_mean":x,
                        "source_original_tag_bootstrap95":ci,
                        "nominal_bootstrap_interval_sign":sign
                    }
                a,b,c=(effects[x]["original_tag_group_mean"] for x in
                         ("q_even_minus_original","q_odd_minus_original","q_even_minus_odd"))
                if a is not None and abs((a-b)-c)>1e-12:
                    raise ValueError("Triply matched original tag/date q contrast is not algebraically closed")
            rows.append({
                "study_doi":s["original_study_doi"],
                "native_taxon":t["taxon"],
                "first_source_status":status,
                "triple_original_tag_support":t.get("triple_original_animal_day_support"),
                "strict_three_phase_same_tag_date_effects":effects if q is not None else None,
            })
    if len(rows)!=12 or len(study_counts)!=2 or study_counts.get("TRIPLE_SOURCE_SUPPORT_EVALUATED")!=7:
        raise ValueError("Original eight studies or native cohorts incomplete")
    supported=[x for x in rows if x["strict_three_phase_same_tag_date_effects"] is not None
               and x["strict_three_phase_same_tag_date_effects"]["q_even_minus_odd"][
                   "nominal_bootstrap_interval_sign"]!="UNAVAILABLE"]
    positive=[{"study_doi":x["study_doi"],"taxon":x["native_taxon"],
               **x["strict_three_phase_same_tag_date_effects"]["q_even_minus_odd"]}
              for x in supported if x["strict_three_phase_same_tag_date_effects"][
                  "q_even_minus_odd"]["nominal_bootstrap_interval_sign"]=="NOMINAL_POSITIVE"]
    if len(supported)!=8 or len(positive)!=3 or {x["study_doi"] for x in positive}!={
        "10.5441/001/1.703","10.5441/001/1.711"}:
        raise ValueError("Original triple-match supported source groups or timing inference shifted")
    bytaxon={(x["study_doi"],x["native_taxon"]):x for x in rows}
    wyoming={}
    for taxon in ("Odocoileus hemionus","Ovis canadensis"):
        x=bytaxon[("10.5441/001/1.703",taxon)]
        s=x["strict_three_phase_same_tag_date_effects"]
        if not (s["q_even_minus_original"]["nominal_bootstrap_interval_sign"]=="NOMINAL_POSITIVE"
                and s["q_odd_minus_original"]["nominal_bootstrap_interval_sign"]=="NOMINAL_NEGATIVE"
                and s["q_even_minus_odd"]["nominal_bootstrap_interval_sign"]=="NOMINAL_POSITIVE"):
            raise ValueError("Full/EVEN/ODD original same-date Wyoming source phase signs changed")
        wyoming[taxon]={
            "original_triply_matched_tag_count":
                s["q_even_minus_odd"]["same_original_tags_and_calendar_dates_with_estimable_axis"],
            "q_even_minus_original":s["q_even_minus_original"],
            "q_odd_minus_original":s["q_odd_minus_original"],
            "q_even_minus_odd":s["q_even_minus_odd"]
        }
    return {
        "version":"RSDS_EIGHT_STUDY_ORIGINAL_GPS_TRIPLE_TAG_DATE_ECOLOGICAL_INTERPRETATION_V2",
        "date":"2026-10-28",
        "status":"SOURCE_ORIGINAL_TAG_DATE_MATCHED_OBSERVATION_PROCESS_PHASE_BIAS",
        "original_field_studies":8,
        "original_study_taxon_rows":12,
        "qualifying_same_tag_original_three_phase_source_taxon_rows":8,
        "qualified_original_source_taxa_with_positive_even_minus_odd_q":{
            "n":3,"source_taxa":positive
        },
        "original_Wyoming_two_native_taxa_exact_three_phase_recurrence":wyoming,
        "source_state_counts":dict(study_counts),
        "source_axis_nominal_interval_counts":{k:dict(v) for k,v in summary.items()},
        "source_original_first_result_git_blobs":PIN,
        "all_eight_source_study_native_taxa_group_effects":rows,
        "previous_12_ERA5_source_json_version_HOLDS_UNCHANGED":True,
        "previous_48_source_signed_sunlight_external_NONREPLICATED_UNCHANGED":NEG,
        "interpretation_ceiling":[
            "The exact same original tagged animals, calendar years and astronomy-qualified dates are used in the ORIGINAL/EVEN/ODD contrast. This rules out different admissible animal/date sets as the sole explanation of the observed Wyoming source sign reversal.",
            "It does not identify biological responses to light, sun azimuth, heat or UTC clock parity. The GPS travel segment length and acquisition phase define different OBSERVATION estimands.",
            "Three-group supported positive even-minus-odd nominal q differences come from only TWO independent original field study DOIs; they do not represent globally or phylogenetically general animal behavior.",
            "Native source animals only conditionally sampled; 2000 paired animal-bootstrap 95% intervals are post-result descriptive and uncorrected for multiple axes and taxa.",
            "Unavailable triple-match original groups remain UNAVAILABLE, including Turkish turkey and Thai fishing cat temporal/date support and old great egret 200MB source HOLD.",
            "Original source 'hourly' GPS reflects at most one location within UTC-hour bins, not a guarantee of 1h actual between-fix intervals. Utah native fix timestamps overwhelmingly favor odd UTC hours.",
            "The q response estimates log1p NIGHT−DAY straight-line observed GPS step speed, not latent duration/activity or identical true movement paths at different fix intervals.",
            "Source study DOI and native taxon rows are NOT interchangeable independent world ecosystem units.",
            "Old ERA5 source response 12/12 serialized JSON SHA version holds and previous 48-source radiation-direction external nonreplication are unchanged."
        ],
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,default=ROOT/DEST)
    a=cli.parse_args()
    d=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_TRIPLE_TAG_DATE_PHASE_ECOLOGY",json.dumps({
        "studies":d["original_field_studies"],
        "source_taxa":d["original_study_taxon_rows"],
        "with_triple_n_ge_5":d["qualifying_same_tag_original_three_phase_source_taxon_rows"],
        "positive_even_minus_odd":d["qualified_original_source_taxa_with_positive_even_minus_odd_q"]["n"],
        "wyoming":{taxon:{"n":q["original_triply_matched_tag_count"],
                "even_minus_odd":q["q_even_minus_odd"]["original_tag_group_mean"]}
                for taxon,q in d["original_Wyoming_two_native_taxa_exact_three_phase_recurrence"].items()},
        "historical_ERA5_HOLDS":d["previous_12_ERA5_source_json_version_HOLDS_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
