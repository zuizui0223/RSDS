"""Original GPS ecological macro atlas: 8 archived studies / 12 native source taxa.

Integrates the previously checked Utah/Wyoming study-taxon time×Sun-source
estimands with the separately pre-outcome frozen six new Movebank studies.
All source input Git blob identities are pinned. No raw telemetry used;
all missing/HOLDs retained. Descriptive source conditional, not causal/global.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
"research/RSDS_CROSS_AXIS_GENERALITY_FIRST_RESULT_20261026.json":
    "951d27f4b772507b0c26af3118186639e1591bec",
"research/RSDS_NEW_SIX_ORIGINAL_GPS_ECOLOGY_FIRST_RESULT_20261026.json":
    "e2759ddb1120cbcd067dd5189775d685aea65a32",
"research/MOVEBANK_UTAH_DIEL_PHOTOPERIOD_VERIFIED_RESULT_V0_20261019.json":
    "5870c6920d83258c5b62b36dc7d043802c494fef",
"research/wyoming_two_ungulates_20261022/wyoming_original_result.json":
    "b9428dd18f2556e1bed245894546d6399da045dc",
"research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_FIRST_RESULT_20261026.json":
    "1ec684797f0b7aae0ec441163e1bbeec48d3bd61",
}
NEGATIVE = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def frozen_load(root, path):
    raw=(root/path).read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
    if sha != SOURCES[path]:
        raise ValueError("One original scientific first result changed unexpectedly: "+path)
    return json.loads(raw.decode("utf-8"))


def observed(value, interval, n, point_only=False):
    if value is None:
        if interval is not None:
            raise ValueError("Undefined ecological source q or direction must not have a confidence interval")
        return {"status":"UNAVAILABLE","mean":None,"source_animals":n,"animal_bootstrap95":None}
    if not isinstance(value,(float,int)) or not math.isfinite(value):
        raise ValueError("Nonfinite native original field-study estimate")
    if point_only:
        if interval is not None:
            raise ValueError("Original point-only joint direction cannot manufacture a CI")
        return {"status":"POINT_ONLY_NO_ORIGINAL_CI","mean":value,"source_animals":n,
                "animal_bootstrap95":None}
    if interval is None or len(interval)!=2 or n<5:
        raise ValueError("Supported q/daylength or raw solar effect lacks original tag uncertainty")
    low,high=interval
    if not all(isinstance(v,(int,float)) and math.isfinite(v) for v in interval):
        raise ValueError("Nonfinite original bootstrap interval")
    if not low<=value<=high:
        raise ValueError("Original point estimate beyond its declared tag bootstrap limits")
    status=("POSITIVE_NOMINAL_SOURCE_TAG_INTERVAL" if low>0 else
            "NEGATIVE_NOMINAL_SOURCE_TAG_INTERVAL" if high<0 else "INCLUDES_ZERO")
    return {"status":status,"mean":value,"source_animals":n,"animal_bootstrap95":interval}


def normalize_new_source_axis(axis):
    """Unify labels without reclassifying any first-original source result.

    Previous and newly frozen studies used different *names* for the same
    nominal original-tag bootstrap interval signs. Recheck all point estimates
    and intervals, and fail rather than silently disagree with first results.
    """
    names = {
        "POSITIVE_NOMINAL_BOOTSTRAP": "POSITIVE_NOMINAL_SOURCE_TAG_INTERVAL",
        "NEGATIVE_NOMINAL_BOOTSTRAP": "NEGATIVE_NOMINAL_SOURCE_TAG_INTERVAL",
        "INCLUDES_ZERO": "INCLUDES_ZERO",
        "UNAVAILABLE": "UNAVAILABLE",
    }
    if axis["status"] not in names:
        raise ValueError("New source original bootstrapped axis status unknown")
    rec = observed(axis["mean"], axis["bootstrap95"], axis["n_animals"])
    if rec["status"] != names[axis["status"]]:
        raise ValueError("New original source outcome inconsistent with its frozen sign interval")
    return rec


def macro(root=ROOT):
    root=Path(root).resolve()
    old,new,utah,wyoming,geo=(
        frozen_load(root,p) for p in SOURCES
    )
    if (len(old["joint_time_direction_source_taxon_atlas"])!=6 or
            new["frozen_source_original_study_count"]!=6 or
            new["number_of_studies_by_source_status"].get("ORIGINAL_METADATA_HOLD")!=1 or
            new["previous_external_48_study_solar_radiation_nonreplication_UNCHANGED"]!=NEGATIVE or
            geo["external_48_original_source_shortwave_nonreplication_UNCHANGED"]!=NEGATIVE):
        raise ValueError("Previous original study records or source HOLD denominator were changed")
    if old["sunlight_weather_temporal_observation_overlap"]["preregistered_heldout_decision_UNCHANGED"]!=NEGATIVE:
        raise ValueError("Original prior frozen external radiation nonreplication modified")

    utah_by={z["taxon"]:z for z in utah["species"]}
    wyoming_by={z["taxon"]:z for z in wyoming["source_taxa"]}
    paired_old={(x["region"],x["taxon"]):x for x in geo["all_original_source_taxon_results"]}
    rows=[]
    for x in old["joint_time_direction_source_taxon_atlas"]:
        region, taxon=x["source_region"],x["taxon"]
        if region=="Utah":
            r=utah_by[taxon]
            if r["eligible_diel_animals"]!=x["q_n_animals"] or abs(r["diel_mean"]-x["night_day_q"])>1e-12:
                raise ValueError("Utah original q and source/native taxon identity mismatch")
            ci=r["diel_bootstrap_95"]
        elif region=="Wyoming":
            r=wyoming_by[taxon]["day_night_speed_contrast"]
            if r["n_individuals"]!=x["q_n_animals"] or abs(r["mean"]-x["night_day_q"])>1e-12:
                raise ValueError("Wyoming original q and source native taxon identity mismatch")
            ci=r["animal_bootstrap_95_percentile"]
        else:
            raise ValueError("Unknown original Utah/Wyoming fixed study")
        p=paired_old[(region,taxon)]["source_geometry_paired"]
        if p is None or p["same_original_tagged_animals_all_axes"]<5:
            raise ValueError("First two source GPS same-tag context frame changed")
        rows.append({
            "source_study_doi":x["source_doi"],"source_region":region,
            "taxon":taxon,"source_frame":"PRE_EXISTING_TWO_ORIGINAL_STUDIES",
            "source_status":"ORIGINAL_SOURCE_GPS_ECOLOGY_EVALUATED",
            "q":observed(x["night_day_q"],ci,x["q_n_animals"]),
            "daylength_slope":observed(x["daylength_slope_per_hour"],
                                       x["daylength_ci"],x["daylength_n_animals"]),
            "sun_raw":observed(x["mean_sunward_cosine"],x["sunward_ci"],x["sunward_n_animals"]),
            "sun_joint":observed(x["solar_heading_joint_component"],None,x["sunward_n_animals"],
                                 point_only=True),
            "same_tag_contexts":{"status":"PAIRED_TAG_ESTIMABLE",
                "paired_source_tagged_animals":p["same_original_tagged_animals_all_axes"]},
        })
    if len(rows)!=6:
        raise ValueError("Original six Utah/Wyoming taxa incomplete")
    for s in new["six_full_source_study_taxon_results"]:
        axes=s["within_same_original_tag_two_q_coverage_contexts"]
        rows.append({
            "source_study_doi":s["study_doi"],
            "source_region":s["region_metadata_hint"],
            "taxon":s["native_taxon"],
            "source_frame":"PROSPECTIVELY_FIXED_SIX_NEW_ORIGINAL_STUDIES",
            "source_status":s["original_GPS_outcome_status"],
            "q":normalize_new_source_axis(s["source_q_night_minus_day"]),
            "daylength_slope":normalize_new_source_axis(s["source_q_daylength_per_hour"]),
            "sun_raw":normalize_new_source_axis(s["source_solar_relative_heading_cos"]),
            "sun_joint":normalize_new_source_axis(s["source_joint_heading_sun_covariance"]),
            "same_tag_contexts":{"status":axes["status"],
                "paired_source_tagged_animals":axes["actual_paired_tagged_animal_support"]},
        })
    if len(rows)!=12 or len({r["source_study_doi"] for r in rows})!=8:
        raise ValueError("Macro cross-study source frame must include 2 old + 6 prospective studies")
    q=Counter(r["q"]["status"] for r in rows)
    day=Counter(r["daylength_slope"]["status"] for r in rows)
    raw=Counter(r["sun_raw"]["status"] for r in rows)
    joint=Counter(r["sun_joint"]["status"] for r in rows)
    paired=Counter(r["same_tag_contexts"]["status"] for r in rows)
    pair_total=sum(r["same_tag_contexts"]["paired_source_tagged_animals"] or 0 for r in rows)
    if (q["UNAVAILABLE"] != 1 or sum(v for k,v in q.items() if k!="UNAVAILABLE") != 11
            or sum(v for k,v in raw.items() if k!="UNAVAILABLE") != 10
            or paired["PAIRED_TAG_ESTIMABLE"] != 9 or pair_total != 288):
        raise ValueError("Source outcome supports differ from frozen 8 original study macro cohort")
    return {
        "version":"RSDS_8_ORIGINAL_GPS_STUDY_12_TAXON_TIME_SOLAR_ATLAS_V1",
        "date":"2026-10-26",
        "stage":"POST_NEW_SIX_SOURCE_RESULTS_EXPLORATORY_ORIGINAL_STUDY_COMPARISON",
        "original_study_denominator":8,
        "source_taxon_rows":12,
        "original_studies_with_admissible_q":len({r["source_study_doi"] for r in rows if r["q"]["status"]!="UNAVAILABLE"}),
        "source_taxa_with_q":12-q["UNAVAILABLE"],
        "source_taxa_with_solar_raw":12-raw["UNAVAILABLE"],
        "source_taxa_with_paired_original_tags":paired["PAIRED_TAG_ESTIMABLE"],
        "source_tag_context_pair_occurrences_NOT_global_biomes":pair_total,
        "source_equal_unit_not_12_independent_ecosystems":True,
        "original_source_file_blob_receipts":SOURCES,
        "all_sources_including_egret_original_metadata_HOLD":rows,
        "observed_nominal_source_animal_bootstrap_sign_counts":{
            "night_day_q":dict(q),"q_association_with_daylength_per_hour":dict(day),
            "mean_sun_relative_heading_cos":dict(raw),
            "joint_heading_sun_covariance":dict(joint),
            "same_original_tagged_context":dict(paired),
        },
        "old_48_radiation_study_external_negative_direction_test_UNCHANGED":NEGATIVE,
        "old_12_ERA5_provider_raw_SHA_HOLDS_UNCHANGED":True,
        "interpretation_limits":[
            "Even though 8 study DOIs are distinct, taxa share original field studies and ecological recruitment history; 12 study×taxon rows are not 12 independent ecosystem samples.",
            "Six new studies were chosen purposively on metadata and guild diversity, not randomly or as a preregistered ecosystem-trait hypothesis. The 200MB egret original file HOLD stays in denominator.",
            "q is nightly minus daytime log1p minimum observed straight-line GPS movement speed, not fraction or duration of latent activity.",
            "A negative q is more DAY-displacement; positive q more NIGHT-displacement only under source GPS sampling. Tracking duty cycle and speed sampling can confound guild differences.",
            "Photoperiod slopes are observational within-animal-year associations and can proxy resource, temperature, route and seasonal biology.",
            "Raw Sun-relative cos can vary entirely due to geographical route heading marginal; joint circular covariance is a distinct estimator.",
            "Old six joint heading covariance values do not have source animal-bootstrap intervals in the selected original cross-axis archive: keep POINT_ONLY, never infer zero-spanning uncertainty from new four source studies.",
            "All four new supported source joint solar cosine bootstrap intervals span zero, not proof of no solar response.",
            "Original same-tag context ranks select GPS observational support, not a repeated standardized ecological exposure or high/low irradiance.",
            "Unestimated new source GPS sun direction for fishing cats, 200MB egret HOLD and original ERA5 version holds are NOT ecological zero.",
            "Original tagged animal identities, locations, UTC event times, original restricted source bytes and per-animal original step effects are not published.",
            "No causal sunlight/radiation/solar-compass use, ecological competition, guild adaptation or world-population generality has been identified."
        ]
    }


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out",type=Path,
                     default=ROOT/"research/RSDS_8_STUDY_TIME_SOLAR_HEADING_MACRO_ATLAS_V1_20261026.json")
    args=cli.parse_args()
    result=macro()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                        encoding="utf-8")
    print("RSDS_8_STUDY_TIME_SOLAR_ATLAS",json.dumps({
        "original_studies":result["original_study_denominator"],
        "native_taxon_rows":result["source_taxon_rows"],
        "q_support":result["source_taxa_with_q"],
        "solar_heading_support":result["source_taxa_with_solar_raw"],
        "same_tag_pair_groups":result["source_taxa_with_paired_original_tags"],
        "nominal_source_animal_signs":result["observed_nominal_source_animal_bootstrap_sign_counts"],
        "old_negative":result["old_48_radiation_study_external_negative_direction_test_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
