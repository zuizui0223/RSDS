"""Evidence-locked within-one-Wyoming-study seasonal source/taxon comparison.

Reads NO original GPS. Pin the original newly estimated source result and
both former 8-study and triple-source results by their Git blob identity.
This distinguishes observation-calendar overlap from the causal biology:
Wyoming mule deer and bighorn were coobserved on the same study-local dates,
but different individual tags, habitats, and collars remain possible.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CONTRACT="research/RSDS_WYOMING_COOCCURRENT_DATES_PHOTOPERIOD_CONTRAST_V1_CONTRACT.json"
SOURCE="research/RSDS_WYOMING_COOCCURRENT_DATES_PHOTOPERIOD_V1_FIRST_RESULT_20261029.json"
PREVIOUS="research/RSDS_TRIPLE_DATE_ECOLOGICAL_Q_DAYLENGTH_INVARIANCE_V1_FIRST_RESULT_20261028.json"
EXPECTED_SOURCE_SHA="9b04b18ee73227b4e7a1715ebfb97f610af66952"
EXPECTED_PREVIOUS_SHA="157af4cacbe4359b07f4cedbed2793f8c49972b3"
OUT=ROOT/"research/RSDS_WYOMING_SHARED_SEASON_DIVERGENT_PHOTOPERIOD_VALIDATION_V1_20261029.json"
TAXA=("Odocoileus hemionus","Ovis canadensis")
CLOCKS=("original","even","odd")
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"


def checked(root,path,digest):
    raw=(root/path).read_bytes()
    calculated=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
    if calculated!=digest:
        raise ValueError("Original evidence Git blob is not the frozen source: "+path)
    return json.loads(raw.decode("utf-8"))


def validate(root=ROOT):
    root=Path(root)
    contract=json.loads((root/CONTRACT).read_text(encoding="utf-8"))
    first=checked(root,SOURCE,EXPECTED_SOURCE_SHA)
    prior=checked(root,PREVIOUS,EXPECTED_PREVIOUS_SHA)
    if (contract["version"]!="RSDS_WYOMING_COOCCURRENT_DATES_TWO_TAXON_PHOTOPERIOD_V1"
            or first["version"]!="RSDS_WYOMING_COOCCURRENT_CALENDAR_DAY_PHOTOPERIOD_V1_ORIGINAL_RESULT"
            or first["study_Doi"]!="10.5441/001/1.703"
            or first["status"]!="STUDY_COMMON_DATE_EVALUATED"
            or first["original_provider_GPS_MD5_verified"] is not True
            or first["both_native_taxa_original_all_date_and_original_triple_source_estimates_replayed"] is not True
            or set(first["native_taxa"])!=set(TAXA)
            or first["historical_12_ERA5_raw_JSON_HOLDS_UNCHANGED"] is not True
            or first["historical_48_source_signed_radiation_external_NONREPLICATED_UNCHANGED"]!=NEGATIVE):
        raise ValueError("Source first original GPS provenance/negative result changed")
    s=first["shared_support"]
    if (s["original_study_calendar_dates_with_both_taxa_observed_on_all_three_schedules"]!=193
            or s["distinct_original_calendar_years_with_coobserved_dates"]!=2
            or s["common_daylight_standard_reference_latitude_degrees"]!=43
            or s["original_tags_with_10_days_on_both_taxa_shared_calendar"]!={TAXA[0]:108,TAXA[1]:41}):
        raise ValueError("Original two species same-calendar observation support changed")
    earlier=next(x for x in prior["all_fixed_original_studies_and_native_taxa"]
                 if x["original_study_DOI"]==first["study_Doi"])
    previous={x["taxon"]:x for x in earlier["native_source_results"]}
    current={x["taxon"]:x for x in first["source_taxon_results"]}
    if set(current)!=set(TAXA):
        raise ValueError("Missing either original Wyoming native taxon")
    rows={}
    for taxon in TAXA:
        row=current[taxon]
        count=row["original_study_coobserved_date_tagged_animals_daylength"]
        old=previous[taxon]["strict_identical_original_tag_date_support"]
        if (count>old["original_tags_with_common_photoperiod_slope_support"] or
                count!=s["original_tags_with_10_days_on_both_taxa_shared_calendar"][taxon]):
            raise ValueError("Cross-species co-observed tag count exceeds frozen source support")
        axes=row["common_date_photoperiod_q_slope_by_schedule"]
        if axes["original_tag_count"]!=count or axes["status"]!="ESTIMABLE":
            raise ValueError("Shared-season daylength source estimator unsupported")
        expected_sign="POSITIVE" if taxon==TAXA[0] else "NEGATIVE"
        for phase in CLOCKS:
            d=axes["schedules"][phase]
            if (d["nominal_sign"]!=expected_sign or
                    not isinstance(d["animal_bootstrap95"],list) or
                    len(d["animal_bootstrap95"])!=2 or
                    not d["animal_bootstrap95"][0]<=d["mean"]<=d["animal_bootstrap95"][1]):
                raise ValueError("Source weekday/calendar native taxon sign not replicated by three GPS schedules")
        rows[taxon]={"source_original_tag_n":count,
                     "common_date_reference_photoperiod_coefficients":{
                          k:axes["schedules"][k] for k in CLOCKS}}
    effects=first["same_study_between_native_taxon"]["q_daylength_slope_mule_minus_bighorn_by_schedule"]
    if effects["native_taxon_tag_counts"]!={TAXA[0]:108,TAXA[1]:41}:
        raise ValueError("Wrong native source tagged comparison denominator")
    paired=first["calendar_date_paired_source_difference"]
    if paired["coobserved_original_source_dates"]!=193 or paired["shared_year_month_blocks"]!=8:
        raise ValueError("Original same source dated block support changed")
    contrast={}
    calendar={}
    for phase in CLOCKS:
        z=effects["schedules"][phase]
        calculated=rows[TAXA[0]]["common_date_reference_photoperiod_coefficients"][phase]["mean"] - \
                   rows[TAXA[1]]["common_date_reference_photoperiod_coefficients"][phase]["mean"]
        if (z["status"]!="NOMINAL_POSITIVE" or
                abs(z["mule_minus_bighorn_slope"]-calculated)>1e-12 or
                z["descriptive_original_tag_bootstrap95"][0]<=0):
            raise ValueError("Within-study coobserved taxa slope contrast lost original positive sign")
        b=paired["slopes"][phase]
        if (b["status"]!="NOMINAL_POSITIVE" or
                b["year_month_block_descriptive95"][0]<=0 or
                not math.isfinite(b["mule_q_minus_bighorn_q_vs_shared_reference_daylength"])):
            raise ValueError("Calendar-date paired group contrast not supported in source")
        contrast[phase]=z
        calendar[phase]=b
    return {
      "version":"RSDS_WYOMING_SAME_STUDY_CALENDAR_COHORT_PHOTOPERIOD_V1_VALIDATED",
      "date":"2026-10-29",
      "stage":"POST_OUTCOME_ORIGINAL_SOURCE_CALENDAR_SUPPORT_DIAGNOSTIC_NOT_CAUSAL_SPECIES_ADAPTATION",
      "original_study_DOI":"10.5441/001/1.703",
      "original_field_study_replications":1,
      "original_native_taxa":list(TAXA),
      "source_coobserved_calendar_dates":193,
      "distinct_original_calendar_years":2,
      "distinct_year_month_blocks_for_date_resampling":8,
      "astronomical_daylength_common_reference_latitude_degrees":43.0,
      "original_taxon_source_tagged_95pct_results":rows,
      "within_study_mule_minus_bighorn_tag_bootstrap_slopes":contrast,
      "paired_same_calendar_date_year_month_bootstrap_slopes":calendar,
      "all_three_GPS_schedules_show_mule_positive_bighorn_negative_and_positive_contrast":True,
      "prior_eight_study_triple_native_sample_tag_lower_bound_preserved":True,
      "previous_twelve_ERA5_raw_source_version_HOLDS_UNCHANGED":True,
      "previous_forty_eight_source_negative_radiation_direction_replicated":False,
      "previous_forty_eight_source_negative_direction_verdict_UNCHANGED":NEGATIVE,
      "original_input_git_blobs":{SOURCE:EXPECTED_SOURCE_SHA,PREVIOUS:EXPECTED_PREVIOUS_SHA},
      "interpretation_ceiling":[
          "This is one Wyoming original study DOI, two non-identical sets of original tagged animals, not two independently sampled field studies or an externally replicated worldwide photoperiod effect.",
          "The taxa share 193 original source calendar-year/local solar dates at the study level, not necessarily every tag observed on all 193 dates. Individual-year/date support and astronomical daylight are equalized only within the declared source support rule.",
          "Astronomical daylight for this diagnostic was standardized at 43 degrees north as specified before opening Wyoming GPS. It is not actual latitude-specific photoperiod or measured surface shortwave received by each animal.",
          "Date-paired per-day taxon q means can have unequal tagged-animal compositions and date-dependent observation probability; eight original month-year clusters are few for uncertainty and cannot remove climate, habitat, resource or study instrument confounding.",
          "All 2000 tag-bootstrap and year-month block-bootstrap 95% intervals are exploratory and unadjusted for multiplicity. Study selection was informed by previously observed opposite native slopes.",
          "q is NIGHT−DAY log1p observed straight-line GPS displacement speed, not actual total activity duration or internal photoperiod memory.",
          "UTC original/even/odd GPS schedules measure different step chords and source observation opportunities; standardized source calendar does not establish an experiment.",
          "Historical 12 old raw ERA5 response SHA HOLDs and separate independent 48-source signed radiation-direction nonreplication are not reclassified.",
          "No original GPS tags, actual locations, animal-event timestamps, years/dates or individual coefficients were published."
      ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,default=OUT)
    a=p.parse_args()
    d=validate()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                     encoding="utf-8")
    print("RSDS_WYOMING_SHARED_SOURCE_SEASON_VALIDATED",json.dumps({
        "source_DOIs":d["original_field_study_replications"],
        "shared_dates":d["source_coobserved_calendar_dates"],
        "tag_counts":{k:x["source_original_tag_n"] for k,x in d["original_taxon_source_tagged_95pct_results"].items()},
        "all_three_OPPOSITE_native_daylength_signs":d["all_three_GPS_schedules_show_mule_positive_bighorn_negative_and_positive_contrast"],
        "historic_shortswave_nonreplicated":d["previous_forty_eight_source_negative_direction_verdict_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
