"""Audit source-study GPS UTC parity alias and phase-specific q conclusions.

Uses only independently finished original study group-only GPS result.
The source result is Git-BLOB-pinned, no fresh animal GPS or exposed tags.
This is an OBSERVATION PROCESS sensitivity, not a behavioral phase
effect, solar mechanism or a "1h -> 2h" experimental thinning in sources
already acquired on predominantly one of two alternate UTC hours.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"research/RSDS_EIGHT_STUDY_GPS_CADENCE_MATCHED_DAY_FIRST_RESULT_20261028.json"
OUTPUT=ROOT/"research/RSDS_EIGHT_STUDY_GPS_CADENCE_PHASE_ALIAS_DIAGNOSTIC_20261028.json"
EXPECTED_SHA="8b93fe361eccbdc13647e465ddfaa24ad10c7c98"
NEGATIVE="NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
PHASES=("EVEN_UTC_HOURS","ODD_UTC_HOURS")


def source():
    raw=INPUT.read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode("ascii")+raw).hexdigest()
    if sha!=EXPECTED_SHA:
        raise ValueError("First eight-original-study observed GPS cadence result changed")
    z=json.loads(raw.decode("utf-8"))
    if (z["original_field_studies"]!=8
        or z["original_source_taxon_rows"]!=12
        or z["historical_48_source_negative_radiation_replication_UNCHANGED"]!=NEGATIVE
        or z["historic_ERA5_old_12_source_raw_SHA_HOLDS_UNCHANGED"] is not True):
        raise ValueError("Original science result denominators or history changed")
    return z


def check_group(t):
    states=t.get("even_odd_cadence_effects")
    if not isinstance(states,dict):
        return None
    if set(states)!=set(PHASES):
        raise ValueError("Original even/odd source result incomplete")
    base=t["original_hourly_GPS_data_support"]
    initial=base["original_selected_1hour_source_fixes"]
    if initial <=0 or base["UTC_phase_retained_GPS_fixes"]!=initial:
        raise ValueError("Original one-record-per-hour source denominator changed")
    events=[]
    for phase in PHASES:
        p=states[phase]
        a=p["data_support"]
        n=a["UTC_phase_retained_GPS_fixes"]
        fraction=n/initial
        if not math.isclose(fraction,a["fraction_of_original_hourly_GPS_fixes_retained"],rel_tol=0,abs_tol=1e-12):
            raise ValueError("Original GPS observation UTC parity shares disagree")
        if n<0 or n>initial:
            raise ValueError("Impossible original UTC parity fix counts")
        q=p["matched_q_velocity_proxy_difference"]
        slope=p["matched_q_daylength_slope_difference_per_hour"]
        if q["n_original_tagged_animals_paired_between_cadences"]>base["tags_with_some_q_observable"]:
            raise ValueError("The original GPS paired tags exceed source qualified animals")
        if slope["n_original_tagged_animals_paired_between_cadences"]>q["n_original_tagged_animals_paired_between_cadences"]:
            raise ValueError("Original q-daylight effect without original q support")
        if (q["status"]=="INSUFFICIENT_MATCHED_TAGS") != (q["paired_tag_mean_2h_minus_same_day_1h"] is None):
            raise ValueError("Unestimated original q cannot be converted to ecological zero")
        if q["status"]=="SOURCE_TAGGED_ANIMAL_DESCRIPTIVE" and q["n_original_tagged_animals_paired_between_cadences"]<5:
            raise ValueError("Original q effect without five paired original GPS tags")
        events.append({"phase":phase,
                       "original_UTC_hour_fix_share":fraction,
                       "original_UTC_hour_selected_fixes":n,
                       "matched_original_q_animal_count":q["n_original_tagged_animals_paired_between_cadences"],
                       "matched_q_days":p["identical_original_tag_date_support"]["original_tagged_days_identical_between_1h_and_2h"],
                       "q_two_hour_minus_onehour":q["paired_tag_mean_2h_minus_same_day_1h"],
                       "q_nominal_animal_bootstrap_sign":q["nominal_interval_sign"],
                       "q_animal_bootstrap95":q["original_tag_resample95"],
                       "daylength_slope_twohour_minus_onehour":slope["paired_tag_mean_2h_minus_same_day_1h"],
                       "daylength_nominal_animal_bootstrap_sign":slope["nominal_interval_sign"],
                       "daylength_animal_bootstrap95":slope["original_tag_resample95"],
                       "source_support_original_tag_years":p["identical_original_tag_date_support"]["source_tag_years_with_ge10_matched_daypart_days"],
                       "q_source_tagged_animals_with_some_q":a["tags_with_some_q_observable"],
                      })
    if sum(x["original_UTC_hour_selected_fixes"] for x in events)!=initial:
        raise ValueError("Even plus odd does not partition exactly the original selected GPS fixes")
    even=events[0]["original_UTC_hour_fix_share"]
    odd=events[1]["original_UTC_hour_fix_share"]
    imbalance=abs(even-odd)
    classification=("PARITY_ALIASED_ORIGINAL_COLLAR_SCHEDULE" if min(even,odd)<.2
                    else "NEAR_BALANCED_UTC_FIX_PARITY" if .4<=even<=.6
                    else "INTERMEDIATE_PHASE_UNBALANCE")
    signs=[x["q_nominal_animal_bootstrap_sign"] for x in events]
    return {
        "original_q_daily_tags_before_alternate_cadence":base["tags_with_some_q_observable"],
        "total_selected_source_original_hourly_bin_fixes":initial,
        "fixes_even_plus_odd_exact_partition":True,
        "source_original_UTC_parity_imbalance_abs_difference":imbalance,
        "source_original_fix_schedule_phase_class":classification,
        "phase_effects":events,
        "BOTH_PHASES_SAME_SIGNED_NONZERO_NOMINAL_q_EFFECT":signs[0]==signs[1] and signs[0] in ("POSITIVE","NEGATIVE"),
        "OPPOSITE_SIGNED_NONZERO_NOMINAL_q_EFFECT_ACROSS_PHASE_SPECIFIC_MATCHES":set(signs)=={"POSITIVE","NEGATIVE"},
        "warning":("The two phase effects may use different animal/date intersections; "
                   "a sign reversal between separately matched parity cohorts is NOT yet "
                   "an exactly triple-matched original animal/date contrast.")
    }


def analyze():
    first=source()
    rows=[]
    stats=Counter()
    opposite=[]
    same=[]
    study_units=defaultdict(list)
    for sr in first["all_original_source_study_and_native_taxa"]:
        doi=sr["original_study_doi"]
        study_units[doi].append(sr["source_status"])
        for taxon in sr["source_native_taxon_results"]:
            r=check_group(taxon)
            status=sr["source_status"]
            stats[status]+=1
            if r is None:
                rows.append({"original_study_doi":doi,"taxon":taxon["taxon"],
                             "source_status":status,
                             "source_original_fix_schedule_phase_class":"UNAVAILABLE",
                             "supported_q_effects":None})
                continue
            stats[r["source_original_fix_schedule_phase_class"]]+=1
            if r["OPPOSITE_SIGNED_NONZERO_NOMINAL_q_EFFECT_ACROSS_PHASE_SPECIFIC_MATCHES"]:
                opposite.append({"study_doi":doi,"taxon":taxon["taxon"]})
            if r["BOTH_PHASES_SAME_SIGNED_NONZERO_NOMINAL_q_EFFECT"]:
                same.append({"study_doi":doi,"taxon":taxon["taxon"],
                             "common_source_signed_q_bias":r["phase_effects"][0]["q_nominal_animal_bootstrap_sign"]})
            rows.append({"original_study_doi":doi,"taxon":taxon["taxon"],
                         "source_status":status,**r})
    if (len(rows)!=12 or len(study_units)!=8 or
            stats["PARITY_ALIASED_ORIGINAL_COLLAR_SCHEDULE"]!=4 or
            stats["NEAR_BALANCED_UTC_FIX_PARITY"]!=7 or
            stats["PREDECLARED_FILE_METADATA_HOLD"]!=1 or
            len(opposite)!=2):
        raise ValueError("Original 12 cohort GPS UTC parity/observation coverage changed")
    source_equal={mode:[] for mode in PHASES}
    for doi in sorted(study_units):
        studyrows=[x for x in rows if x["original_study_doi"]==doi]
        for phase in PHASES:
            m=[p["q_two_hour_minus_onehour"]
               for sr in studyrows if sr.get("supported_q_effects") is not None
               for p in sr["phase_effects"] if p["phase"]==phase
               and p["q_two_hour_minus_onehour"] is not None]
            source_equal[phase].append({
                "original_study_doi":doi,
                "native_taxon_study_mean_q_measurement_bias":
                    sum(m)/len(m) if m else None,
                "native_taxa_with_original_paired_q_effect":len(m)
            })
    return {
        "version":"RSDS_ORIGINAL_GPS_UTC_PARITY_ALIASED_CADENCE_DIAGNOSTIC_V1",
        "date":"2026-10-28",
        "stage":"POST_FIRST_PARITY_OUTCOME_OBSERVATION_PROCESS_AUDIT",
        "original_source_group_first_result_git_blob_SHA":EXPECTED_SHA,
        "original_field_studies":8,
        "source_taxon_rows":12,
        "source_cohorts_with_native_UTC_phase_fix_schedule_imbalance_below_20pct":4,
        "source_cohorts_with_near_balanced_UTC_phase_fix_schedule":7,
        "source_cohorts_metadata_HOLD":1,
        "native_taxa_with_opposite_nominal_signed_q_cadence_effects_between_even_odd_UTC_phases":opposite,
        "native_taxa_with_same_nominal_nonzero_q_cadence_effects_between_even_odd":same,
        "all_original_source_taxon_phase_records":rows,
        "source_original_study_equal_phase_q_bias":source_equal,
        "main_interpretation":("The original source one-fix-per-UTC-hour selection is NOT a guaranteed "
                "one-hour GPS fix collection schedule. Four Utah source×taxon "
                "cohorts have overwhelmingly odd UTC hour bins. Both Wyoming ungulates "
                "show opposite signed nominal parity-specific measurement effects, "
                "but calendar-day/animal intersections differ between phases, "
                "so exact triple matched support is required before attributing "
                "a difference to UTC phase alone."),
        "historical_12_ERA5_source_version_HOLDS_UNCHANGED":True,
        "historic_48_source_sunlight_direction_external_NONREPLICATED_UNCHANGED":
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE",
        "claim_ceiling":[
            "Observational UTC hour parity labels are not two experimental movement treatments or measured solar clock response.",
            "The comparator '1h' refers to at most one source fix per UTC-hour bin; original accepted GPS steps could span 0.5 to 2.0h, not necessarily exactly one hour.",
            "If only 4-13% of original fix bins occur in the EVEN phase, EVEN restriction does NOT represent retaining half the original sampling opportunities.",
            "Source q proxy uses GPS displacement chord speed and can change with interval length, trip tortuosity, daypart boundary, missingness and actual original collar duty cycle.",
            "Individual calendar years and dates overlap between old 1h and EACH phase separately, NOT necessarily between even and odd phase results. Opposite effect signs cannot be interpreted as paired same-day phase reversal without triple intersection.",
            "Original 12 source-taxon units include four Utah and two Wyoming same-study taxa; count original study DOI before global scaling.",
            "The source HOLD and all source-conditioned nominal animal CIs, including unresolved values, remain unchanged. Not a climate mediation, direct animal solar response, true activity time, or population representative inference."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,default=OUTPUT)
    args=p.parse_args()
    out=analyze()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("RSDS_GPS_CADENCE_PHASE_ALIAS",json.dumps({
        "study_DOI":out["original_field_studies"],
        "source_taxa":out["source_taxon_rows"],
        "strong_source_UTC_parity_imbalance":out["source_cohorts_with_native_UTC_phase_fix_schedule_imbalance_below_20pct"],
        "near_balanced":out["source_cohorts_with_near_balanced_UTC_phase_fix_schedule"],
        "two_source_taxa_apparent_phase_flip":out["native_taxa_with_opposite_nominal_signed_q_cadence_effects_between_even_odd_UTC_phases"],
        "historic_HOLD":out["historical_12_ERA5_source_version_HOLDS_UNCHANGED"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
