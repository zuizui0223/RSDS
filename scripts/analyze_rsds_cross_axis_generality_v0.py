"""RSDS cross-axis ecology: when time-use generality differs from heading generality.

Computes from frozen source-original GPS-derived summaries, NOT paper effects.
The six matched original-study/taxon time×direction records are matched by
exact DOI + native taxon. The 48-study radiation sample is kept separate.
Everything is post-outcome descriptive; never a causal sunlight or a new
independent blind replication test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = "provenance/ODSP_SOLAR_SEASON_MIGRATION_MANIFEST_20261026.json"
TIME = "research/solar_ecological_calendar_cross_source_20261025/source_level_daily_climate_daylight.csv"
HEADING = "research/GLOBAL_SOLAR_AZIMUTH_COMBINED_V2_ALL_ORIGINAL_SOURCE_TAXON_MATRIX_20261025.csv"
REPEATS = "research/GLOBAL_SOLAR_AZIMUTH_COMBINED_V2_MACRO_SUMMARY_20261025.json"
RAD_A = "research/SOLAR_AZIMUTH_24_SOURCE_MULTI_POPULATION_ERA5_V0_FIRST_ALL_RESULTS_20261025.json"
RAD_B = "research/SOLAR_RADIATION_EXTERNAL_24_SOURCE_FIRST_ALL_RESULTS_20261025.json"
RAD_DECISION = "research/SOLAR_RADIATION_TWO_24_SOURCE_NONREPLICATION_DIAGNOSTIC_FIRST_RESULT_20261025.json"
SEASON_A = "research/MOVEBANK_SOLAR_MOVEMENT_V0_RESULT_20261019.json"
SEASON_B = "research/MOVEBANK_SOLAR_MOVEMENT_EXTENSION_V0_RESULT_20261019.json"
RESULT_PATH = "research/RSDS_CROSS_AXIS_GENERALITY_FIRST_RESULT_20261026.json"
INPUTS = (TIME, HEADING, REPEATS, RAD_A, RAD_B, RAD_DECISION, SEASON_A, SEASON_B)
DECISION_FROZEN = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"

EXPECTED_COARSE_STUDY_DOIS = {
    "10.5441/001/1.711": "Utah UDWR",
    "10.5441/001/1.703": "Wyoming mule bighorn",
}
NATIVE_TAXA = ("Antilocapra americana", "Cervus elaphus",
               "Odocoileus hemionus", "Ovis canadensis")


def git_sha(data: bytes) -> str:
    return hashlib.sha1(("blob " + str(len(data)) + "\0").encode("ascii") + data).hexdigest()


def normalize_doi(doi: str) -> str:
    d = str(doi).strip().lower().removeprefix("doi:")
    if d.startswith("https://doi.org/"):
        d = d.removeprefix("https://doi.org/")
    if not d.startswith("10.5441/"):
        raise ValueError("Unknown original-study DOI namespace; do not fuzzy-match fields")
    return d


def value(x):
    if x is None or str(x).strip() == "":
        return None
    y = float(x)
    if not math.isfinite(y):
        raise ValueError("Non-finite original ecological source estimate")
    return y


def json_read(root, path):
    return json.loads((root / path).read_text(encoding="utf-8"))


def csv_read(root, path):
    with (root / path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def validate_source_assets(root):
    manifest = json_read(root, MANIFEST)
    if (manifest["source_commit"] != "d58fb759416fc458f69ef61ea79b7b1d4164f2c5"
            or manifest["text_items"] != 390):
        raise ValueError("Original scientific source domain is not the frozen RSDS migration")
    ledger = {x["target"]: x for x in manifest["files"]}
    for path in INPUTS:
        if path not in ledger:
            raise ValueError("Required original study source is absent from frozen migration: " + path)
        raw = (root / path).read_bytes()
        record = ledger[path]
        if len(raw) != record["bytes"] or git_sha(raw) != record["blob_sha"]:
            raise ValueError("Original study source was modified after migration: " + path)
    return {p: ledger[p]["blob_sha"] for p in INPUTS}


def support_sign(mean_value, low, high):
    if mean_value is None or low is None or high is None:
        return "UNAVAILABLE"
    if low > high or not low <= mean_value <= high:
        raise ValueError("Source bootstrap bounds are invalid")
    if low > 0:
        return "POSITIVE_SOURCE_ANIMAL_BOOTSTRAP"
    if high < 0:
        return "NEGATIVE_SOURCE_ANIMAL_BOOTSTRAP"
    return "UNRESOLVED_SIGN"


def choose_six_originals(time, heading):
    if len(time) != 6:
        raise ValueError("Original Utah/Wyoming time+meteorology roster must contain 6 native taxa")
    if len(heading) != 69 or sum(x["status"] == "descriptive_source_estimated" for x in heading) != 41:
        raise ValueError("Original 69 source taxon rows and 41 admissible GPS direction outcomes changed")
    mapped = defaultdict(list)
    for x in heading:
        mapped[(normalize_doi(x["doi"]), x["taxon"])].append(x)
    records = []
    seen = set()
    for t in time:
        k = (normalize_doi(t["source_doi"]), t["taxon"])
        if k in seen or k[0] not in EXPECTED_COARSE_STUDY_DOIS or t["original_source"] != EXPECTED_COARSE_STUDY_DOIS[k[0]]:
            raise ValueError("Original independent field study and taxon cannot be uniquely matched")
        seen.add(k)
        matches = mapped.get(k, [])
        if len(matches) != 1 or matches[0]["status"] != "descriptive_source_estimated":
            raise ValueError("Diel and solar-azimuth outputs did not match a single original source")
        h = matches[0]
        fields = {
            "night_day_q": value(t["source_mean_q_night_minus_day_log1p_kmh"]),
            "daylength_slope_per_hour": value(t["raw_daylight_slope_per_hour"]),
            "daylength_ci": [value(t["raw_daylight_95CI_low"]), value(t["raw_daylight_95CI_high"])],
            "daylength_given_era5_Tmax_and_30day_water_proxy": value(t["daylight_given_Tmax_30d_PminusET0_per_hour"]),
            "adjusted_ci": [value(t["conditional_daylight_95CI_low"]), value(t["conditional_daylight_95CI_high"])],
            "mean_sunward_cosine": value(h["raw_solar_relative_alignment"]),
            "sunward_ci": [value(h["raw_animal_95_lower"]), value(h["raw_animal_95_upper"])],
            "solar_geographic_marginal": value(h["geographic_marginal_sun_term"]),
            "solar_heading_joint_component": value(h["joint_direction_solar_term"]),
            "within_solar_clock_covariance": value(h["clock_conditioned_covariance"]),
            "clock_ci": [value(h["clock_95_lower"]), value(h["clock_95_upper"])],
            "q_n_animals": int(t["night_day_qualified_animals"]),
            "daylength_n_animals": int(t["raw_daylight_slope_animals"]),
            "adjusted_n_animals": int(t["daylight_ERA5_joint_slope_animals"]),
            "sunward_n_animals": int(h["n_qualified_animals"]),
        }
        if any(fields[x] is None for x in ("night_day_q", "daylength_slope_per_hour",
                                           "mean_sunward_cosine", "solar_geographic_marginal",
                                           "solar_heading_joint_component")):
            raise ValueError("Frozen same-source original time or direction effect missing")
        if abs(fields["mean_sunward_cosine"] -
               fields["solar_geographic_marginal"] -
               fields["solar_heading_joint_component"]) > 1e-10:
            raise ValueError("Solar geographic marginal plus joint covariance decomposition changed")
        fields["daylength_sign"] = support_sign(fields["daylength_slope_per_hour"], *fields["daylength_ci"])
        fields["solar_raw_sign"] = support_sign(fields["mean_sunward_cosine"], *fields["sunward_ci"])
        fields["solar_clock_sign"] = support_sign(fields["within_solar_clock_covariance"], *fields["clock_ci"])
        fields["adjusted_sign"] = support_sign(fields["daylength_given_era5_Tmax_and_30day_water_proxy"], *fields["adjusted_ci"])
        fields["taxon"] = t["taxon"]
        fields["source_doi"] = k[0]
        fields["source_region"] = "Utah" if k[0].endswith("/1.711") else "Wyoming"
        fields["original_source_alias"] = h["source_archive"]
        # Same source, not necessarily same biological individuals/intervals across axes.
        fields["cross_axis_inferential_unit"] = "SOURCE_TAXON_OBSERVATIONAL_SUMMARY"
        records.append(fields)
    if len(records) != 6 or len({r["source_doi"] for r in records}) != 2:
        raise ValueError("Source number or independent field study identity mismatch")
    if set(x["taxon"] for x in records) != set(NATIVE_TAXA):
        raise ValueError("Native original source taxa missing")
    return sorted(records, key=lambda x: (x["taxon"], x["source_region"]))


def pair_time_and_direction(records):
    by_taxon = defaultdict(list)
    for r in records:
        by_taxon[r["taxon"]].append(r)
    out = []
    for taxon, rows in sorted(by_taxon.items()):
        if len(rows) != 2:
            continue
        if {r["source_region"] for r in rows} != {"Utah", "Wyoming"}:
            raise ValueError("Duplicate source is not an independent second region")
        x = next(r for r in rows if r["source_region"] == "Utah")
        y = next(r for r in rows if r["source_region"] == "Wyoming")
        raw_delta = x["mean_sunward_cosine"] - y["mean_sunward_cosine"]
        marginal_delta = x["solar_geographic_marginal"] - y["solar_geographic_marginal"]
        joint_delta = x["solar_heading_joint_component"] - y["solar_heading_joint_component"]
        if abs(raw_delta - marginal_delta - joint_delta) > 1e-10:
            raise ValueError("Original cross-regional solar contrast does not decompose")
        out.append({
            "taxon": taxon,
            "original_studies": 2,
            "regions": ["Utah", "Wyoming"],
            "two_regions_daylength_observed_signs": [x["daylength_sign"], y["daylength_sign"]],
            "two_regions_raw_solar_direction_signs": [x["solar_raw_sign"], y["solar_raw_sign"]],
            "daylength_slope_difference_Utah_minus_Wyoming": x["daylength_slope_per_hour"] - y["daylength_slope_per_hour"],
            "raw_sunward_mean_difference_Utah_minus_Wyoming": raw_delta,
            "solar_marginal_component_difference": marginal_delta,
            "solar_joint_component_difference": joint_delta,
            "signed_raw_solar_direction_point_estimate_agrees": (x["mean_sunward_cosine"] * y["mean_sunward_cosine"] > 0),
            "both_positive_daylength_CIs": all(r["daylength_sign"] == "POSITIVE_SOURCE_ANIMAL_BOOTSTRAP" for r in rows),
            "both_raw_direction_CIs_have_opposite_signs": (
                {x["solar_raw_sign"], y["solar_raw_sign"]} ==
                {"POSITIVE_SOURCE_ANIMAL_BOOTSTRAP", "NEGATIVE_SOURCE_ANIMAL_BOOTSTRAP"}
            ),
            "night_day_q_difference_mule_minus_bighorn": None,
            "limitation": "Source/taxon contrasts cannot separate field ecosystem, latitude, migration route, individual GPS eligibility, or telemetry protocol; CIs are original-animal descriptive.",
        })
    if len(out) != 2:
        raise ValueError("Original two species in two regional source design modified")
    for region in ("Utah", "Wyoming"):
        q = {r["taxon"]: r["night_day_q"] for r in records if r["source_region"] == region}
        diff = q["Odocoileus hemionus"] - q["Ovis canadensis"]
        for x in out:
            x[f"night_day_q_mule_minus_bighorn_{region}"] = diff
    return out


def nine_taxon_repeatability(summary):
    original = summary["taxa_with_multiple_original_archives"]
    if len(original) != 9:
        raise ValueError("Previously frozen nine native taxa with two original GPS archives changed")
    comparisons = []
    for t in original:
        pair = t["sources"]
        if len(pair) != 2 or t["source_count"] != 2 or pair[0]["source"] == pair[1]["source"]:
            raise ValueError("Original within-taxon source identity/rank ambiguous")
        for axis in ("raw", "marginal", "joint", "clock"):
            if not all(math.isfinite(float(s[axis])) for s in pair):
                raise ValueError("Original circular decomposition nonfinite")
        comparisons.append({
            "taxon": t["taxon"],
            "two_distinct_archive_ids": True,
            "raw_source_signs_agree": pair[0]["raw"] * pair[1]["raw"] > 0,
            "joint_circular_covariance_source_signs_agree": pair[0]["joint"] * pair[1]["joint"] > 0,
            "solar_clock_covariance_source_signs_agree": pair[0]["clock"] * pair[1]["clock"] > 0,
        })
    return {
        "n_taxa_with_exactly_two_archived_sources": 9,
        "same_raw_cosine_sign": sum(x["raw_source_signs_agree"] for x in comparisons),
        "opposite_raw_cosine_sign": sum(not x["raw_source_signs_agree"] for x in comparisons),
        "same_joint_circular_covariance_sign": sum(x["joint_circular_covariance_source_signs_agree"] for x in comparisons),
        "opposite_joint_circular_covariance_sign": sum(not x["joint_circular_covariance_source_signs_agree"] for x in comparisons),
        "same_within_solar_clock_covariance_sign": sum(x["solar_clock_covariance_source_signs_agree"] for x in comparisons),
        "comparisons": comparisons,
        "nonclaim": "Two original source archives per taxon, possibly overlapping biological populations, are not a world-random taxonomic replication sample. Small clock effects and sign changes alone are not validated clock mechanisms.",
    }


def radiation_time_overlap(rad_a, rad_b, decision, time, season_a, season_b):
    first = rad_a["source_groups"]
    second = rad_b["all_external_first_source_results"]
    if len(first) != 24 or len(second) != 24:
        raise ValueError("Predeclared radiation source frame changed")
    if decision["all_48_previously_frozen_source_studies_included"] != 48:
        raise ValueError("Source denominator changed")
    if decision["prior_frozen_external_replication_decision_UNCHANGED"] != DECISION_FROZEN:
        raise ValueError("Previously failed radiation replication was reclassified")
    first_dois = set()
    second_dois = set()
    for s in first:
        first_dois.update(map(normalize_doi, s["source_dois"]))
    for s in second:
        second_dois.update(map(normalize_doi, s["original_source_doIs"]))
    if len(first_dois) != 24 or len(second_dois) != 24 or first_dois & second_dois:
        raise ValueError("Original radiation cohort studies overlap on DOI")
    six_dois = {normalize_doi(r["source_doi"]) for r in time}
    if len(six_dois) != 2:
        raise ValueError("Time/heading source study cardinality changed")
    seasonal_sources = season_a["sources"] + season_b["sources"]
    seasonal_source_dois = {normalize_doi(r["doi"]) for r in seasonal_sources}
    seasonal_defined_dois = {normalize_doi(r["doi"]) for r in seasonal_sources
                             if r.get("relative_effect") is not None}
    shared = (first_dois | second_dois) & seasonal_source_dois
    matched_defined = (first_dois | second_dois) & seasonal_defined_dois
    valid = sum(x["source_mean_rad_slope_per_100W"] is not None for x in
                decision["all_original_studies_anonymized_source_descriptions"])
    if valid != 20 or len(decision["all_original_studies_anonymized_source_descriptions"]) != 48:
        raise ValueError("Original shortwave slope sample/hold count changed")
    return {
        "radiation_original_source_archives": 48,
        "radiation_estimable_shortwave_slopes": valid,
        "radiation_undefined_source_slopes_NOT_ZERO": 48 - valid,
        "fixed_diel_daylength_comparator_study_DOIs": 2,
        "shared_DOIs_with_diel_daylength_studies": len(six_dois & (first_dois | second_dois)),
        "fixed_seasonal_photoperiod_source_studies": len(seasonal_source_dois),
        "shared_DOIs_with_any_seasonal_branch_archive": len(shared),
        "shared_DOIs_with_estimable_seasonal_branch_effect": len(matched_defined),
        "earlier_study_source_mean_slope": decision["study_equal_original_shortwave_slope_two_cohort_difference"]["equal_source_earlier_mean"],
        "external_study_source_mean_slope": decision["study_equal_original_shortwave_slope_two_cohort_difference"]["equal_source_external_mean"],
        "preregistered_heldout_decision_UNCHANGED": DECISION_FROZEN,
        "interpretation": "The 48-source radiation sample does not supply the same-source identifiable temporal daylength niche estimates. A joint individual-level irradiation × temporal allocation effect cannot be inferred by cross-joining unrelated source studies.",
    }


def run(root=ROOT):
    root = Path(root).resolve()
    provenance = validate_source_assets(root)
    time = csv_read(root, TIME)
    heading = csv_read(root, HEADING)
    records = choose_six_originals(time, heading)
    return {
        "version": "rsds-astronomical-environmental-temporal-direction-macro-ecology-v0",
        "date": "2026-10-26",
        "status": "POST_OUTCOME_ORIGINAL_SOURCE_OBSERVATIONAL_CROSS_AXIS_ATLAS",
        "objective": "Identify source-level temporal niche repeatability and direction niche nonrepeatability without treating study/taxon/animal GPS estimates as identical or independently sampled world biomes.",
        "source_provenance_git_blob_sha": provenance,
        "joint_time_direction_source_taxon_atlas": records,
        "independent_study_cross_regional_taxon_pairs": pair_time_and_direction(records),
        "nine_taxon_solar_heading_repeatability": nine_taxon_repeatability(json_read(root, REPEATS)),
        "sunlight_weather_temporal_observation_overlap": radiation_time_overlap(
            json_read(root, RAD_A), json_read(root, RAD_B), json_read(root, RAD_DECISION),
            time, json_read(root, SEASON_A), json_read(root, SEASON_B)),
        "claim_ceiling": [
            "This is a cross-axis source-study match of different GPS observational estimands and subsets. It does not establish within-animal temporal and direction niche statistical independence.",
            "Original q is NIGHT minus DAY log1p observed GPS step speed, NOT latent occupation time or nightly activity proportion.",
            "Sun-relative raw cos is direction relative to instant solar azimuth, not solar-compass use; solar clock-adjusted covariance is a separate estimand.",
            "Cross-regional raw direction sign changes can arise from differing source geographic heading marginals or migration corridors, not radiation selection.",
            "Both daylight and azimuth associations are noncausal; original animal bootstraps quantify conditional source uncertainty, not global biological population generality.",
            "Weather-conditional daylight slopes are post-outcome sensitivities, can use different individual/month support, and have archived ERA5 response-version replay discrepancy. Do not interpret differences as mediation.",
            "Zero common spatial 1-degree cell-month support for Wyoming two-taxon comparison: repeated source-level q rank is not sympatric competition.",
            "Radiation's external predeclared negative-direction replication failure remains intact; the atlas cannot overwrite it.",
            "No original animal identifier, exact GPS time or coordinate or precise original flight track is imported/exported in this comparison.",
            "Camera detection records require instrument uptime and observation bias models, and are not treated as true individual latent temporal niche in this GPS atlas.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / RESULT_PATH)
    args = parser.parse_args()
    result = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    m = result["independent_study_cross_regional_taxon_pairs"]
    print("RSDS ORIGINAL_SOURCE_CROSS_AXIS_ATLAS 6 study-taxon rows in two original field studies")
    print("RSDS SOURCES repeated taxa:", ", ".join(x["taxon"] for x in m))
    for x in m:
        print("TAXON", x["taxon"], "daylength_two_positive",
              x["both_positive_daylength_CIs"], "raw_solar_opposite_signed",
              x["both_raw_direction_CIs_have_opposite_signs"])
    print("SOLAR_9_TAXA", result["nine_taxon_solar_heading_repeatability"]["same_raw_cosine_sign"],
          result["nine_taxon_solar_heading_repeatability"]["same_joint_circular_covariance_sign"])
    print("RADIATION_PREDECLARED_VERDICT", result["sunlight_weather_temporal_observation_overlap"]["preregistered_heldout_decision_UNCHANGED"])
    print("RSDS CROSS_AXIS ATLAS COMPLETE")


if __name__ == "__main__":
    main()
