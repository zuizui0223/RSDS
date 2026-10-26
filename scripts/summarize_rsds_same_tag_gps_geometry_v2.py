"""Audit new RSDS GPS-only paired-tag ecological time vs Sun-heading scale effects.

Pure standard library; uses immutable original-source GROUP summaries only.
This is post-result synthesis, not a new original GPS analysis, causal test,
multiplicity-corrected hypothesis family or world-population meta-analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = "research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_FIRST_RESULT_20261026.json"
WEATHER = "research/RSDS_SAME_TAG_CONTEXT_PAIR_V1_FIRST_RESULT_20261026.json"
FIRST = "research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
EXPECTED_BLOBS = {
    GEOMETRY: "1ec684797f0b7aae0ec441163e1bbeec48d3bd61",
    WEATHER: "60f30797fd240e9b773777b2b53de943e5d8627c",
    FIRST: "d8d84ca3583cc26ff373fd929f3f98156ca9b0cb",
}
PRESERVED_VERDICT = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
CORE_AXES = ("q", "sun_cos", "sun_sin", "geographic_marginal", "joint_circular")
ALL_AXES = CORE_AXES + ("sun_altitude", "apparent_solar_hour", "daylength_hours")
EXPECT = {
    "Utah": {"Antilocapra americana", "Cervus elaphus",
             "Odocoileus hemionus", "Ovis canadensis"},
    "Wyoming": {"Odocoileus hemionus", "Ovis canadensis"}
}


def load_and_lock(root, path):
    payload = (root / path).read_bytes()
    git_blob = hashlib.sha1(("blob " + str(len(payload)) + "\0").encode() + payload).hexdigest()
    if git_blob != EXPECTED_BLOBS[path]:
        raise ValueError("An original scientific first result was replaced or updated: " + path)
    return json.loads(payload.decode("utf-8"))


def interval_support(axis):
    n = axis["common_original_tagged_animals"]
    if axis["status"] not in ("PAIRED_TAG_ESTIMABLE", "PAIRED_TAG_UNAVAILABLE"):
        raise ValueError("Unrecognized original tag result status")
    if n < 5:
        if axis["status"] != "PAIRED_TAG_UNAVAILABLE" or axis["rank2_minus_rank1_tag_mean"] is not None:
            raise ValueError("Unidentifiable original tag change was interpreted as zero")
        return "UNAVAILABLE"
    if axis["status"] != "PAIRED_TAG_ESTIMABLE":
        raise ValueError("Supported source tag contrast was not reported")
    effect = axis["rank2_minus_rank1_tag_mean"]
    low, high = axis["paired_tag_bootstrap95"]
    if not all(isinstance(v, (float, int)) and math.isfinite(v) for v in (effect, low, high)):
        raise ValueError("Nonfinite original tag bootstrap contrast")
    if not low <= effect <= high:
        raise ValueError("Original source interval and point summary conflict")
    if low > 0:
        return "POSITIVE_UNADJUSTED_SOURCE_TAG_INTERVAL"
    if high < 0:
        return "NEGATIVE_UNADJUSTED_SOURCE_TAG_INTERVAL"
    return "INCLUDES_ZERO"


def run(root=ROOT):
    root = Path(root).resolve()
    gps = load_and_lock(root, GEOMETRY)
    old = load_and_lock(root, WEATHER)
    first = load_and_lock(root, FIRST)
    if (gps["original_source_count"] != 2 or gps["source_taxon_rows"] != 6
            or gps["source_taxa_with_paired_tag_estimates"] != 6
            or not gps["no_ERA5_retrieved_or_reinterpreted"]
            or not gps["previous_raw_ERA5_SHA_HOLD_unchanged"]):
        raise ValueError("Original GPS-only six-taxon paired design was modified")
    if any(z != PRESERVED_VERDICT for z in (
        gps["external_48_original_source_shortwave_nonreplication_UNCHANGED"],
        old["heldout_negative_radiation_nonreplication_UNCHANGED"],
        first["frozen_external_radiation_nonreplication_UNCHANGED"]
    )):
        raise ValueError("External frozen negative directional replication changed")
    if (old["native_study_taxon_rows"] != 6
            or old["number_of_source_taxa_with_estimate_by_axis"] !=
            {"q": 0, "sun_cos": 0, "sun_sin": 0}):
        raise ValueError("Earlier strict ERA5 response holds were discarded")
    weather_contexts = [c for s in old["all_native_source_taxon_results"]
                        for c in s["original_contexts"]]
    if len(weather_contexts) != 12 or any(
            c.get("status") != "ERA5_VERSION_HOLD" or
            c.get("ERA5_source_raw_SHA256_replayed") is not False
            for c in weather_contexts):
        raise ValueError("12 source weather-version holds are not preserved")
    by_region = {k: set() for k in EXPECT}
    rows = []
    axis_status = {axis: Counter() for axis in ALL_AXES}
    source_tag_total = 0
    for record in gps["all_original_source_taxon_results"]:
        region, taxon = record["region"], record["taxon"]
        if region not in by_region or taxon in by_region[region] or taxon not in EXPECT[region]:
            raise ValueError("The fixed published original study taxonomy/denominator changed")
        by_region[region].add(taxon)
        if (record["new_context_status"] != "PAIRED_TAG_ESTIMABLE"
                or record["source_original_q_daylight_exact_replay"] is not True
                or record["all_contexts_selected_by_original_q_coverage"] is not True):
            raise ValueError("GPS original-source or paired-tag qualification changed")
        tagged = record["source_geometry_paired"]
        n = tagged["same_original_tagged_animals_all_axes"]
        if n < 5 or set(tagged["source_effects"]) != set(ALL_AXES):
            raise ValueError("Full original tag axes not present or inadequate sample")
        source_tag_total += n
        effects = {}
        for axis in ALL_AXES:
            z = tagged["source_effects"][axis]
            if z["common_original_tagged_animals"] != n:
                raise ValueError("GPS time and direction used mismatched tagged animals")
            support = interval_support(z)
            axis_status[axis][support] += 1
            effects[axis] = {
                "paired_tag_count": n,
                "context2_minus_context1_mean": z["rank2_minus_rank1_tag_mean"],
                "descriptive_tag_bootstrap95": z["paired_tag_bootstrap95"],
                "interval_classification": support,
                "signed_tag_changes": {
                    "positive": z["positive_tag_deltas"],
                    "negative": z["negative_tag_deltas"],
                },
                "within_same_tag_sign_flips": z["sign_flip_within_same_tag"],
            }
            if (z["positive_tag_deltas"] + z["negative_tag_deltas"] > n or
                    z["sign_flip_within_same_tag"] > n):
                raise ValueError("Source tagged animal sign-change tally impossible")
        diff = (effects["sun_cos"]["context2_minus_context1_mean"]
                - effects["geographic_marginal"]["context2_minus_context1_mean"]
                - effects["joint_circular"]["context2_minus_context1_mean"])
        if abs(diff) > 1e-10:
            raise ValueError("Original Sun-relative cosine vs marginal+joint identity failed")
        pair_rho = tagged["within_tag_rho_delta_q_vs_delta_joint_circular"]
        if pair_rho is not None and not (-1 <= pair_rho <= 1):
            raise ValueError("Invalid same-tag time/direction rank correlation")
        rows.append({
            "region": region,
            "taxon": taxon,
            "original_paired_tagged_animals": n,
            "q_coverage_selected_original_contexts": 2,
            "results": effects,
            "same_tag_rank_rho_delta_q_vs_delta_joint_solar_circular": pair_rho,
            "same_tag_rank_rho_delta_q_vs_delta_raw_sun_cos":
                tagged["within_tag_rho_delta_q_vs_delta_raw_sun_cos"],
        })
    if by_region != EXPECT or len(rows) != 6 or source_tag_total != 269:
        raise ValueError("Preselected original source frame/paired tag denominators changed")

    repeated = []
    for taxon in ("Odocoileus hemionus", "Ovis canadensis"):
        a, b = [next(r for r in rows if r["region"] == reg and r["taxon"] == taxon)
                for reg in ("Utah", "Wyoming")]
        repeated.append({
            "taxon": taxon,
            "source_study_count": 2,
            "q_coverage_rank2_minus_rank1_in_each_original_study": {
                "Utah": a["results"]["q"]["context2_minus_context1_mean"],
                "Wyoming": b["results"]["q"]["context2_minus_context1_mean"],
                "Utah_interval_classification": a["results"]["q"]["interval_classification"],
                "Wyoming_interval_classification": b["results"]["q"]["interval_classification"],
            },
            "no_common_high_low_environmental_exposure_assigned": True,
        })

    # These are descriptive counts of one-at-a-time, source-tag bootstraps.
    return {
        "version": "RSDS_ORIGINAL_GPS_TIME_AND_SOLAR_HEADING_V2_POSTHOC_INTERPRETATION",
        "date": "2026-10-26",
        "status": "SOURCE_ORIGINAL_PAIR_SUPPORT_CONFIRMED_NONCAUSAL",
        "study_units": 2,
        "study_taxon_rows": len(rows),
        "summed_source_paired_tag_counts_NOT_independent_biomes": source_tag_total,
        "original_source_data_git_blobs_checked": EXPECTED_BLOBS,
        "each_axis_nominal_tag_bootstrap95_classification_counts": {
            k: dict(v) for k, v in axis_status.items()
        },
        "source_taxon_ecological_results": rows,
        "repeated_native_taxa_with_two_original_studies": repeated,
        "earlier_weather_version_holds_preserved": 12,
        "earlier_weather_response_comparison_supported_taxa": 0,
        "old_external_shortwave_nonreplication_UNCHANGED": PRESERVED_VERDICT,
        "interpretation_ceiling": [
            "New paired original GPS ecological responses involve exactly the same tagged animal in two coverage-selected field-site/year contexts, not experimental adaptation or photic treatment.",
            "Day/night q is log1p observed nighttime GPS displacement minus daytime GPS displacement, not total nocturnal activity share, duration or time budget.",
            "An original animal with q-qualified nights is paired with SUN-QUALIFIED DAYTIME directional steps on eligible days; individual-night vs daylight heading responses need not be the same GPS steps.",
            "Original q on qualified day covers all valid source-day intervals, possibly outside the single one-degree geographic grid used to qualify daytime direction segments; site specificity is not perfect.",
            "Raw cosine Sun-relative bearing has an exact geographic heading marginal plus joint circular component decomposition; a nonsignificant joint coefficient DOES NOT mean the geographic marginal caused all raw variation.",
            "All nominal original tag-resampled 95 percent intervals cover separate source/taxon/axes without multiple comparison adjustment or source/biome population sampling confidence.",
            "The two source-context ranks reflect original GPS q-day and shared tag support, not chronological order, temperature, radiation exposure level, or randomized ecological treatment.",
            "The old 12 ERA5 source JSON SHA replay holds are unchanged; no new weather measurement has entered this GPS-only analysis.",
            "The original separately preregistered 48-study solar-radiation negative-direction replication failure remains untouched.",
            "No precise original animal identities, raw GPS coordinates, exact dates/times or individual measurements have been published."
        ]
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out", type=Path, default=ROOT / "research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_ECOLOGY_SYNTHESIS_20261026.json")
    args = cli.parse_args()
    d = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(d, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    print("RSDS_SOURCE_PAIR_SUPPORT", json.dumps({
        "studies": d["study_units"], "taxa": d["study_taxon_rows"],
        "tag_group_counts": d["summed_source_paired_tag_counts_NOT_independent_biomes"],
        "bootstrap_supported_axis_summary":
            d["each_axis_nominal_tag_bootstrap95_classification_counts"],
        "ERA5_HOLD": d["earlier_weather_version_holds_preserved"]
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
