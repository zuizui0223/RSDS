"""RSDS paired original-tag repeated context audit (POST v0 exploratory).

Replays native published Utah/Wyoming q/daylength source outcomes and then
the FIRST frozen ERA5 weather source per original context before opening any
post-result paired-animal contrasts. Context rank2-minus-rank1 is defined by
ORIGINAL q-animal support, NOT by higher radiation/temperature or effect sign.

No original tagged-animal identifiers, coordinates, per-animal effects,
precise event times, or source original GPS payloads are emitted.
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

from scripts.movebank_native_ecology_20261019 import (
    HERE, checksum_and_download, original_files
)
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.solar_azimuth_era5_wind_v0 import (
    get_original_hourly_response, validate_and_parse_weather, year_request_parameters
)
from scripts.run_rsds_same_animal_solar_weather_v0 import (
    evaluate_context, fixed_sources, original_q_gridyears, replay_prior_taxon,
    selected_daylight_bearings, source_meta
)

CONTRACT = HERE / "research/RSDS_SAME_TAG_CONTEXT_PAIR_V1_CONTRACT.json"
BASELINE = HERE / "research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_FIRST_RESULTS_20261026.json"
BOOT = 2000
SEED = 20261026


def first_sources():
    frozen = json.loads(BASELINE.read_text(encoding="utf-8"))
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract["version"] != "RSDS_SAME_TAG_CONTEXT_PAIR_V1":
        raise ValueError("Not the post-v0 frozen same-tag context-pair route")
    if frozen["fixed_original_study_count"] != 2 or frozen["total_source_native_taxon_rows"] != 6:
        raise ValueError("Source denominator changed")
    if frozen["frozen_external_radiation_nonreplication_UNCHANGED"] != \
            contract["external_nonreplication_immutable"]:
        raise ValueError("Original source nonreplication outcome was rewritten")
    by_region = {source["region"]: source for source in frozen["sources"]}
    if set(by_region) != {"Utah", "Wyoming"}:
        raise ValueError("Original source registry changed")
    return contract, by_region


def paired_percentiles(diffs, draws=BOOT, seed=SEED):
    """Resample original *paired individuals*, never interval/fix records."""
    if len(diffs) < 5:
        return None
    if any(not math.isfinite(float(z)) for z in diffs):
        raise ValueError("Original paired tag difference is nonfinite")
    rng = random.Random(seed)
    n = len(diffs)
    vals = sorted(mean(diffs[rng.randrange(n)] for _ in range(n))
                  for _ in range(draws))
    return [vals[int(.025 * draws)],
            vals[min(draws - 1, int(.975 * draws) - 1)]]


def ranked_spearman(a, b):
    if len(a) != len(b) or len(a) < 5:
        return None
    def ranks(vec):
        ranked = [0.0] * len(vec)
        order = sorted(range(len(vec)), key=vec.__getitem__)
        i = 0
        while i < len(vec):
            j = i + 1
            while j < len(vec) and vec[order[j]] == vec[order[i]]:
                j += 1
            for k in order[i:j]:
                ranked[order[k]] = (i + j + 1) / 2
            i = j
        return ranked
    aa, bb = ranks(a), ranks(b)
    ma, mb = mean(aa), mean(bb)
    num = sum((x - ma) * (y - mb) for x, y in zip(aa, bb))
    va = sum((x - ma) ** 2 for x in aa)
    vb = sum((y - mb) ** 2 for y in bb)
    return num / math.sqrt(va * vb) if va > 0 and vb > 0 else None


def paired_context_effects(first, second):
    """Private source-tag keys are consumed but never returned."""
    keys = set(first) & set(second)
    axes = ("sun_cos", "sun_sin", "q")
    stats = {}
    deltas = defaultdict(dict)
    for i, axis in enumerate(axes):
        common = sorted(tag for tag in keys
                        if first[tag].get(axis) is not None and second[tag].get(axis) is not None)
        diffs = [second[tag][axis] - first[tag][axis] for tag in common]
        if any(not math.isfinite(float(v)) for v in diffs):
            raise ValueError("Nonfinite repeated original-tag coefficient")
        for tag, value in zip(common, diffs):
            deltas[tag][axis] = value
        n = len(diffs)
        stats[axis] = {
            "paired_source_animals": n,
            "status": "PAIRED_ESTIMABLE" if n >= 5 else "INSUFFICIENT_PAIRED_ANIMALS",
            "rank2_minus_rank1_source_tag_mean": mean(diffs) if n >= 5 else None,
            "paired_animal_bootstrap95": paired_percentiles(diffs, seed=SEED + i),
            "positive_tag_differences": sum(x > 0 for x in diffs),
            "negative_tag_differences": sum(x < 0 for x in diffs),
            "zero_tag_differences": sum(x == 0 for x in diffs),
            "within_tag_sign_flip": sum(
                first[tag][axis] * second[tag][axis] < 0 for tag in common
            ),
            "interpretation": "rank2-minus-rank1, original q-support-ranked contexts; not a high-versus-low irradiance effect."
        }
    fully_paired = sorted(tag for tag, v in deltas.items() if all(a in v for a in axes))
    q = [deltas[tag]["q"] for tag in fully_paired]
    cos = [deltas[tag]["sun_cos"] for tag in fully_paired]
    sin = [deltas[tag]["sun_sin"] for tag in fully_paired]
    return {
        "identical_tag_ids_in_both_original_contexts": len(keys),
        "per_axis": stats,
        "all_three_axis_paired_original_tags": len(fully_paired),
        "within_paired_tag_rho_delta_q_vs_delta_sun_cos": ranked_spearman(q, cos),
        "within_paired_tag_rho_delta_q_vs_delta_sun_sin": ranked_spearman(q, sin),
        "release_contains_individual_keys": False,
        "causal_or_independent_world_replication": False
    }


def check_context_before_paired(context, prior):
    """Do not compare tagged effects until original context support replays."""
    original = {
        "selected_gridyear_original_coverage_rank": context["rank"],
        "original_qualified_q_days": context["source_q_days"],
        "original_q_qualified_animals": context["source_q_animals"],
        "original_moving_daylight_steps": context["preweather_solar_steps"],
    }
    if any(prior.get(k) != v for k, v in original.items()):
        raise ValueError("Frozen original q-selected source-gridyear replay failed")
    if prior.get("ecology_state") != "ORIGINAL_SAME_ANIMAL_WEATHER_EVALUATED":
        raise ValueError("First v0 context does not have a source-identifiable original endpoint")
    return original


def verify_first_context_effects(new, old):
    if new.get("ecology_state") != old.get("ecology_state"):
        raise ValueError("First contextual estimand support changed")
    mandatory = ("hourly_ERA5_movement_join_fraction",
                 "hourly_weather_matched_solar_steps",
                 "direction_animals_with_identifiable_slope",
                 "joint_daypart_and_direction_animals")
    for k in mandatory:
        a, b = new.get(k), old.get(k)
        if a is None or b is None or not math.isclose(a, b, abs_tol=1e-10, rel_tol=1e-10):
            raise ValueError("Previous contextual support could not be reproduced: " + k)
    for response in ("direction_sun_cos", "direction_sun_sin", "daypart_q_vs_irradiance"):
        a, b = new.get(response), old.get(response)
        if not isinstance(a, dict) or not isinstance(b, dict):
            raise ValueError("Missing original response for direct tag matching")
        for k in ("n_source_animals", "source_mean", "positive_animals", "negative_animals"):
            x, y = a.get(k), b.get(k)
            if (x is None) != (y is None):
                raise ValueError("Original response support no longer matches")
            if x is not None and (not math.isfinite(float(x)) or
                                  not math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-10)):
                raise ValueError("The frozen original contextual weather association did not replay")
    return True


def analyze_taxon(taxon, animals, first_taxon, cache):
    selected, counts = original_q_gridyears(animals)
    c = selected_daylight_bearings(animals, selected)
    old = first_taxon["selected_contexts"]
    if len(c) != len(old) or len(c) != 2:
        return {
            "status": "MISSING_CONTEXT",
            "source_q_qualified_contexts": len(c),
            "frozen_first_contexts": len(old),
            "effect": None
        }
    contexts = []
    per_tag = []
    for ctx, previous in zip(c, old):
        block = {
            "original_selection_rank": ctx["rank"],
            "status": "ORIGINAL_REPLAY_HOLD",
            "original_qualified_q_days": ctx["source_q_days"],
            "original_qualified_q_animals": ctx["source_q_animals"],
        }
        try:
            check_context_before_paired(ctx, previous)
            grid = ctx["_grid"]
            if grid not in cache:
                try:
                    raw = get_original_hourly_response(year_request_parameters(*grid))
                    wx, receipt = validate_and_parse_weather(raw, grid[2])
                    cache[grid] = (wx, receipt, "SOURCE_WEATHER_VERIFIED")
                except Exception as exc:
                    cache[grid] = (None, None, type(exc).__name__)
            weather, receipt, wx_status = cache[grid]
            if weather is None:
                block.update(status="ERA5_VERSION_HOLD",hold_error_type=wx_status)
                contexts.append(block)
                per_tag.append(None)
                continue
            verified_sha = receipt["original_raw_json_sha256"]
            block["ERA5_source_raw_SHA256_replayed"] = verified_sha == previous.get("ERA5_original_provider_response_sha256")
            block["ERA5_source_original_hours"] = receipt["source_hourly_records"]
            if not block["ERA5_source_raw_SHA256_replayed"]:
                block["status"] = "ERA5_VERSION_HOLD"
                contexts.append(block)
                per_tag.append(None)
                continue
            public, effects = evaluate_context(ctx, weather)
            verify_first_context_effects(public, previous)
            block["status"] = "FIRST_CONTEXT_NUMERICAL_REPLAY_PASS"
            block["source_direction_animals"] = public["direction_animals_with_identifiable_slope"]
            block["source_joint_q_heading_animals"] = public["joint_daypart_and_direction_animals"]
            contexts.append(block)
            per_tag.append(effects)
        except Exception as exc:
            block["status"] = "ORIGINAL_REPLAY_HOLD"
            block["hold_error_type"] = type(exc).__name__
            contexts.append(block)
            per_tag.append(None)
    if len(per_tag) != 2 or any(z is None for z in per_tag):
        outcome = None
        state = "PAIRED_UNAVAILABLE_CONTEXT_OR_WEATHER_HOLD"
    else:
        outcome = paired_context_effects(per_tag[0], per_tag[1])
        state = ("PAIRED_ESTIMABLE" if any(z["status"] == "PAIRED_ESTIMABLE"
                                        for z in outcome["per_axis"].values())
                 else "INSUFFICIENT_PAIRED_ANIMALS")
    return {
        "status": state,
        "source_q_context_selection": counts,
        "frozen_first_contexts_verified": all(x["status"] == "FIRST_CONTEXT_NUMERICAL_REPLAY_PASS"
                                              for x in contexts),
        "context_coverage": contexts,
        "paired_tag_effect": outcome,
        "original_first_source_effects_RETAINED_UNCHANGED": True
    }


def run_source(region, directory):
    _, baseline = first_sources()
    source, taxa, refs = source_meta(region)
    first = baseline[region]
    if first["original_source_doi"] != source["doi"] or first["original_native_taxa"] != taxa:
        raise ValueError("Original study identity changed")
    output = {
        "version": "RSDS_SAME_TAG_CONTEXT_PAIR_V1",
        "date": "2026-10-26",
        "region": region,
        "original_source_doi": source["doi"],
        "source_original_q_daylight_replay_before_ERA5": False,
        "state": "SOURCE_HOLD",
        "source_native_taxa": list(taxa),
        "source_taxa": [],
        "paired_context_comparison_is_not_high_low_radiation": True,
        "original_48_study_radiation_nonreplication_UNCHANGED":
            "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE",
    }
    temp = directory / ("private_rsds_source_" + region.lower() +
                        (".zip" if region == "Utah" else ".csv"))
    try:
        matches = [v for v in original_files(source["package_uuid"])
                   if v["bitstream_uuid"] == source["original_bitstream_uuid"]]
        if (len(matches) != 1 or matches[0]["bytes"] != source["original_bytes"]
                or matches[0]["checksum"]["value"].lower() != source["original_md5"].lower()):
            raise ValueError("Frozen provider original GPS identity not preserved")
        output["source_original_integrity"] = checksum_and_download(matches[0], temp)
        if region == "Utah":
            tracks, counts, other = extract_time_only_original_source(temp, taxa)
        else:
            tracks, counts, other = streaming_time_selected_csv(temp, taxa)
        output["original_taxonomy_unclassified_rows"] = other
        for taxon in taxa:
            original = replay_prior_taxon(taxon, tracks[taxon], refs[taxon])
            reference = next(x for x in first["taxon_results"] if x["taxon"] == taxon)
            for k in ("source_original_q", "source_original_q_animals",
                      "source_original_daylength_slope", "source_original_daylength_animals"):
                if not math.isclose(original[k], reference[k], rel_tol=1e-12, abs_tol=1e-12):
                    raise ValueError("Original source previously frozen q/daylight not reproduced")
            output["source_taxa"].append({"taxon": taxon, **original,
                                          "new_paired_sensitivity_status": "NOT_EVALUATED"})
        output["source_original_q_daylight_replay_before_ERA5"] = True
        print("RSDS_PAIRED_SOURCE_ORIGINAL_GPS_REPLAY", region, len(taxa), flush=True)
        cache = {}
        for taxon, row in zip(taxa, output["source_taxa"]):
            old = next(x for x in first["taxon_results"] if x["taxon"] == taxon)
            try:
                row.update(analyze_taxon(taxon, tracks[taxon], old, cache))
                row["new_paired_sensitivity_status"] = row["status"]
            except Exception as exc:
                row["new_paired_sensitivity_status"] = "ORIGINAL_REPLAY_HOLD"
                row["hold_error_type"] = type(exc).__name__
            tracks[taxon].clear()
            print("RSDS_WITHIN_TAG_CONTEXT", json.dumps(row, sort_keys=True), flush=True)
        output["unique_original_weather_gridyears_refetched"] = len(cache)
        output["source_weather_contexts_verified_SHA"] = sum(
            row.get("frozen_first_contexts_verified", False)
            for row in output["source_taxa"]
        )
        output["state"] = "SOURCE_PAIRED_CONTEXT_REASSESSMENT_COMPLETE"
    except Exception as exc:
        output["state"] = "SOURCE_HOLD"
        output["failure_type"] = type(exc).__name__
        if not output["source_taxa"]:
            output["source_taxa"] = [{"taxon": t, "new_paired_sensitivity_status": "SOURCE_HOLD"}
                                     for t in taxa]
    finally:
        temp.unlink(missing_ok=True)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / ("source_" + region.lower() + "_paired_context.json")).write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("RSDS_SAME_TAG_ORIGINAL_SOURCE", json.dumps({
        "region": region, "status": output["state"],
        "qualified_source_taxa": sum(row.get("new_paired_sensitivity_status") == "PAIRED_ESTIMABLE"
                                     for row in output["source_taxa"]),
        "weather_context_SHA_checks_succeeded": output.get("source_weather_contexts_verified_SHA", 0),
    }), flush=True)
    return output


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--region", choices=("Utah", "Wyoming"), required=True)
    cli.add_argument("--out", type=Path, required=True)
    args = cli.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    run_source(args.region, args.out)


if __name__ == "__main__":
    main()
