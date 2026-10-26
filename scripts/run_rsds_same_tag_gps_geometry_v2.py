"""Weather-independent original tagged-animal time and Sun-direction geometry.

Published Movebank original GPS provider bitstreams are MD5-validated. All
native source q/NIGHT−DAY and q~daylength effects replay exactly BEFORE
calculating this NEW coverage-selected same-tagged-animal two-context
geometric niche comparison. NO new ERA5 request is made; original weather
version holds from the separate source route remain unchanged.

Private individual IDs, WGS84 GPS observations, grid identity and dates stay
inside the ephemeral CI runner. Publish only original source study summaries.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from scripts.movebank_native_ecology_20261019 import (
    HERE, original_files, checksum_and_download, daylength
)
from scripts.movebank_solar_movement_20261019 import valid_interval
from scripts.movebank_utah_four_taxa_solar_20261019 import extract_time_only_original_source
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.run_rsds_same_animal_solar_weather_v0 import (
    source_meta, replay_prior_taxon, selected_daylight_bearings
)
from scripts.run_rsds_same_tag_context_pair_v1 import paired_percentiles, ranked_spearman

CONTRACT = HERE / "research/RSDS_SAME_TAG_GPS_GEOMETRY_V2_CONTRACT.json"
ERA5_FIRST_HOLD = HERE / "research/RSDS_SAME_TAG_CONTEXT_PAIR_V1_FIRST_RESULT_20261026.json"
INITIAL_NEGATIVE = "NEGATIVE_DIRECTION_NOT_REPLICATED_IN_HELDOUT_SOURCE_SAMPLE"
BOOTSEED = 20261026
AXES = ("q", "sun_cos", "sun_sin", "geographic_marginal", "joint_circular",
        "sun_altitude", "apparent_solar_hour", "daylength_hours")


def check_contract():
    spec = json.loads(CONTRACT.read_text(encoding="utf-8"))
    earlier = json.loads(ERA5_FIRST_HOLD.read_text(encoding="utf-8"))
    if spec["version"] != "RSDS_SAME_TAG_GPS_GEOMETRY_V2":
        raise ValueError("Not the source-frozen weather-independent GPS route")
    if (earlier["native_study_taxon_rows"] != 6 or
            earlier["heldout_negative_radiation_nonreplication_UNCHANGED"] != INITIAL_NEGATIVE):
        raise ValueError("Original first study roster or preregistered negative result changed")
    assert len(spec["source_original_frames"]) == 2
    return {x["region"]: x for x in spec["source_original_frames"]}


def q_source_grid_contexts(animals):
    """Build ALL original q-only gridyear candidates before angle outcomes."""
    by_grid = defaultdict(list)
    for animal, hourly_fix in animals.items():
        original = sorted(hourly_fix.values())
        days = defaultdict(lambda: {"day": [], "night": [], "lat": [], "lon": []})
        for first, second in zip(original, original[1:]):
            step = valid_interval(first[0], first[1], first[2],
                                  second[0], second[1], second[2])
            if step is None:
                continue
            d, part, speed, hours, midlat = step
            z = days[d]
            z[part].append(speed)
            z["lat"].append(midlat)
            z["lon"].append((first[2] + second[2]) / 2)
        by_year = defaultdict(list)
        for d, z in days.items():
            if len(z["day"]) < 3 or len(z["night"]) < 3:
                continue
            la, lo = mean(z["lat"]), mean(z["lon"])
            by_year[d.year].append((
                d, (math.floor(la), math.floor(lo), d.year),
                math.log1p(median(z["night"])) - math.log1p(median(z["day"])),
                daylength(la, d)
            ))
        for yr, year_days in by_year.items():
            if len(year_days) < 10:
                continue
            for day, grid, q, light in year_days:
                by_grid[grid].append({
                    "animal": animal, "date": day, "q": q, "daylength": light
                })
    selected = []
    for key, rows in by_grid.items():
        count = defaultdict(int)
        for row in rows:
            count[row["animal"]] += 1
        keep = {a: n for a, n in count.items() if n >= 8}
        source_rows = [row for row in rows if row["animal"] in keep]
        if len(keep) >= 5 and len(source_rows) >= 40:
            selected.append({
                "_key": key,
                "_rows": source_rows,
                "_count": keep,
                "_n": len(keep),
                "_days": len(source_rows)
            })
    return selected, {"all_original_q_gridyears": len(by_grid),
                      "eligible_q_coverage_gridyears": len(selected)}


def select_matched_tag_contexts(candidates):
    """Choose max shared ORIGINAL tagged-animal support; no ecological signs."""
    if len(candidates) < 2:
        return [], {
            "eligible_context_pairs": 0,
            "original_q_tags_in_both_selected_contexts": 0,
            "selection_is_coverage_only": True,
        }
    combos = []
    for a, b in itertools.combinations(candidates, 2):
        overlap = set(a["_count"]) & set(b["_count"])
        matched_q = sum(min(a["_count"][tag], b["_count"][tag]) for tag in overlap)
        # Public exposure and q direction are deliberately NOT in this key.
        rank = (-len(overlap), -matched_q, -(a["_days"] + b["_days"]),
                min(a["_key"], b["_key"]), max(a["_key"], b["_key"]))
        combos.append((rank, a, b, len(overlap), matched_q))
    combos.sort(key=lambda z: z[0])
    _, a, b, overlap, support_days = combos[0]
    duo = sorted((a, b), key=lambda e: (-e["_n"], -e["_days"], e["_key"]))
    selected = [(d["_key"], d["_rows"], d["_n"]) for d in duo]
    return selected, {
        "eligible_context_pairs": len(combos),
        "original_q_tags_in_both_selected_contexts": overlap,
        "sum_of_minimum_tag_q_days_in_overlap": support_days,
        "selected_context_qualified_q_days": [d["_days"] for d in duo],
        "selected_context_qualified_q_animals": [d["_n"] for d in duo],
        "selection_is_coverage_only": True,
        "two_context_rank_means_only_q_sampling_coverage": True
    }


def individual_geometry(ctx):
    """Within each original tag use days as the equal unit for q and direction."""
    source_q = defaultdict(dict)
    for row in ctx["_qdays"]:
        key = (row["animal"], row["date"])
        if key in source_q:
            raise ValueError("Duplicate original tagged-animal/day in source context")
        source_q[row["animal"]][row["date"]] = row
    day_vectors = defaultdict(list)
    for e in ctx["_events"]:
        day_vectors[(e["animal"], e["date"])].append(e)
    result = {}
    for animal, dates in source_q.items():
        measurements = []
        clock_bins = set()
        for d, qrow in sorted(dates.items()):
            steps = day_vectors.get((animal, d), [])
            if len(steps) < 3:
                continue
            bearings = [math.radians(e["bearing"]) for e in steps]
            suns = [math.radians(e["sun_azimuth"]) for e in steps]
            cos_b = mean(math.cos(b) for b in bearings)
            sin_b = mean(math.sin(b) for b in bearings)
            cos_a = mean(math.cos(a) for a in suns)
            sin_a = mean(math.sin(a) for a in suns)
            raw_cos = mean(math.cos(b-a) for b, a in zip(bearings, suns))
            raw_sin = mean(math.sin(b-a) for b, a in zip(bearings, suns))
            measurements.append({
                "q": qrow["q"], "cos_b": cos_b, "sin_b": sin_b,
                "cos_a": cos_a, "sin_a": sin_a, "raw_cos": raw_cos,
                "raw_sin": raw_sin,
                "h": mean(e["elevation"] for e in steps),
                "hour": mean(e["solar_hour"] for e in steps),
                "L": qrow["daylength"],
                "steps": len(steps),
            })
            clock_bins.update(int(e["solar_hour"] // 1) for e in steps)
        if len(measurements) < 8 or len(clock_bins) < 2:
            continue
        def days_mean(key): return mean(x[key] for x in measurements)
        raw = days_mean("raw_cos")
        marginal = (days_mean("cos_b") * days_mean("cos_a") +
                    days_mean("sin_b") * days_mean("sin_a"))
        z = {
            "q": days_mean("q"),
            "sun_cos": raw,
            "sun_sin": days_mean("raw_sin"),
            "geographic_marginal": marginal,
            "joint_circular": raw - marginal,
            "sun_altitude": days_mean("h"),
            "apparent_solar_hour": days_mean("hour"),
            "daylength_hours": days_mean("L"),
            "_qdays": len(measurements),
            "_sun_steps": sum(x["steps"] for x in measurements),
        }
        if any(not math.isfinite(z[key]) for key in AXES):
            raise ValueError("Nonfinite original tagged geometry estimate")
        result[animal] = z
    return result


def compare_private_tags(context1, context2):
    tags = sorted(set(context1) & set(context2))
    result = {}
    diffs = defaultdict(dict)
    for i, axis in enumerate(AXES):
        delta = [context2[tag][axis] - context1[tag][axis] for tag in tags]
        for tag, diff in zip(tags, delta):
            diffs[tag][axis] = diff
        n = len(delta)
        result[axis] = {
            "common_original_tagged_animals": n,
            "status": "PAIRED_TAG_ESTIMABLE" if n >= 5 else "PAIRED_TAG_UNAVAILABLE",
            "rank2_minus_rank1_tag_mean": mean(delta) if n >= 5 else None,
            "paired_tag_bootstrap95": paired_percentiles(delta, seed=BOOTSEED + i),
            "positive_tag_deltas": sum(v > 0 for v in delta),
            "negative_tag_deltas": sum(v < 0 for v in delta),
            "sign_flip_within_same_tag": sum(
                context1[tag][axis] * context2[tag][axis] < 0 for tag in tags
            ),
        }
    return {
        "same_original_tagged_animals_all_axes": len(tags),
        "source_effects": result,
        "within_tag_rho_delta_q_vs_delta_joint_circular": ranked_spearman(
            [diffs[t]["q"] for t in tags], [diffs[t]["joint_circular"] for t in tags]
        ),
        "within_tag_rho_delta_q_vs_delta_raw_sun_cos": ranked_spearman(
            [diffs[t]["q"] for t in tags], [diffs[t]["sun_cos"] for t in tags]
        ),
        "original_private_animal_ids_exported": False
    }


def analyze_taxon(taxon, original):
    candidate, scope = q_source_grid_contexts(original)
    selected, selection = select_matched_tag_contexts(candidate)
    result = {
        "taxon": taxon,
        **scope,
        **selection,
        "original_provider_only_no_ERA5": True,
        "source_status": "PAIRED_TAG_UNAVAILABLE",
        "context_support": [],
        "paired_context_effect": None
    }
    if len(selected) != 2:
        return result
    contexts = selected_daylight_bearings(original, selected)
    per_context = []
    for c in contexts:
        tagged = individual_geometry(c)
        result["context_support"].append({
            "source_context_coverage_rank": c["rank"],
            "original_q_qualified_source_animals": c["source_q_animals"],
            "original_q_qualified_days": c["source_q_days"],
            "original_daylight_heading_steps": c["preweather_solar_steps"],
            "original_tags_with_complete_time_and_heading": len(tagged)
        })
        per_context.append(tagged)
    paired = compare_private_tags(per_context[0], per_context[1])
    result["paired_context_effect"] = paired
    result["source_status"] = (
        "PAIRED_TAG_ESTIMABLE" if paired["same_original_tagged_animals_all_axes"] >= 5
        else "PAIRED_TAG_UNAVAILABLE"
    )
    return result


def run_source(region, out):
    frame = check_contract()[region]
    src, taxa, prior = source_meta(region)
    if src["doi"] != frame["doi"] or taxa != frame["taxa"]:
        raise ValueError("Provider study/taxon identity differs from v2 frozen frame")
    result = {
        "version": "RSDS_SAME_TAG_GPS_GEOMETRY_V2",
        "date": "2026-10-26",
        "region": region,
        "source_doi": src["doi"],
        "frozen_study_native_taxa": list(taxa),
        "status": "SOURCE_HOLD",
        "original_provider_and_q_daylight_replay_before_heading": False,
        "source_taxa": [],
        "old_weather_SHA_holds_UNCHANGED": True,
        "historic_negative_shortwave_replication_UNCHANGED":
            INITIAL_NEGATIVE
    }
    filename = out / ("private_original_rsds_" + region.lower() +
                      (".zip" if region == "Utah" else ".csv"))
    try:
        candidates = [x for x in original_files(src["package_uuid"])
                      if x["bitstream_uuid"] == src["original_bitstream_uuid"]]
        if (len(candidates) != 1 or candidates[0]["bytes"] != src["original_bytes"]
                or candidates[0]["checksum"]["value"].lower() != src["original_md5"].lower()):
            raise ValueError("Original archived GPS provider bitstream integrity mismatch")
        result["original_provider_verified"] = checksum_and_download(candidates[0], filename)
        if region == "Utah":
            animals, counts, unknown = extract_time_only_original_source(filename, taxa)
        else:
            animals, counts, unknown = streaming_time_selected_csv(filename, taxa)
        result["source_nonmatching_native_taxonomy_rows"] = unknown
        for taxon in taxa:
            prior_result = replay_prior_taxon(taxon, animals[taxon], prior[taxon])
            result["source_taxa"].append({"taxon": taxon, **prior_result,
                                          "geometry_source_status": "PENDING"})
        result["original_provider_and_q_daylight_replay_before_heading"] = True
        print("RSDS_GEOMETRY_ORIGINAL_Q_DAYLIGHT_REPLAY", region, len(taxa), flush=True)
        for taxon, group in zip(taxa, result["source_taxa"]):
            try:
                group.update(analyze_taxon(taxon, animals[taxon]))
                group["geometry_source_status"] = group["source_status"]
            except Exception as exc:
                group["geometry_source_status"] = "SOURCE_HOLD"
                group["hold_error_type"] = type(exc).__name__
            animals[taxon].clear()
            print("RSDS_GEOMETRY_TAG_TAXON", json.dumps(group, sort_keys=True), flush=True)
        result["status"] = "ORIGINAL_SOURCE_GPS_GEOMETRY_COMPLETE"
    except Exception as exc:
        result["status"] = "SOURCE_HOLD"
        result["source_error_type"] = type(exc).__name__
        if not result["source_taxa"]:
            result["source_taxa"] = [{"taxon": t, "geometry_source_status":"SOURCE_HOLD"}
                                     for t in taxa]
    finally:
        filename.unlink(missing_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    (out / ("source_" + region.lower() + "_geometry_v2.json")).write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("RSDS_GPS_ONLY_ORIGINAL_TAG_CONTEXT_RESULT", json.dumps({
        "region": region,
        "status": result["status"],
        "original_source_taxa": len(result["source_taxa"]),
        "with_paired_tags": sum(
            r.get("geometry_source_status") == "PAIRED_TAG_ESTIMABLE"
            for r in result["source_taxa"])
    }, sort_keys=True), flush=True)
    return result


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--region", choices=("Utah","Wyoming"), required=True)
    cli.add_argument("--out", required=True, type=Path)
    args = cli.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    run_source(args.region, args.out)


if __name__ == "__main__":
    main()
