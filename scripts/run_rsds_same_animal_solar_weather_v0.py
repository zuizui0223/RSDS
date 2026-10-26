"""Original tagged-animal solar/daypart/ERA5 ecological coupling, two fixed studies.

This is a new exploratory original-data route AFTER previous source outcomes.
The original GPS source/MD5 plus native q/daylength summary MUST replay before
any new hourly ERA5 query. No raw animal tracks, individual effects, exact
times or grid identifiers appear in the saved public result.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import random
from collections import defaultdict
from datetime import date, datetime, timedelta
from pathlib import Path
from statistics import mean, median

import numpy as np

from scripts.movebank_native_ecology_20261019 import HERE, original_files, checksum_and_download, daylength
from scripts.movebank_utah_four_taxa_solar_20261019 import (
    extract_time_only_original_source
)
from scripts.movebank_wyoming_two_ungulates_20261022 import streaming_time_selected_csv
from scripts.movebank_utah_diel_photoperiod_20261019 import one_taxon
from scripts.movebank_solar_movement_20261019 import valid_interval
from scripts.solar_azimuth_macro_v0 import directional_step
from scripts.solar_azimuth_era5_wind_v0 import (
    year_request_parameters, get_original_hourly_response,
    validate_and_parse_weather, wind_flow_to
)

CONTRACT = HERE / "research/RSDS_SAME_ANIMAL_SOLAR_WEATHER_V0_CONTRACT.json"
UTAH = HERE / "research/MOVEBANK_UTAH_FOUR_TAXA_SOLAR_V0_20261019.json"
UTAH_PRIOR = HERE / "research/MOVEBANK_UTAH_DIEL_PHOTOPERIOD_VERIFIED_RESULT_V0_20261019.json"
WYOMING = HERE / "research/MOVEBANK_WYOMING_TWO_UNGULATES_SOURCE_V0_20261022.json"
WYOMING_PRIOR = HERE / "research/wyoming_two_ungulates_20261022/wyoming_original_result.json"

BOOT = 2000
SEED = 20261026


def fixed_sources():
    spec = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if spec["version"] != "rsds-original-individual-temporal-heading-irradiance-v0":
        raise ValueError("Wrong separately frozen same-tagged-animal route")
    sources = spec["original_sources"]
    if [s["region"] for s in sources] != ["Utah", "Wyoming"]:
        raise ValueError("Original two field study sources changed")
    return {s["region"]: s for s in sources}


def source_meta(region):
    item = json.loads((UTAH if region == "Utah" else WYOMING).read_text(encoding="utf-8"))["source"]
    if region == "Utah":
        prior = json.loads(UTAH_PRIOR.read_text(encoding="utf-8"))
        references = {
            z["taxon"]: (z["eligible_diel_animals"], z["diel_mean"],
                          z["eligible_slope_animals"], z["slope_mean"])
            for z in prior["species"]
        }
        taxa = item["catalog_taxa"]
    else:
        prior = json.loads(WYOMING_PRIOR.read_text(encoding="utf-8"))
        references = {
            z["taxon"]: (
                z["day_night_speed_contrast"]["n_individuals"],
                z["day_night_speed_contrast"]["mean"],
                z["daylength_association"]["n_individuals"],
                z["daylength_association"]["mean"]
            ) for z in prior["source_taxa"]
        }
        taxa = item["published_species"]
    registered = fixed_sources()[region]
    if item["doi"] != registered["doi"] or taxa != registered["native_taxa"]:
        raise ValueError("Original DOI/taxa roster is not frozen")
    if item["original_md5"] != registered["expected_md5"] or set(references) != set(taxa):
        raise ValueError("Original source provider or taxon replay changed")
    return item, taxa, references


def replay_prior_taxon(taxon, animals, expected):
    result, _spatial = one_taxon(taxon, animals)
    q, light = result["day_night_speed_contrast"], result["daylength_association"]
    if q is None or light is None:
        raise ValueError("Original taxon lost its previously qualified original GPS endpoint")
    expected_n, expected_q, expected_light_n, expected_light = expected
    if (q["n_individuals"] != expected_n or light["n_individuals"] != expected_light_n
            or abs(q["mean"] - expected_q) > 1e-12
            or abs(light["mean"] - expected_light) > 1e-12):
        raise ValueError("Original per-taxon q/daylength GPS endpoint did not exactly replay")
    return {
        "source_original_q_animals": q["n_individuals"],
        "source_original_q": q["mean"],
        "source_original_daylength_animals": light["n_individuals"],
        "source_original_daylength_slope": light["mean"],
        "original_source_replay_success": True,
    }


def original_q_gridyears(animals):
    """Count-blind ORIGINAL animal-day gridyear candidate selection.

    Original animals used for source q need >=10 eligible q dates/year. A
    selected gridyear must have >=5 tagged animals with >=8 source q dates
    in that gridyear. No heading sign or weather is touched here.
    """
    by_grid = defaultdict(list)
    for animal, fixes in animals.items():
        data = sorted(fixes.values())
        days = defaultdict(lambda: {"day": [], "night": [], "lat": [], "lon": []})
        for a, b in zip(data, data[1:]):
            v = valid_interval(a[0], a[1], a[2], b[0], b[1], b[2])
            if v is None:
                continue
            d, part, speed, _elapsed, latitude = v
            z = days[d]
            z[part].append(speed)
            z["lat"].append(latitude)
            z["lon"].append((a[2] + b[2]) / 2)
        per_year = defaultdict(list)
        for d, z in days.items():
            if len(z["day"]) < 3 or len(z["night"]) < 3:
                continue
            lat, lon = mean(z["lat"]), mean(z["lon"])
            per_year[d.year].append((
                d, (math.floor(lat), math.floor(lon), d.year),
                math.log1p(median(z["night"])) - math.log1p(median(z["day"])),
                daylength(lat, d)
            ))
        for yr, eligible in per_year.items():
            if len(eligible) < 10:
                continue
            for day, grid, q, light in eligible:
                by_grid[grid].append({
                    "animal": animal, "date": day, "q": q, "daylength": light
                })
    candidates = []
    for grid, rows in by_grid.items():
        animals_by_day = defaultdict(int)
        for r in rows:
            animals_by_day[r["animal"]] += 1
        keep = {a for a, n in animals_by_day.items() if n >= 8}
        supporting = [r for r in rows if r["animal"] in keep]
        if len(keep) >= 5 and len(supporting) >= 40:
            candidates.append((grid, supporting, len(keep)))
    candidates.sort(key=lambda a: (-a[2], -len(a[1]), a[0]))
    return candidates[:2], {
        "all_source_original_gridyear_candidates": len(by_grid),
        "eligible_source_gridyear_candidates": len(candidates),
        "selected_source_gridyears": min(len(candidates), 2),
        "selection_blind_to_heading_and_weather": True,
    }


def selected_daylight_bearings(animals, selected):
    """Revisit only source q-qualified animal-days; never export their identity."""
    context = []
    lookups = defaultdict(list)
    for rank, (grid, rows, n_animals) in enumerate(selected, 1):
        c = {
            "_grid": grid,
            "_qdays": rows,
            "_events": [],
            "rank": rank,
            "source_q_days": len(rows),
            "source_q_animals": n_animals,
            "preweather_solar_steps": 0,
        }
        context.append(c)
        for r in rows:
            lookups[(r["animal"], r["date"])].append(c)
    for animal, fixes in animals.items():
        if not any(k[0] == animal for k in lookups):
            continue
        obs = sorted(fixes.values())
        for a, b in zip(obs, obs[1:]):
            v = valid_interval(a[0], a[1], a[2], b[0], b[1], b[2])
            if v is None or v[1] != "day":
                continue
            d = v[0]
            possible = lookups.get((animal, d), [])
            if not possible:
                continue
            step = directional_step(a, b)
            if step is None or step.solar_date != d:
                continue
            midpoint = a[0] + (b[0] - a[0]) / 2
            lat = (a[1] + b[1]) / 2
            lon = (a[2] + ((b[2] - a[2] + 180) % 360 - 180) / 2 + 180) % 360 - 180
            grid = (math.floor(lat), math.floor(lon), d.year)
            for c in possible:
                if grid != c["_grid"]:
                    continue
                c["_events"].append({
                    "animal": animal, "date": d, "utc_hour": midpoint.replace(minute=0, second=0, microsecond=0),
                    "solar_hour": step.solar_hour, "elevation": step.sun_elevation,
                    "sun_azimuth": step.sun_azimuth, "bearing": step.bearing,
                    "daylength": daylength(lat, d)
                })
                c["preweather_solar_steps"] += 1
    return context


def associate_weather(c, weather):
    events = []
    for row in c["_events"]:
        wx = weather.get(row["utc_hour"])
        if wx is None:
            continue
        speed, from_az, shortwave, temp = wx
        if not all(math.isfinite(x) for x in (speed, from_az, shortwave, temp)):
            continue
        wind_to = wind_flow_to(from_az)
        delta = math.radians(row["bearing"] - row["sun_azimuth"])
        e = dict(row)
        e.update({
            "R": shortwave, "temp": temp, "wind": speed,
            "wind_sun": math.cos(math.radians(wind_to - row["sun_azimuth"])),
            "cosy": math.cos(delta), "siny": math.sin(delta),
        })
        events.append(e)
    return events, len(events) / len(c["_events"]) if c["_events"] else 0.0


def demean_grouped(rows, field_names, blockkey):
    """Return within-animal original calendar-month × 1h solar-clock contrasts."""
    buckets = defaultdict(list)
    for i, e in enumerate(rows):
        buckets[blockkey(e)].append(i)
    keep = []
    used_keys = []
    for key, indexes in sorted(buckets.items()):
        if len(indexes) >= 6 and len({rows[k]["date"] for k in indexes}) >= 3:
            used_keys.append(key)
            keep.extend(indexes)
    if not keep:
        return None, [], []
    frame = [rows[i] for i in keep]
    a = np.array([[e[name] for name in field_names] for e in frame], dtype=float)
    bucket_idx = defaultdict(list)
    for j, e in enumerate(frame):
        bucket_idx[blockkey(e)].append(j)
    for idx in bucket_idx.values():
        a[idx] -= np.mean(a[idx], axis=0)
    return a, frame, list(bucket_idx)


def partial_weather_direction(rows):
    if (len(rows) < 40 or len({x["date"] for x in rows}) < 6
            or len({(x["date"].year, x["date"].month) for x in rows}) < 2):
        return None
    events = []
    for x in rows:
        z = dict(x)
        theta = math.radians(x["sun_azimuth"])
        z.update({
            "h": x["elevation"] / 45.0, "h2": (x["elevation"] / 45.0) ** 2,
            "solar_x": math.cos(theta), "solar_y": math.sin(theta),
            "L": x["daylength"] / 12, "T": x["temp"] / 20,
            "W": x["wind"] / 10, "WC": x["wind_sun"],
        })
        events.append(z)
    names = ["R", "cosy", "siny", "h", "h2", "solar_x", "solar_y", "L", "T", "W", "WC"]
    arr, good, blocks = demean_grouped(
        events, names, lambda e: (e["date"].year, e["date"].month, int(e["solar_hour"] // 1)))
    if arr is None or len(good) < 40 or len(blocks) < 2:
        return None
    nuisance = arr[:, 3:]
    if len(good) <= nuisance.shape[1] + 10 or np.linalg.matrix_rank(nuisance) < 3:
        return None
    beta, *_ = np.linalg.lstsq(nuisance, arr[:, :3], rcond=1e-10)
    residuals = arr[:, :3] - nuisance @ beta
    r, yc, ys = residuals.T
    if float(np.std(r)) < 20 or float(np.sum(r * r)) <= 1e-9:
        return None
    cos_beta = 100 * float(np.dot(r, yc) / np.dot(r, r))
    sin_beta = 100 * float(np.dot(r, ys) / np.dot(r, r))
    if not (math.isfinite(cos_beta) and math.isfinite(sin_beta)):
        return None
    for e, x in zip(good, r):
        e["_radiation_anomaly"] = float(x)
    return {
        "sun_cos_slope_per_100W": cos_beta,
        "sun_sin_slope_per_100W": sin_beta,
        "retained_solar_steps": len(good),
        "fine_solar_clock_blocks": len(blocks),
        "residual_shortwave_sd_Wm2": float(np.std(r)),
        "_steps": good,
    }


def q_vs_weather_anomaly(original_days, direction):
    """Same-tag/date original source q and daytime event-level ERA5 anomaly."""
    if direction is None:
        return None
    by_date = defaultdict(list)
    for s in direction["_steps"]:
        by_date[s["date"]].append(s)
    rows = []
    for d in original_days:
        steps = by_date.get(d["date"], [])
        if len(steps) < 3:
            continue
        rows.append({
            "date": d["date"], "q": d["q"], "R": mean(z["_radiation_anomaly"] for z in steps),
            "L": d["daylength"], "T": mean(z["temp"] for z in steps),
            "W": mean(z["wind"] for z in steps),
        })
    if len(rows) < 20 or len({(r["date"].year, r["date"].month) for r in rows}) < 2:
        return None
    names = ("R", "q", "L", "T", "W")
    # The original-source within-month q estimator is pre-specified as
    # demeaning days within each 1-month block with >=4 observed paired days.
    month = defaultdict(list)
    for i, row in enumerate(rows):
        month[(row["date"].year, row["date"].month)].append(i)
    indices = [j for k, arr in month.items() if len(arr) >= 4 for j in arr]
    if len(indices) < 20:
        return None
    good = [rows[i] for i in indices]
    frame = np.asarray([[row[k] for k in names] for row in good], dtype=float)
    cells = defaultdict(list)
    for i, row in enumerate(good):
        cells[(row["date"].year, row["date"].month)].append(i)
    if len(cells) < 2:
        return None
    for index in cells.values():
        frame[index] -= np.mean(frame[index], axis=0)
    controls = frame[:, 2:]
    if np.linalg.matrix_rank(controls) < 1 or len(good) < controls.shape[1] + 10:
        return None
    beta, *_ = np.linalg.lstsq(controls, frame[:, :2], rcond=1e-10)
    rr, qq = (frame[:, :2] - controls @ beta).T
    if np.std(rr) < 10 or np.std(qq) < .001 or np.dot(rr, rr) <= 1e-9:
        return None
    out = 100 * float(np.dot(rr, qq) / np.dot(rr, rr))
    if not math.isfinite(out):
        return None
    return {"q_slope_per_100W": out, "source_paired_q_days": len(good),
            "source_paired_q_months": len(cells),
            "q_shortwave_anomaly_sd_Wm2": float(np.std(rr))}


def boot_summary(values):
    v = [float(x) for x in values if x is not None and math.isfinite(x)]
    if len(v) < 5:
        return {"n_source_animals": len(v), "status": "UNSUPPORTED", "source_mean": None,
                "source_animal_bootstrap95": None,
                "positive_animals": sum(x > 0 for x in v),
                "negative_animals": sum(x < 0 for x in v)}
    rng = random.Random(SEED)
    n = len(v)
    results = sorted(mean(v[rng.randrange(n)] for _ in range(n)) for _ in range(BOOT))
    return {
        "n_source_animals": n, "status": "SOURCE_DESCRIPTIVE_OBSERVATIONAL",
        "source_mean": mean(v), "source_animal_bootstrap95": [results[49], results[1949]],
        "positive_animals": sum(x > 0 for x in v),
        "negative_animals": sum(x < 0 for x in v)
    }


def spearman(x, y):
    if len(x) != len(y) or len(x) < 5:
        return None
    def rank(a):
        order = sorted(range(len(a)), key=a.__getitem__)
        output = [0.0] * len(a)
        i = 0
        while i < len(a):
            j = i + 1
            while j < len(a) and a[order[j]] == a[order[i]]:
                j += 1
            for k in order[i:j]:
                output[k] = (i + j + 1) / 2
            i = j
        return output
    return float(np.corrcoef(rank(x), rank(y))[0, 1]) if len(set(x)) > 1 and len(set(y)) > 1 else None


def evaluate_context(context, weather):
    result = {
        "selected_gridyear_original_coverage_rank": context["rank"],
        "original_qualified_q_days": context["source_q_days"],
        "original_q_qualified_animals": context["source_q_animals"],
        "original_moving_daylight_steps": context["preweather_solar_steps"],
        "ecology_state": "SOURCE_HOLD",
    }
    events, match = associate_weather(context, weather)
    result["hourly_ERA5_movement_join_fraction"] = match
    result["hourly_weather_matched_solar_steps"] = len(events)
    if context["preweather_solar_steps"] < 60 or match < .9:
        result.update(ecology_state="WEATHER_OR_SOLAR_SUPPORT_HOLD")
        return result, {}
    q_lookup = defaultdict(list)
    for z in context["_qdays"]:
        q_lookup[z["animal"]].append(z)
    grouped = defaultdict(list)
    for z in events:
        grouped[z["animal"]].append(z)
    animal_effects = {}
    weak = 0
    for animal, arr in grouped.items():
        angle = partial_weather_direction(arr)
        if angle is None:
            weak += 1
            continue
        q = q_vs_weather_anomaly(q_lookup[animal], angle)
        animal_effects[animal] = {
            "sun_cos": angle["sun_cos_slope_per_100W"],
            "sun_sin": angle["sun_sin_slope_per_100W"],
            "q": q["q_slope_per_100W"] if q else None,
        }
    result.update({
        "ecology_state": "ORIGINAL_SAME_ANIMAL_WEATHER_EVALUATED" if animal_effects else "INSUFFICIENT_ANIMAL_FINE_CLOCK_SUPPORT",
        "individuals_without_qualified_fine_solar_clock": weak,
        "direction_animals_with_identifiable_slope": len(animal_effects),
        "joint_daypart_and_direction_animals": sum(x["q"] is not None for x in animal_effects.values()),
        "direction_sun_cos": boot_summary([x["sun_cos"] for x in animal_effects.values()]),
        "direction_sun_sin": boot_summary([x["sun_sin"] for x in animal_effects.values()]),
        "daypart_q_vs_irradiance": boot_summary([x["q"] for x in animal_effects.values()]),
    })
    return result, animal_effects


def analyze_original_taxon(taxon, animals, cache):
    choices, coverage = original_q_gridyears(animals)
    ctxs = selected_daylight_bearings(animals, choices)
    result = {"taxon": taxon, **coverage, "selected_contexts": [],
              "source_animal_identity_deduplicated_across_contexts": True,
              "source_primary_original_GPS_not_released": True}
    all_effects = defaultdict(lambda: defaultdict(list))
    for ctx in ctxs:
        key = ctx["_grid"]
        if key not in cache:
            try:
                params = year_request_parameters(*key)
                original = get_original_hourly_response(params)
                weather, receipt = validate_and_parse_weather(original, key[2])
                cache[key] = (weather, receipt, "VERIFIED_HOURLY_ERA5")
            except Exception as err:
                cache[key] = (None, None, type(err).__name__)
        weather, receipt, status = cache[key]
        if weather is None:
            result["selected_contexts"].append({
                "selected_gridyear_original_coverage_rank": ctx["rank"],
                "original_qualified_q_days": ctx["source_q_days"],
                "original_q_qualified_animals": ctx["source_q_animals"],
                "ecology_state": "ERA5_SOURCE_HOLD",
                "hold_error_type": status
            })
            continue
        report, individuals = evaluate_context(ctx, weather)
        report["ERA5_original_provider_response_sha256"] = receipt["original_raw_json_sha256"]
        report["ERA5_full_calendar_hour_count"] = receipt["source_hourly_records"]
        result["selected_contexts"].append(report)
        for name, vals in individuals.items():
            for k, v in vals.items():
                if v is not None:
                    all_effects[name][k].append(v)
    collapsed = {
        animal: {k: mean(values) for k, values in byvar.items()}
        for animal, byvar in all_effects.items()
    }
    paired = [(x["q"], x["sun_cos"], x["sun_sin"]) for x in collapsed.values()
              if all(k in x for k in ("q", "sun_cos", "sun_sin"))]
    result["source_animal_summary"] = {
        "sun_cos_per100W": boot_summary([v["sun_cos"] for v in collapsed.values() if "sun_cos" in v]),
        "sun_sin_per100W": boot_summary([v["sun_sin"] for v in collapsed.values() if "sun_sin" in v]),
        "q_daypart_per100W": boot_summary([v["q"] for v in collapsed.values() if "q" in v]),
        "paired_source_animals": len(paired),
        "within_source_animal_rho_q_vs_sun_cos": spearman([x[0] for x in paired], [x[1] for x in paired]),
        "within_source_animal_rho_q_vs_sun_sin": spearman([x[0] for x in paired], [x[2] for x in paired]),
        "status": "POST_RESULT_DESCRIPTIVE_NOT_CAUSAL"
    }
    return result


def run_source(region, out):
    src, taxa, references = source_meta(region)
    result = {"version": "rsds-original-individual-temporal-heading-irradiance-v0",
              "region": region, "original_source_doi": src["doi"],
              "status": "ORIGINAL_SOURCE_HOLD",
              "source_original_replayed_before_NEW_weather": False,
              "original_native_taxa": list(taxa), "taxon_results": [],
              "claim_limit": "Same original tags/days but different original GPS daytime/nighttime steps. ERA5 gridded observational anomalies not causal sunlight exposure. Group/bootstrap study-specific."}
    private = out / ("private_original_" + region.lower() + (".zip" if region == "Utah" else ".csv"))
    try:
        matches = [x for x in original_files(src["package_uuid"])
                   if x["bitstream_uuid"] == src["original_bitstream_uuid"]]
        if (len(matches) != 1 or matches[0]["bytes"] != src["original_bytes"]
                or matches[0]["checksum"]["value"].lower() != src["original_md5"].lower()):
            raise ValueError("Original Movebank source bitstream/MD5/size differs from frozen source")
        result["original_provider_integrity"] = checksum_and_download(matches[0], private)
        if region == "Utah":
            groups, counts, other = extract_time_only_original_source(private, taxa)
        else:
            groups, counts, other = streaming_time_selected_csv(private, taxa)
        result["source_taxonomy_unclassified_rows"] = other
        for taxon in taxa:
            replay = replay_prior_taxon(taxon, groups[taxon], references[taxon])
            result["taxon_results"].append({"taxon": taxon, **replay})
        result["source_original_replayed_before_NEW_weather"] = True
        print("RSDS_PREWEATHER_ORIGINAL_GPS_REPLAY_COMPLETE", region, len(taxa), flush=True)
        cache = {}
        for row in result["taxon_results"]:
            taxon = row["taxon"]
            try:
                row.update(analyze_original_taxon(taxon, groups[taxon], cache))
            except Exception as err:
                row.update(status="NEW_ECOLOGICAL_CONTEXT_HOLD", error_type=type(err).__name__)
            groups[taxon].clear()
            row["status"] = row.get("status", "POST_RESULT_ECOLOGICAL_SOURCE_EVALUATED")
            print("RSDS_SAME_ANIMAL_SOURCE_TAXON", json.dumps(row, sort_keys=True), flush=True)
        result["weather_requested_unique_original_gridyears"] = len(cache)
        result["weather_verified_unique_original_gridyears"] = sum(v[0] is not None for v in cache.values())
        result["status"] = "SAME_ORIGINAL_ANIMAL_WEATHER_FRAME_COMPLETE"
    except Exception as err:
        result["status"] = "ORIGINAL_SOURCE_HOLD"
        result["failure_type"] = type(err).__name__
        if not result["taxon_results"]:
            result["taxon_results"] = [{"taxon": t, "status": "ORIGINAL_SOURCE_HOLD"} for t in taxa]
    finally:
        private.unlink(missing_ok=True)
    filename = out / ("same_animal_original_" + region.lower() + ".json")
    filename.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("RSDS_ORIGINAL_SAME_ANIMAL_FIRST_RESULT", json.dumps(result, sort_keys=True), flush=True)
    return result


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--source", choices=["Utah", "Wyoming"], required=True)
    cli.add_argument("--out", type=Path, required=True)
    args = cli.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    run_source(args.source, args.out)


if __name__ == "__main__":
    main()
