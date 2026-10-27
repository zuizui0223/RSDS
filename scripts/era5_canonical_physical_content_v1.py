"""Canonical UTC-hourly ERA5 physical content identity (no GPS or clock effects).

A raw Open-Meteo response SHA256 is NOT a climate-data version identifier:
transient generationtime_ms, JSON whitespace, key order, or number serialization
may change even with equal physical values. This source-independent canonical
version includes the ordered UTC-hour grid, variable names/units, and exact
numeric arrays rounded to six physical-decimal places. Weather model grid
metadata is checked SEPARATELY, never dismissed if it changes.

This module emits only hashes and aggregate differences. For real animal GPS
source grids, callers must NEVER emit source coordinates, timestamps or URLs.
"""
from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timedelta

from scripts.solar_azimuth_era5_wind_v0 import (
    FIELDS, UNITS, ERA5OriginalSourceHold, validate_and_parse_weather
)

MODEL = "era5"
SCHEMA = "rsds-era5-hourly-utc-physical-v1"
GRID_FIELDS = ("latitude", "longitude", "elevation", "utc_offset_seconds")


def numeric_value(v):
    if v is None:
        return None
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ERA5OriginalSourceHold("Unexpected nonfinite or nonnumeric ERA5 response element")
    if float(v) == 0:
        return "0.000000"
    return f"{float(v):.6f}"


def physical_identity(raw, year):
    """Canonical physical data SHA, plus in-memory arrays for pair comparison.

    Exact continuous UTC calendar-hour alignment is required even though the
    legacy parser accepted a unique but not necessarily gapless same-year list.
    """
    if not isinstance(raw, (bytes, bytearray)) or len(raw) > 6_000_000:
        raise ERA5OriginalSourceHold("Original provider raw payload type or size invalid")
    data = json.loads(raw)
    verified, receipt = validate_and_parse_weather(raw, year)
    original_hourly = data["hourly"]
    labels = original_hourly["time"]
    expected = datetime(year, 1, 1)
    for label in labels:
        if label != expected.strftime("%Y-%m-%dT%H:%M"):
            raise ERA5OriginalSourceHold("Canonical provider UTC-hour continuity mismatch")
        expected += timedelta(hours=1)
    if expected != datetime(year + 1, 1, 1):
        raise ERA5OriginalSourceHold("Canonical provider UTC year incomplete")
    values = {}
    missing = {}
    for k in FIELDS:
        c = [numeric_value(v) for v in original_hourly[k]]
        if len(c) != len(labels):
            raise ERA5OriginalSourceHold("Canonical physical vector length mismatch")
        values[k] = c
        missing[k] = sum(v is None for v in c)
    payload = {
        "schema": SCHEMA, "model": MODEL, "utc_calendar_year": year,
        "units": {k: UNITS[k] for k in FIELDS},
        "time": labels, "values": values
    }
    serialized = json.dumps(payload,sort_keys=True,separators=(",",":"),
                            ensure_ascii=False,allow_nan=False).encode("utf-8")
    meta = {k: numeric_value(data.get(k)) for k in GRID_FIELDS}
    return {
        "raw_SHA256": hashlib.sha256(raw).hexdigest(),
        "physical_UTC_hourly_SHA256": hashlib.sha256(serialized).hexdigest(),
        "hour_count":len(labels),
        "missing_by_field":missing,
        "_physical_values":values,
        "_hour_labels":labels,
        "_grid_meta":meta,
        "previous_legacy_original_raw_sha_receipt":receipt["original_raw_json_sha256"],
        "provider_hourly_invalid_as_original_parser":receipt["climate_missing_or_invalid_hours"],
    }


def compare_two_same_request_responses(a, b, year):
    one, two = physical_identity(a,year), physical_identity(b,year)
    if one["_hour_labels"] != two["_hour_labels"]:
        raise ERA5OriginalSourceHold("Two original responses do not cover identical UTC hour labels")
    changed = {}
    maxabs = {}
    for k in FIELDS:
        x,y=one["_physical_values"][k],two["_physical_values"][k]
        count = sum(u != v for u,v in zip(x,y))
        vals=[abs(float(u)-float(v)) for u,v in zip(x,y)
              if u is not None and v is not None and u != v]
        changed[k]=count
        maxabs[k]=max(vals) if vals else 0.0
    grid_same = one["_grid_meta"] == two["_grid_meta"]
    raw_same = one["raw_SHA256"] == two["raw_SHA256"]
    physical_same = one["physical_UTC_hourly_SHA256"] == two["physical_UTC_hourly_SHA256"]
    return {
        "raw_response_SHA_identical":raw_same,
        "physical_hourly_numeric_SHA_identical":physical_same,
        "returned_weather_grid_metadata_identical":grid_same,
        "raw_SHA_disagreement_with_identical_physical_values_and_grid":
            not raw_same and physical_same and grid_same,
        "UTC_hours_verified":one["hour_count"],
        "first_new_raw_SHA256":one["raw_SHA256"],
        "second_new_raw_SHA256":two["raw_SHA256"],
        "first_new_hourly_physical_SHA256":one["physical_UTC_hourly_SHA256"],
        "second_new_hourly_physical_SHA256":two["physical_UTC_hourly_SHA256"],
        "hours_differing_by_field":changed,
        "maximum_absolute_numeric_change_by_field":maxabs,
        "first_missing_hourly_values_by_field":one["missing_by_field"],
        "second_missing_hourly_values_by_field":two["missing_by_field"],
        "physical_hourly_schema":SCHEMA,
        "hours_with_any_invalid_legacy_original_row_first":
            one["provider_hourly_invalid_as_original_parser"],
        "hours_with_any_invalid_legacy_original_row_second":
            two["provider_hourly_invalid_as_original_parser"],
        "does_not_compare_to_unarchived_original_v0_physical_values":True,
    }
