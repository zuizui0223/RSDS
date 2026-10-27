"""Synthetic known-truth canonical physical weather version tests (NO GPS)."""
import copy
import json
import unittest
from datetime import datetime, timedelta

from scripts.era5_canonical_physical_content_v1 import (
    compare_two_same_request_responses, physical_identity
)
from scripts.solar_azimuth_era5_wind_v0 import ERA5OriginalSourceHold


def example(year=2021):
    n=8784 if year%4==0 and (year%100!=0 or year%400==0) else 8760
    start=datetime(year,1,1)
    times=[(start+timedelta(hours=i)).strftime("%Y-%m-%dT%H:%M")
           for i in range(n)]
    return {
        "latitude":40.5,"longitude":-104.5,
        "elevation":1500.0,"utc_offset_seconds":0,
        "timezone":"GMT","generationtime_ms":0.11,
        "hourly_units":{
            "wind_speed_100m":"m/s","wind_direction_100m":"°",
            "shortwave_radiation":"W/m²","temperature_2m":"°C",
            "time":"iso8601"
        },
        "hourly":{
            "time":times,
            "wind_speed_100m":[1.0]*n,
            "wind_direction_100m":[250.0]*n,
            "shortwave_radiation":[0.0]*n,
            "temperature_2m":[8.0]*n,
        }
    }


def dump(x,sort=False,indent=None):
    return json.dumps(x,ensure_ascii=False,sort_keys=sort,
                      indent=indent,separators=None if indent else (",",":")).encode()


class ERA5CanonicalHourlyPhysicalContentTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload=example()

    def test_raw_metadata_serialization_changes_are_not_physical_changes(self):
        a=copy.deepcopy(self.payload)
        b=copy.deepcopy(a)
        b["generationtime_ms"]=9289.3
        b["timezone"]="UTC"
        c=compare_two_same_request_responses(dump(a),dump(b,sort=True,indent=2),2021)
        self.assertFalse(c["raw_response_SHA_identical"])
        self.assertTrue(c["physical_hourly_numeric_SHA_identical"])
        self.assertTrue(c["returned_weather_grid_metadata_identical"])
        self.assertTrue(c["raw_SHA_disagreement_with_identical_physical_values_and_grid"])
        self.assertTrue(all(n==0 for n in c["hours_differing_by_field"].values()))

    def test_raw_numeric_int_float_serialization_same_canonical_hash(self):
        a=copy.deepcopy(self.payload)
        b=copy.deepcopy(a)
        b["hourly"]["wind_speed_100m"]=[1]*8760
        b["hourly"]["temperature_2m"]=[8]*8760
        c=compare_two_same_request_responses(dump(a),dump(b),2021)
        self.assertFalse(c["raw_response_SHA_identical"])
        self.assertTrue(c["physical_hourly_numeric_SHA_identical"])

    def test_real_numeric_drift_is_detected(self):
        a=copy.deepcopy(self.payload)
        b=copy.deepcopy(a)
        b["hourly"]["shortwave_radiation"][200]=280.75
        b["hourly"]["shortwave_radiation"][201]=281.00
        c=compare_two_same_request_responses(dump(a),dump(b),2021)
        self.assertFalse(c["physical_hourly_numeric_SHA_identical"])
        self.assertEqual(c["hours_differing_by_field"]["shortwave_radiation"],2)
        self.assertEqual(c["maximum_absolute_numeric_change_by_field"]["shortwave_radiation"],281)

    def test_grid_metadata_change_is_separately_flagged_even_same_weather(self):
        a=copy.deepcopy(self.payload)
        b=copy.deepcopy(a)
        b["latitude"]=41.0
        c=compare_two_same_request_responses(dump(a),dump(b),2021)
        self.assertTrue(c["physical_hourly_numeric_SHA_identical"])
        self.assertFalse(c["returned_weather_grid_metadata_identical"])
        self.assertFalse(c["raw_SHA_disagreement_with_identical_physical_values_and_grid"])

    def test_hourly_missingness_changes_are_detected(self):
        a=copy.deepcopy(self.payload)
        b=copy.deepcopy(a)
        b["hourly"]["wind_speed_100m"][30]=None
        c=compare_two_same_request_responses(dump(a),dump(b),2021)
        self.assertFalse(c["physical_hourly_numeric_SHA_identical"])
        self.assertEqual(c["second_missing_hourly_values_by_field"]["wind_speed_100m"],1)

    def test_timezone_and_hourly_unit_fail_closed(self):
        a=copy.deepcopy(self.payload)
        a["utc_offset_seconds"]=3600
        with self.assertRaises(ERA5OriginalSourceHold):
            physical_identity(dump(a),2021)
        b=copy.deepcopy(self.payload)
        b["hourly_units"]["shortwave_radiation"]="kW/m²"
        with self.assertRaises(ERA5OriginalSourceHold):
            physical_identity(dump(b),2021)

    def test_hourly_gap_and_reordered_time_fail_closed(self):
        a=copy.deepcopy(self.payload)
        a["hourly"]["time"][1]=a["hourly"]["time"][2]
        with self.assertRaises(ERA5OriginalSourceHold):
            physical_identity(dump(a),2021)
        b=copy.deepcopy(self.payload)
        b["hourly"]["time"][1]="2021-01-01T02:30"
        with self.assertRaises(ERA5OriginalSourceHold):
            physical_identity(dump(b),2021)

    def test_full_leap_year_hours(self):
        a=example(2020)
        r=physical_identity(dump(a),2020)
        self.assertEqual(r["hour_count"],8784)
        self.assertEqual(len(r["physical_UTC_hourly_SHA256"]),64)


if __name__=="__main__":
    unittest.main()
