"""Original GPS observation-gate separation: known truth, replays, privacy."""
from __future__ import annotations

from datetime import datetime
import copy
import unittest
from unittest.mock import patch

from scripts.rsds_35_original_gps_observation_funnel_v1 import (
    profile_native,redact_native_replay
)
from scripts import run_rsds_35_original_gps_observation_funnel_v1 as runner


def fixed_paced_tags(n=6):
    """Private *synthetic* tracks: 12 independent 24h sampled days in a year."""
    records={}
    for tag in range(n):
        source={}
        for month in range(1,13):
            for hour in range(24):
                t=datetime(2025,month,11,hour,0,0)
                # Constant spatial position yields q=0 without breaking
                # original physical-segment acceptance. Date geometry true.
                source[t]=(t,43.0,0.0)
        records["PRIVATE_SYNTHETIC_TAG_"+str(tag)]=source
    return records


class ObservationFunnelKnownTruth(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sample=profile_native("Synthetic verified native",fixed_paced_tags(6))
        cls.counts=cls.sample["source_original_support_counts"]

    def test_original_independent_tag_day_and_time_support_replays(self):
        c=self.counts
        self.assertEqual(c["source_original_tagged_animals_in_hourly_parsed_track"],6)
        self.assertEqual(c["original_selected_UTC_hour_fix_count"],6*12*24)
        self.assertEqual(c["UTC_even_hour_selected_fix_count"],6*12*12)
        self.assertEqual(c["UTC_odd_hour_selected_fix_count"],6*12*12)
        self.assertEqual(c["accepted_native_tags_with_ge10_daypart_dates"],6)
        self.assertEqual(c["accepted_native_tags_with_within_year_daylength_slope"],6)
        self.assertEqual(c["accepted_native_tag_years_with_ge10_daypart_dates"],6)
        self.assertEqual(c["accepted_native_tag_years_with_ge10_dates_and_daylength_span_ge0p5h"],6)
        self.assertEqual(c["accepted_native_tag_solar_dates_with_ge3_each_daypart"],6*12)
        d=self.sample["replayed_archived_native_q_support"]
        self.assertTrue(d["q_estimable"])
        self.assertTrue(d["daylength_estimable"])
        self.assertAlmostEqual(d["source_q_mean_for_private_replay_only"],0.0,places=12)

    def test_time_then_solar_then_geometry_gates_are_nested(self):
        c=self.counts
        a=c["original_adjacent_source_fix_intervals"]
        b=c["accepted_original_0p5_to_2h_time_intervals"]
        s=c["accepted_same_solar_date_and_daypart_intervals"]
        g=c["accepted_finite_speed_max60kmh_intervals"]
        self.assertGreater(a,b)
        self.assertGreater(b,s)
        self.assertEqual(s,g)
        x=self.sample["original_observation_gate_ratios"]
        self.assertAlmostEqual(x["original_time_0p5_to_2h_over_adjacent"],b/a)
        self.assertAlmostEqual(x["same_solar_date_daypart_over_time_eligible"],s/b)
        self.assertAlmostEqual(x["total_accepted_over_adjacent"],g/a)

    def test_sufficient_GPS_movement_speed_rejection_is_distinct(self):
        a=datetime(2025,1,11,2,0)
        b=datetime(2025,1,11,3,0)
        # ~222 km in one hour; same solar night, time valid.
        native={"PRIVATE_TAG":{
            a:(a,43.0,0.0),
            b:(b,45.0,0.0)
        }}
        c=profile_native("Synthetic speedy jump",native)["source_original_support_counts"]
        self.assertEqual(c["original_adjacent_source_fix_intervals"],1)
        self.assertEqual(c["accepted_original_0p5_to_2h_time_intervals"],1)
        self.assertEqual(c["accepted_same_solar_date_and_daypart_intervals"],1)
        self.assertEqual(c["accepted_finite_speed_max60kmh_intervals"],0)
        self.assertEqual(c["accepted_native_tags_with_ge10_daypart_dates"],0)

    def test_sparse_interval_observation_loss(self):
        a=datetime(2025,1,11,1)
        b=datetime(2025,1,11,4)
        data={"PRIVATE_TAG":{a:(a,43,0),b:(b,43,0)}}
        c=profile_native("Synthetic sparse",data)["source_original_support_counts"]
        self.assertEqual(c["original_adjacent_source_fix_intervals"],1)
        self.assertEqual(c["accepted_original_0p5_to_2h_time_intervals"],0)
        self.assertEqual(c["accepted_finite_speed_max60kmh_intervals"],0)
        self.assertEqual(c["accepted_native_tags_with_ge10_daypart_dates"],0)

    def test_private_original_coordinates_and_tag_names_are_not_exported(self):
        out=redact_native_replay(self.sample)
        self.assertNotIn("source_q_mean_for_private_replay_only",str(out))
        self.assertNotIn("PRIVATE_SYNTHETIC_TAG",str(out))
        self.assertNotIn("2025-01-11",str(out))
        self.assertNotIn("latitude",str(out))
        self.assertTrue(out["no_animal_ID_date_location_GPS_row_or_movement_values_exported"])

    def test_35_original_source_identity_and_HOLDs(self):
        c,a,b,rows,w1=runner.frozen()
        self.assertEqual(len(a),16)
        self.assertEqual(len(b),19)
        self.assertEqual(len(rows),35)
        self.assertEqual(sum(s["registered_native_count"] for s in rows),50)
        self.assertEqual(sum(s["source_MD5_confirmed"] for s in rows),28)
        self.assertEqual(sum(s["complete_source_native_q"] for s in rows),10)
        self.assertEqual(c["original_GPS_source_wide_HOLD"],3)
        with patch.object(runner,"original_files",side_effect=RuntimeError("EXPOSED")) as api:
            x=runner.one_source("wave2",10)  # Homo sapiens, no source access
            y=runner.one_source("wave1",0)  # original byte-budget HOLD
            z=runner.one_source("wave2",2)  # original source-wide parsing HOLD
            api.assert_not_called()
        self.assertEqual(x["source_status"],"HUMAN_LOCATION_PRIVACY_HOLD_NO_ACCESS")
        self.assertEqual(x["original_native_funnel"],[])
        self.assertFalse(x["original_MD5_replayed"])
        self.assertEqual(y["source_status"],"ARCHIVED_ORIGINAL_SOURCE_HOLD_UNCHANGED")
        self.assertEqual(z["source_status"],"ARCHIVED_ORIGINAL_SOURCE_HOLD_UNCHANGED")

    def test_disallowed_source_index_and_frozen_historic_nulls(self):
        with self.assertRaisesRegex(ValueError,"indices"):
            runner.one_source("wave1",16)
        with self.assertRaisesRegex(ValueError,"indices"):
            runner.one_source("wave2",19)
        with self.assertRaisesRegex(ValueError,"two source waves"):
            runner.one_source("source_good_sign_only",0)
        with patch.object(runner,"NEGATIVE","REPLICATED"):
            with self.assertRaisesRegex(ValueError,"diagnostic contract changed"):
                runner.frozen()


if __name__=="__main__":
    unittest.main()
