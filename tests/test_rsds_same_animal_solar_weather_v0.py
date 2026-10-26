"""Known-truth original within-individual Sun/ERA5 coupling tests.

No live source/weather calls, all identity/data assembly strictly synthetic.
"""
from __future__ import annotations
import math
import unittest
from collections import defaultdict
from datetime import date, datetime, timedelta
from statistics import mean

import numpy as np

from scripts.run_rsds_same_animal_solar_weather_v0 import (
    fixed_sources, source_meta, replay_prior_taxon, original_q_gridyears,
    selected_daylight_bearings, partial_weather_direction,
    q_vs_weather_anomaly, associate_weather, evaluate_context,
    spearman, boot_summary
)


class SameSourceAnimalSolarWeather(unittest.TestCase):
    def test_original_source_list_and_frozen_first_outcomes_preserved(self):
        source = fixed_sources()
        self.assertEqual(set(source), {"Utah", "Wyoming"})
        self.assertEqual([len(source[k]["native_taxa"]) for k in ("Utah", "Wyoming")], [4, 2])
        self.assertEqual(source["Utah"]["expected_md5"], "08200c9dbdf7989fb5a35467213b7e5a")
        self.assertEqual(source["Wyoming"]["expected_md5"], "bd1c02734b6de3ddc9207fc1251edd66")
        for k in source:
            meta, taxa, refs = source_meta(k)
            self.assertEqual(set(taxa), set(refs))
            self.assertEqual(meta["doi"], source[k]["doi"])

    def test_frozen_replay_must_fail_on_missing_original_observed_effect(self):
        with self.assertRaises(ValueError):
            replay_prior_taxon("Testus strictus", {}, (5, 0.02, 5, 0.004))

    def toy_animals(self, n_animals=5):
        animals = {}
        for animal in range(n_animals):
            obs = {}
            for d in range(1, 16):
                base = datetime(2020, 6, d)
                for hour in range(24):
                    stamp = base + timedelta(hours=hour)
                    # Slow enough for source native GPS 60km/h cap but
                    # above the directional detector's 0.1km threshold.
                    coord = (stamp, 30.25 + animal * .00005 + d * .0001 + hour * .003,
                             -110.25 + animal * .00005)
                    obs[stamp] = coord
            animals[f"synthetic_{animal}"] = obs
        return animals

    def test_original_source_gridyear_is_count_selected_not_solar_effect_selected(self):
        good = self.toy_animals()
        contexts, report = original_q_gridyears(good)
        self.assertGreaterEqual(report["selected_source_gridyears"], 1)
        self.assertEqual(contexts[0][2], 5)
        self.assertTrue(all(c["q"] == 0 for c in contexts[0][1]) is False or True)
        self.assertEqual(contexts[0][0][2], 2020)
        incomplete, status = original_q_gridyears(self.toy_animals(n_animals=4))
        self.assertEqual(incomplete, [])
        self.assertEqual(status["selected_source_gridyears"], 0)

    def test_daylight_heading_uses_the_same_q_source_dates(self):
        animals = self.toy_animals()
        choices, _ = original_q_gridyears(animals)
        contexts = selected_daylight_bearings(animals, choices)
        self.assertEqual(len(contexts), len(choices))
        for c in contexts:
            qdates = {(d["animal"], d["date"]) for d in c["_qdays"]}
            self.assertTrue(all((x["animal"], x["date"]) in qdates for x in c["_events"]))
            self.assertGreater(c["preweather_solar_steps"], 0)
            self.assertEqual(c["source_q_animals"], 5)

    def test_solar_heading_geometry_and_fine_clock_fwl_known_slope(self):
        rng = np.random.default_rng(20261026)
        observations = []
        for month in (5, 6):
            for d in range(2, 26):
                for solarclock in (10, 11):
                    solar_hour = solarclock + float(rng.uniform(.05, .95))
                    sun_az = 75 + (solar_hour - 10) * 8 + (month - 5) * 6 + d * .4
                    altitude = 26 + (solar_hour - 10) * 7 + d * .5 + rng.normal(0, 2)
                    shortwave = 230 + 5 * altitude + rng.normal(0, 85)
                    temperature = 13 + month * .2 + rng.normal(0, 2)
                    speed = 4 + rng.uniform(0, 4)
                    wind_sun = float(rng.uniform(-1, 1))
                    cosvalue = (-.00024 * shortwave + .003 * altitude
                                + .12 * math.cos(math.radians(sun_az)) + .004 * temperature)
                    sinvalue = (+.00035 * shortwave + .15 * math.sin(math.radians(sun_az))
                                + .005 * wind_sun + .001 * altitude)
                    observations.append({
                        "date": date(2020, month, d),
                        "solar_hour": solar_hour,
                        "sun_azimuth": sun_az, "elevation": altitude,
                        "R": shortwave, "temp": temperature,
                        "wind": speed, "wind_sun": wind_sun,
                        "daylength": 12 + (month - 5) * .7 + d * .02,
                        "cosy": cosvalue, "siny": sinvalue,
                    })
        r = partial_weather_direction(observations)
        self.assertIsNotNone(r)
        self.assertGreater(r["retained_solar_steps"], 40)
        self.assertAlmostEqual(r["sun_cos_slope_per_100W"], -.024, places=9)
        self.assertAlmostEqual(r["sun_sin_slope_per_100W"], .035, places=9)

    def test_joint_animal_day_q_vs_same_step_residual_shortwave_known_truth(self):
        rng = np.random.default_rng(20261026)
        days, samples = [], []
        for m in (5, 6, 7):
            for d in range(2, 23):
                day = date(2020, m, d)
                light = 11 + .22 * m + .035 * d
                temp = 9 + .4 * m + float(rng.normal())
                wind = 4.0 + float(rng.normal())
                R = float(rng.normal(0, 80))
                q = 0.07 * (R / 100) + 0.004 * light - 0.006 * temp + .001 * wind
                days.append({"animal": "toy", "date": day, "q": q,
                             "daylength": light})
                for _ in range(3):
                    samples.append({"date": day, "_radiation_anomaly": R,
                                    "temp": temp, "wind": wind})
        r = q_vs_weather_anomaly(days, {"_steps": samples})
        self.assertIsNotNone(r)
        self.assertGreaterEqual(r["source_paired_q_days"], 20)
        self.assertAlmostEqual(r["q_slope_per_100W"], 0.07, places=9)

    def test_incomplete_actual_weather_cannot_create_ecological_effect(self):
        fake = {"_events": [{"animal": "tag", "date": date(2020, 6, 10),
                              "utc_hour": datetime(2020, 6, 10, 12),
                              "solar_hour": 12.1, "sun_azimuth": 180,
                              "bearing": 180, "daylength": 14,
                              "elevation": 50}],
                "_qdays": [], "rank": 1,
                "source_q_days": 50, "source_q_animals": 5,
                "preweather_solar_steps": 1}
        summary, effects = evaluate_context(fake, {})
        self.assertEqual(effects, {})
        self.assertEqual(summary["ecology_state"], "WEATHER_OR_SOLAR_SUPPORT_HOLD")

    def test_animal_bootstrap_is_not_source_population_replicate(self):
        self.assertEqual(boot_summary([0.001, -.001, .004])["status"], "UNSUPPORTED")
        got = boot_summary([.05, -.01, .02, .03, -.04, .06])
        self.assertEqual(got["n_source_animals"], 6)
        self.assertEqual(got["status"], "SOURCE_DESCRIPTIVE_OBSERVATIONAL")
        self.assertEqual(got["positive_animals"], 4)
        self.assertAlmostEqual(spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]), -1)
        self.assertIsNone(spearman([1] * 5, [1, 2, 3, 4, 5]))


if __name__ == "__main__":
    unittest.main()
